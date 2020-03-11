package org.broadinstitute.hellbender.tools.copynumber;

import com.google.common.collect.ImmutableSet;
import org.apache.commons.collections.ListUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.tools.copynumber.arguments.*;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.*;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.*;
import org.broadinstitute.hellbender.tools.copynumber.models.AlleleFractionModeller;
import org.broadinstitute.hellbender.tools.copynumber.models.CopyRatioModeller;
import org.broadinstitute.hellbender.tools.copynumber.segmentation.MultidimensionalKernelSegmenter;
import org.broadinstitute.hellbender.tools.copynumber.segmentation.MultisampleMultidimensionalKernelSegmenter;
import org.broadinstitute.hellbender.tools.copynumber.utils.genotyping.NaiveHeterozygousPileupGenotypingUtils;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenter;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Models segmented copy ratios from denoised read counts and segmented minor-allele fractions from allelic counts.
 *
 * <p>
 *     Possible inputs are: 1) denoised copy ratios for the case sample, 2) allelic counts for the case sample,
 *     and 3) allelic counts for a matched-normal sample.  All available inputs will be used to to perform
 *     segmentation and model inference.
 * </p>
 *
 * <p>
 *     If allelic counts are available, the first step in the inference process is to genotype heterozygous sites,
 *     as the allelic counts at these sites will subsequently be modeled to infer segmented minor-allele fraction.
 *     We perform a relatively simple and naive genotyping based on the allele counts (i.e., pileups), which is
 *     controlled by a small number of parameters ({@code minimum-total-allele-count},
 *     {@code genotyping-homozygous-log-ratio-threshold}, and {@code genotyping-homozygous-log-ratio-threshold}).
 *     If the matched normal is available, its allelic counts will be used to genotype the sites, and
 *     we will simply assume these genotypes are the same in the case sample.  (This can be critical, for example,
 *     for determining sites with loss of heterozygosity in high purity case samples; such sites will be genotyped as
 *     homozygous if the matched-normal sample is not available.)
 * </p>
 *
 * <p>
 *     Next, we segment, if available, the denoised copy ratios and the alternate-allele fractions at the
 *     genotyped heterozygous sites.  This is done using kernel segmentation (see {@link KernelSegmenter}).
 *     Various segmentation parameters control the sensitivity of the segmentation and should be selected
 *     appropriately for each analysis.
 * </p>
 *
 * <p>
 *     If both copy ratios and allele fractions are available, we perform segmentation using a combined kernel
 *     that is sensitive to changes that occur not only in either of the two but also in both.  However, in this case,
 *     we simply discard all allele fractions at sites that lie outside of the available copy-ratio intervals
 *     (rather than imputing the missing copy-ratio data); these sites are filtered out during the genotyping step
 *     discussed above.  This can have implications for analyses involving the sex chromosomes;
 *     see comments in {@link CreateReadCountPanelOfNormals}.
 * </p>
 *
 * <p>
 *     After segmentation is complete, we run Markov-chain Monte Carlo (MCMC) to determine posteriors for
 *     segmented models for the log2 copy ratio and the minor-allele fraction; see {@link CopyRatioModeller}
 *     and {@link AlleleFractionModeller}, respectively.  After the first run of MCMC is complete,
 *     smoothing of the segmented posteriors is performed by merging adjacent segments whose posterior
 *     credible intervals sufficiently overlap according to specified segmentation-smoothing parameters.
 *     Then, additional rounds of segmentation smoothing (with intermediate MCMC optionally performed in between rounds)
 *     are performed until convergence, at which point a final round of MCMC is performed.
 * </p>
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         (Optional) Denoised-copy-ratios file from {@link DenoiseReadCounts}.
 *         If allelic counts are not provided, then this is required.
 *     </li>
 *     <li>
 *         (Optional) Allelic-counts file from {@link CollectAllelicCounts}.
 *         If denoised copy ratios are not provided, then this is required.
 *     </li>
 *     <li>
 *         (Optional) Matched-normal allelic-counts file from {@link CollectAllelicCounts}.
 *         This can only be provided if allelic counts for the case sample are also provided.
 *     </li>
 *     <li>
 *         Output prefix.
 *         This is used as the basename for output files.
 *     </li>
 *     <li>
 *         Output directory.
 *         This will be created if it does not exist.
 *     </li>
 * </ul>
 *
 * <h3>Outputs</h3>
 *
 * <ul>
 *     <li>
 *         Modeled-segments .modelBegin.seg and .modelFinal.seg files.
 *         These are tab-separated values (TSV) files with a SAM-style header containing a read group sample name, a sequence dictionary,
 *         a row specifying the column headers contained in {@link ModeledSegmentCollection.ModeledSegmentTableColumn},
 *         and the corresponding entry rows.
 *         The initial result before segmentation smoothing is output to the .modelBegin.seg file
 *         and the final result after segmentation smoothing is output to the .modelFinal.seg file.
 *     </li>
 *     <li>
 *         Allele-fraction-model global-parameter files (.modelBegin.af.param and .modelFinal.af.param).
 *         These are tab-separated values (TSV) files with a SAM-style header containing a read group sample name,
 *         a row specifying the column headers contained in {@link ParameterDecileCollection.ParameterTableColumn},
 *         and the corresponding entry rows.
 *         The initial result before segmentation smoothing is output to the .modelBegin.af.param file
 *         and the final result after segmentation smoothing is output to the .modelFinal.af.param file.
 *     </li>
 *     <li>
 *         Copy-ratio-model global-parameter files (.modelBegin.cr.param and .modelFinal.cr.param).
 *         These are tab-separated values (TSV) files with a SAM-style header containing a read group sample name,
 *         a row specifying the column headers contained in {@link ParameterDecileCollection.ParameterTableColumn},
 *         and the corresponding entry rows.
 *         The initial result before segmentation smoothing is output to the .modelBegin.cr.param file
 *         and the final result after segmentation smoothing is output to the .modelFinal.cr.param file.
 *     </li>
 *     <li>
 *         Copy-ratio segment file (.cr.seg).
 *         This is a tab-separated values (TSV) file with a SAM-style header containing a read group sample name, a sequence dictionary,
 *         a row specifying the column headers contained in {@link CopyRatioSegmentCollection.CopyRatioSegmentTableColumn},
 *         and the corresponding entry rows.
 *         It contains the segments from the .modelFinal.seg file converted to a format suitable for input to {@link CallCopyRatioSegments}.
 *     </li>
 *     <li>
 *         CBS-formatted .cr.igv.seg and .af.igv.seg files (compatible with IGV).
 *         These are tab-separated values (TSV) files with CBS-format column headers
 *         (see <a href="http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#CBS">
 *             http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#CBS</a>)
 *         and the corresponding entry rows that can be plotted using IGV (see
 *         <a href="https://software.broadinstitute.org/software/igv/SEG">
 *             https://software.broadinstitute.org/software/igv/SEG</a>).
 *         The posterior medians of the log2 copy ratio and minor-allele fraction are given in the SEGMENT_MEAN
 *         columns in the .cr.igv.seg and .af.igv.seg files, respectively.
 *     </li>
 *     <li>
 *         (Optional) Allelic-counts file containing the counts at sites genotyped as heterozygous in the case sample (.hets.tsv).
 *         This is a tab-separated values (TSV) file with a SAM-style header containing a read group sample name, a sequence dictionary,
 *         a row specifying the column headers contained in {@link AllelicCountCollection.AllelicCountTableColumn},
 *         and the corresponding entry rows.
 *         This is only output if allelic counts are provided as input.
 *     </li>
 *     <li>
 *         (Optional) Allelic-counts file containing the counts at sites genotyped as heterozygous in the matched-normal sample (.hets.normal.tsv).
 *         This is a tab-separated values (TSV) file with a SAM-style header containing a read group sample name, a sequence dictionary,
 *         a row specifying the column headers contained in {@link AllelicCountCollection.AllelicCountTableColumn},
 *         and the corresponding entry rows.
 *         This is only output if matched-normal allelic counts are provided as input.
 *     </li>
 * </ul>
 *
 * <h3>Usage examples</h3>
 *
 * <pre>
 *     gatk ModelSegments \
 *          --denoised-copy-ratios tumor.denoisedCR.tsv \
 *          --allelic-counts tumor.allelicCounts.tsv \
 *          --normal-allelic-counts normal.allelicCounts.tsv \
 *          --output-prefix tumor \
 *          -O output_dir
 * </pre>
 *
 * <pre>
 *     gatk ModelSegments \
 *          --denoised-copy-ratios normal.denoisedCR.tsv \
 *          --allelic-counts normal.allelicCounts.tsv \
 *          --output-prefix normal \
 *          -O output_dir
 * </pre>
 *
 * <pre>
 *     gatk ModelSegments \
 *          --allelic-counts tumor.allelicCounts.tsv \
 *          --normal-allelic-counts normal.allelicCounts.tsv \
 *          --output-prefix tumor \
 *          -O output_dir
 * </pre>
 *
 * <pre>
 *     gatk ModelSegments \
 *          --denoised-copy-ratios normal.denoisedCR.tsv \
 *          --output-prefix normal \
 *          -O output_dir
 * </pre>
 *
 * <pre>
 *     gatk ModelSegments \
 *          --allelic-counts tumor.allelicCounts.tsv \
 *          --output-prefix tumor \
 *          -O output_dir
 * </pre>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Models segmented copy ratios from denoised read counts and segmented minor-allele fractions from allelic counts",
        oneLineSummary = "Models segmented copy ratios from denoised read counts and segmented minor-allele fractions from allelic counts",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class SegmentJointSamples extends CommandLineProgram {
    @Argument(
            doc = "Input files containing denoised copy ratios (output of DenoiseReadCounts).  " +
                    "Sample order must match that of input allelic-counts files.",
            fullName = CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME,
            minElements = 1
    )
    private List<File> inputDenoisedCopyRatiosFiles = null;

    @Argument(
            doc = "Input files containing allelic counts (output of CollectAllelicCounts).  " +
                    "Sample order must match that of input denoised-copy-ratios files.",
            fullName = CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME,
            minElements = 1
    )
    private List<File> inputAllelicCountsFiles = null;

    @Argument(
            doc = "Input file containing allelic counts for a matched normal (output of CollectAllelicCounts).",
            fullName = CopyNumberStandardArgument.NORMAL_ALLELIC_COUNTS_FILE_LONG_NAME
    )
    private File inputNormalAllelicCountsFile = null;

    @Argument(
            doc = "Output file for multidimensional (i.e., copy-ratio and allele-fraction) segments.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputSegmentsFile;

    @ArgumentCollection
    private SomaticGenotypingArgumentCollection genotypingArguments = new SomaticGenotypingArgumentCollection();

    @ArgumentCollection
    private SomaticSegmentationArgumentCollection segmentationArguments = new SomaticSegmentationArgumentCollection();
    private final int maxNumSegmentsPerChromosome = segmentationArguments.maxNumSegmentsPerChromosome;
    private final double kernelVarianceCopyRatio = segmentationArguments.kernelVarianceCopyRatio;
    private final double kernelVarianceAlleleFraction = segmentationArguments.kernelVarianceAlleleFraction;
    private final double kernelScalingAlleleFraction = segmentationArguments.kernelScalingAlleleFraction;
    private final int kernelApproximationDimension = segmentationArguments.kernelApproximationDimension;
    private final List<Integer> windowSizes = segmentationArguments.windowSizes;
    private final double numChangepointsPenaltyFactor = segmentationArguments.numChangepointsPenaltyFactor;

    private void logHeapUsage(final String phase) {
        final int mb = 1024 * 1024;
        final Runtime runtime = Runtime.getRuntime();
        logger.info("Used memory after " + phase + ": " + (runtime.totalMemory() - runtime.freeMemory()) / mb);
    }

    @Override
    protected Object doWork() {
        validateArguments();

        //read input files (return null if not available) and validate metadata
        logHeapUsage("read input files");
        final List<CopyRatioCollection> denoisedCopyRatiosList = inputDenoisedCopyRatiosFiles.stream()
                .map(CopyRatioCollection::new)
                .collect(Collectors.toList());
        final List<AllelicCountCollection> allelicCountsList = inputAllelicCountsFiles.stream()
                .map(AllelicCountCollection::new)
                .collect(Collectors.toList());
        final AllelicCountCollection normalAllelicCounts = new AllelicCountCollection(inputNormalAllelicCountsFile);

        //validate metadata

        //genotype hets (return empty collection containing only metadata if no allelic counts available)
        final List<AllelicCountCollection> hetAllelicCountsList = IntStream.range(0, denoisedCopyRatiosList.size()).boxed()
                .map(i -> NaiveHeterozygousPileupGenotypingUtils.genotypeHets(
                        denoisedCopyRatiosList.get(i), allelicCountsList.get(i), normalAllelicCounts, genotypingArguments)
                        .getHetAllelicCounts())
                .collect(Collectors.toList());
        logHeapUsage("genotype hets");

        //if denoised copy ratios are still null at this point, we assign an empty collection containing only metadata

        //at this point, both denoisedCopyRatios and hetAllelicCounts are non-null, but may be empty;
        //perform one-dimensional or multidimensional segmentation as appropriate
//        final MultidimensionalSegmentCollection multidimensionalSegments;
//        if (!denoisedCopyRatios.getRecords().isEmpty() && hetAllelicCounts.getRecords().isEmpty()) {
//            final CopyRatioSegmentCollection copyRatioSegments = performCopyRatioSegmentation(denoisedCopyRatios);
//            multidimensionalSegments = new MultidimensionalSegmentCollection(
//                    copyRatioSegments.getMetadata(),
//                    copyRatioSegments.getRecords().stream()
//                            .map(s -> new MultidimensionalSegment(s.getInterval(), s.getNumPoints(), 0))
//                            .collect(Collectors.toList()));
//        } else if (denoisedCopyRatios.getRecords().isEmpty() && !hetAllelicCounts.getRecords().isEmpty()) {
//            final AlleleFractionSegmentCollection alleleFractionSegments = performAlleleFractionSegmentation(hetAllelicCounts);
//            multidimensionalSegments = new MultidimensionalSegmentCollection(
//                    alleleFractionSegments.getMetadata(),
//                    alleleFractionSegments.getRecords().stream()
//                            .map(s -> new MultidimensionalSegment(s.getInterval(), 0, s.getNumPoints()))
//                            .collect(Collectors.toList()));
//        } else {
//            multidimensionalSegments = new MultidimensionalKernelSegmenter(denoisedCopyRatios, hetAllelicCounts)
//                    .findSegmentation(maxNumSegmentsPerChromosome,
//                            kernelVarianceCopyRatio, kernelVarianceAlleleFraction, kernelScalingAlleleFraction, kernelApproximationDimension,
//                            ImmutableSet.copyOf(windowSizes).asList(),
//                            numChangepointsPenaltyFactor, numChangepointsPenaltyFactor);
//        }
        final SimpleIntervalCollection segments = new MultisampleMultidimensionalKernelSegmenter(denoisedCopyRatiosList, hetAllelicCountsList)
                .findSegmentation(maxNumSegmentsPerChromosome,
                        kernelVarianceCopyRatio, kernelVarianceAlleleFraction, kernelScalingAlleleFraction, kernelApproximationDimension,
                        ImmutableSet.copyOf(windowSizes).asList(),
                        numChangepointsPenaltyFactor, numChangepointsPenaltyFactor);
        logHeapUsage("multidimensional segmentation");

        //write segments to file
        segments.write(outputSegmentsFile);

        logger.info(String.format("%s complete.", getClass().getSimpleName()));

        return null;
    }

    private void validateArguments() {
        Utils.validateArg(inputDenoisedCopyRatiosFiles.size() == inputAllelicCountsFiles.size(),
                "Number of denoised-copy-ratios files and allelic-counts files must be equal.");

        final int maxNumInputFiles = 2 * inputDenoisedCopyRatiosFiles.size() + 1;
        final List<File> inputFiles = new ArrayList<>(maxNumInputFiles);
        inputFiles.addAll(inputDenoisedCopyRatiosFiles);
        inputFiles.addAll(inputAllelicCountsFiles);
        inputFiles.add(inputNormalAllelicCountsFile);
        CopyNumberArgumentValidationUtils.validateInputs(
                inputFiles.toArray(new File[maxNumInputFiles]));
        CopyNumberArgumentValidationUtils.validateOutputFiles(outputSegmentsFile);
    }
}
