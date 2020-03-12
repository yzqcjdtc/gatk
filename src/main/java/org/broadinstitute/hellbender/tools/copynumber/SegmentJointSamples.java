package org.broadinstitute.hellbender.tools.copynumber;

import com.google.common.collect.ImmutableSet;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.arguments.SomaticGenotypingArgumentCollection;
import org.broadinstitute.hellbender.tools.copynumber.arguments.SomaticSegmentationArgumentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.*;
import org.broadinstitute.hellbender.tools.copynumber.models.AlleleFractionModeller;
import org.broadinstitute.hellbender.tools.copynumber.models.CopyRatioModeller;
import org.broadinstitute.hellbender.tools.copynumber.segmentation.MultisampleMultidimensionalKernelSegmenter;
import org.broadinstitute.hellbender.tools.copynumber.utils.genotyping.NaiveHeterozygousPileupGenotypingUtils;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenter;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Finds common segments across multiple case samples using denoised copy ratios and allelic counts.
 * This segmentation can be used as input to subsequent, individual runs of {@link ModelSegments}
 * on each of the case samples.
 *
 * <p>
 *     Possible data inputs are: 1) denoised copy ratios for the case samples, 2) allelic counts for the case samples,
 *     and 3) allelic counts for a matched-normal sample.  All available inputs will be used to to perform
 *     segmentation.
 * </p>
 *
 * <p>
 *     The first step is to genotype heterozygous sites, as the allelic counts at these sites will subsequently
 *     be modeled to infer segmented minor-allele fraction by {@link ModelSegments}.
 *     We perform a relatively simple and naive genotyping based on the allele counts (i.e., pileups), which is
 *     controlled by a small number of parameters ({@code minimum-total-allele-count},
 *     {@code genotyping-homozygous-log-ratio-threshold}, and {@code genotyping-homozygous-log-ratio-threshold}).
 *     If the matched normal is available, its allelic counts will be used to genotype the sites, and
 *     we will simply assume these genotypes are the same in the case samples.  (This can be critical, for example,
 *     for determining sites with loss of heterozygosity in high purity case samples; such sites will be genotyped as
 *     homozygous if the matched-normal sample is not available.)  If no matched normal is available, we will use
 *     the intersection of all heterozygous sites found across all case samples.
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
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         List of denoised-copy-ratios files from {@link DenoiseReadCounts}.
 *     </li>
 *     <li>
 *         List of allelic-counts file from {@link CollectAllelicCounts}.
 *     </li>
 *     <li>
 *         Matched-normal allelic-counts file from {@link CollectAllelicCounts}.
 *     </li>
 * </ul>
 *
 * <h3>Outputs</h3>
 *
 * <ul>
 *     <li>
 *         Segments file.
 *         This is a tab-separated values (TSV) file with a SAM-style header containing a sequence dictionary,
 *         a row specifying the column headers contained in {@link SimpleIntervalCollection.SimpleIntervalTableColumn},
 *         and the corresponding entry rows.
 *         This segmentation can be used as input to subsequent, individual runs of {@link ModelSegments} on each of
 *         the case samples.
 *     </li>
 * </ul>
 *
 * <h3>Usage examples</h3>
 *
 * <pre>
 *     gatk SegmentJointSamples \
 *          --denoised-copy-ratios tumor-1.denoisedCR.tsv \
 *          ...
 *          --denoised-copy-ratios tumor-N.denoisedCR.tsv \
 *          --allelic-counts tumor-1.allelicCounts.tsv \
 *          ...
 *          --allelic-counts tumor-N.allelicCounts.tsv \
 *          --normal-allelic-counts normal.allelicCounts.tsv \
 *          -O tumor.joint.seg
 * </pre>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Finds common segments across multiple samples using denoised copy ratios and allelic counts",
        oneLineSummary = "Finds common segments across multiple samples using denoised copy ratios and allelic counts",
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
        logger.info("Used memory (MB) after " + phase + ": " + (runtime.totalMemory() - runtime.freeMemory()) / mb);
    }

    @Override
    protected Object doWork() {
        validateArguments();

        //read input files (return null if not available) and validate metadata
        logHeapUsage("reading input files");
        final List<CopyRatioCollection> denoisedCopyRatiosList = inputDenoisedCopyRatiosFiles.stream()
                .map(f -> readOptionalFileOrNull(f, CopyRatioCollection::new))
                .collect(Collectors.toList());
        final List<AllelicCountCollection> allelicCountsList = inputAllelicCountsFiles.stream()
                .map(f -> readOptionalFileOrNull(f, AllelicCountCollection::new))
                .collect(Collectors.toList());
        final AllelicCountCollection normalAllelicCounts = new AllelicCountCollection(inputNormalAllelicCountsFile);

        //validate metadata
        //TODO

        //genotype hets (return empty collection containing only metadata if no allelic counts available)
        //TODO could eliminate some redundant computation/logging
        final List<AllelicCountCollection> hetAllelicCountsList = IntStream.range(0, denoisedCopyRatiosList.size()).boxed()
                .map(i -> NaiveHeterozygousPileupGenotypingUtils.genotypeHets(
                        denoisedCopyRatiosList.get(i), allelicCountsList.get(i), normalAllelicCounts, genotypingArguments)
                        .getHetAllelicCounts())
                .collect(Collectors.toList());
        logHeapUsage("genotyping hets");

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
        logHeapUsage("performing multidimensional segmentation");

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

    private <T> T readOptionalFileOrNull(final File file,
                                         final Function<File, T> read) {
        if (file == null) {
            return null;
        }
        logger.info(String.format("Reading file (%s)...", file));
        return read.apply(file);
    }
}
