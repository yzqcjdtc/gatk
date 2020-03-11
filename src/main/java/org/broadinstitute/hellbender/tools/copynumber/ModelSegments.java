package org.broadinstitute.hellbender.tools.copynumber;

import com.google.common.collect.ImmutableSet;
import htsjdk.samtools.util.OverlapDetector;
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
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.LegacySegment;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.MultidimensionalSegment;
import org.broadinstitute.hellbender.tools.copynumber.models.AlleleFractionModeller;
import org.broadinstitute.hellbender.tools.copynumber.models.AlleleFractionPrior;
import org.broadinstitute.hellbender.tools.copynumber.models.CopyRatioModeller;
import org.broadinstitute.hellbender.tools.copynumber.models.MultidimensionalModeller;
import org.broadinstitute.hellbender.tools.copynumber.segmentation.AlleleFractionKernelSegmenter;
import org.broadinstitute.hellbender.tools.copynumber.segmentation.CopyRatioKernelSegmenter;
import org.broadinstitute.hellbender.tools.copynumber.segmentation.MultidimensionalKernelSegmenter;
import org.broadinstitute.hellbender.tools.copynumber.utils.genotyping.NaiveHeterozygousPileupGenotypingUtils;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenter;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

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
public final class ModelSegments extends CommandLineProgram {
    //filename tags for output
    public static final String HET_ALLELIC_COUNTS_FILE_SUFFIX = ".hets.tsv";
    public static final String NORMAL_HET_ALLELIC_COUNTS_FILE_SUFFIX = ".hets.normal.tsv";
    public static final String SEGMENTS_FILE_SUFFIX = ".seg";
    public static final String BEGIN_FIT_FILE_TAG = ".modelBegin";
    public static final String FINAL_FIT_FILE_TAG = ".modelFinal";
    public static final String COPY_RATIO_MODEL_PARAMETER_FILE_SUFFIX = ".cr.param";
    public static final String ALLELE_FRACTION_MODEL_PARAMETER_FILE_SUFFIX = ".af.param";
    public static final String COPY_RATIO_SEGMENTS_FOR_CALLER_FILE_SUFFIX = ".cr" + SEGMENTS_FILE_SUFFIX;
    public static final String COPY_RATIO_LEGACY_SEGMENTS_FILE_SUFFIX = ".cr.igv" + SEGMENTS_FILE_SUFFIX;
    public static final String ALLELE_FRACTION_LEGACY_SEGMENTS_FILE_SUFFIX = ".af.igv" + SEGMENTS_FILE_SUFFIX;

    @Argument(
            doc = "Input file containing denoised copy ratios (output of DenoiseReadCounts).",
            fullName = CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME,
            optional = true
    )
    private File inputDenoisedCopyRatiosFile = null;

    @Argument(
            doc = "Input file containing allelic counts (output of CollectAllelicCounts).",
            fullName = CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME,
            optional = true
    )
    private File inputAllelicCountsFile = null;

    @Argument(
            doc = "Input file containing allelic counts for a matched normal (output of CollectAllelicCounts).",
            fullName = CopyNumberStandardArgument.NORMAL_ALLELIC_COUNTS_FILE_LONG_NAME,
            optional = true
    )
    private File inputNormalAllelicCountsFile = null;

    @Argument(
            doc = "Prefix for output filenames.",
            fullName =  CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME
    )
    private String outputPrefix;

    @Argument(
            doc = "Output directory.  This will be created if it does not exist.",
            fullName =  StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputDir;

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

    @ArgumentCollection
    private SomaticModelingArgumentCollection modelingArguments = new SomaticModelingArgumentCollection();
    private final double minorAlleleFractionPriorAlpha = modelingArguments.minorAlleleFractionPriorAlpha;
    private final int numSamplesCopyRatio = modelingArguments.numSamplesCopyRatio;
    private final int numBurnInCopyRatio = modelingArguments.numBurnInCopyRatio;
    private final int numSamplesAlleleFraction = modelingArguments.numSamplesAlleleFraction;
    private final int numBurnInAlleleFraction = modelingArguments.numBurnInAlleleFraction;
    private final double smoothingCredibleIntervalThresholdCopyRatio = modelingArguments.smoothingCredibleIntervalThresholdCopyRatio;
    private final double smoothingCredibleIntervalThresholdAlleleFraction = modelingArguments.smoothingCredibleIntervalThresholdAlleleFraction;
    private final int maxNumSmoothingIterations = modelingArguments.maxNumSmoothingIterations;
    private final int numSmoothingIterationsPerFit = modelingArguments.numSmoothingIterationsPerFit;

    @Override
    protected Object doWork() {
        validateArguments();

        //read input files (return null if not available) and validate metadata
        CopyRatioCollection denoisedCopyRatios = readOptionalFileOrNull(inputDenoisedCopyRatiosFile, CopyRatioCollection::new);
        final AllelicCountCollection allelicCounts = readOptionalFileOrNull(inputAllelicCountsFile, AllelicCountCollection::new);
        final AllelicCountCollection normalAllelicCounts = readOptionalFileOrNull(inputNormalAllelicCountsFile, AllelicCountCollection::new);
        final SampleLocatableMetadata metadata = CopyNumberArgumentValidationUtils.getValidatedMetadata(denoisedCopyRatios, allelicCounts);

        //genotype hets
        //hetAllelicCounts is set to an empty collection containing only metadata if no allelic counts are available;
        //output allelic-count files containing hets for the case and the matched-normal are only written when available
        final AllelicCountCollection hetAllelicCounts;
        if (allelicCounts == null) {
            hetAllelicCounts = new AllelicCountCollection(metadata, Collections.emptyList());
        } else {
            final NaiveHeterozygousPileupGenotypingUtils.NaiveHeterozygousPileupGenotypingResult genotypingResult =
                    NaiveHeterozygousPileupGenotypingUtils.genotypeHets(
                            denoisedCopyRatios, allelicCounts, normalAllelicCounts, genotypingArguments);
            hetAllelicCounts = genotypingResult.getHetAllelicCounts();
            if (normalAllelicCounts == null) {
                //case-only mode
                final File hetAllelicCountsFile = new File(outputDir, outputPrefix + HET_ALLELIC_COUNTS_FILE_SUFFIX);
                logger.info(String.format("Writing heterozygous allelic counts to %s...", hetAllelicCountsFile.getAbsolutePath()));
                hetAllelicCounts.write(hetAllelicCountsFile);
            } else {
                //matched-normal mode
                final File hetNormalAllelicCountsFile = new File(outputDir, outputPrefix + NORMAL_HET_ALLELIC_COUNTS_FILE_SUFFIX);
                logger.info(String.format("Writing heterozygous allelic counts for matched normal to %s...", hetNormalAllelicCountsFile.getAbsolutePath()));
                genotypingResult.getHetNormalAllelicCounts().write(hetNormalAllelicCountsFile);

                final File hetAllelicCountsFile = new File(outputDir, outputPrefix + HET_ALLELIC_COUNTS_FILE_SUFFIX);
                logger.info(String.format("Writing allelic counts for case sample at heterozygous sites in matched normal to %s...", hetAllelicCountsFile.getAbsolutePath()));
                hetAllelicCounts.write(hetAllelicCountsFile);
            }
        }

        //if denoised copy ratios are still null at this point, we assign an empty collection containing only metadata
        if (denoisedCopyRatios == null) {
            denoisedCopyRatios = new CopyRatioCollection(metadata, Collections.emptyList());
        }

        //at this point, both denoisedCopyRatios and hetAllelicCounts are non-null, but may be empty;
        //perform one-dimensional or multidimensional segmentation as appropriate
        final MultidimensionalSegmentCollection multidimensionalSegments;
        if (!denoisedCopyRatios.getRecords().isEmpty() && hetAllelicCounts.getRecords().isEmpty()) {
            final CopyRatioSegmentCollection copyRatioSegments = performCopyRatioSegmentation(denoisedCopyRatios);
            multidimensionalSegments = new MultidimensionalSegmentCollection(
                    copyRatioSegments.getMetadata(),
                    copyRatioSegments.getRecords().stream()
                            .map(s -> new MultidimensionalSegment(s.getInterval(), s.getNumPoints(), 0))
                            .collect(Collectors.toList()));
        } else if (denoisedCopyRatios.getRecords().isEmpty() && !hetAllelicCounts.getRecords().isEmpty()) {
            final AlleleFractionSegmentCollection alleleFractionSegments = performAlleleFractionSegmentation(hetAllelicCounts);
            multidimensionalSegments = new MultidimensionalSegmentCollection(
                    alleleFractionSegments.getMetadata(),
                    alleleFractionSegments.getRecords().stream()
                            .map(s -> new MultidimensionalSegment(s.getInterval(), 0, s.getNumPoints()))
                            .collect(Collectors.toList()));
        } else {
            multidimensionalSegments = new MultidimensionalKernelSegmenter(denoisedCopyRatios, hetAllelicCounts)
                    .findSegmentation(maxNumSegmentsPerChromosome,
                            kernelVarianceCopyRatio, kernelVarianceAlleleFraction, kernelScalingAlleleFraction, kernelApproximationDimension,
                            ImmutableSet.copyOf(windowSizes).asList(),
                            numChangepointsPenaltyFactor, numChangepointsPenaltyFactor);
        }

        logger.info("Modeling available denoised copy ratios and heterozygous allelic counts...");
        //initial MCMC model fitting performed by MultidimensionalModeller constructor
        final AlleleFractionPrior alleleFractionPrior = new AlleleFractionPrior(minorAlleleFractionPriorAlpha);
        final MultidimensionalModeller modeller = new MultidimensionalModeller(
                multidimensionalSegments, denoisedCopyRatios, hetAllelicCounts, alleleFractionPrior,
                numSamplesCopyRatio, numBurnInCopyRatio,
                numSamplesAlleleFraction, numBurnInAlleleFraction);

        //write initial segments and parameters to file
        writeModeledSegmentsAndParameterFiles(modeller, BEGIN_FIT_FILE_TAG);

        //segmentation smoothing
        modeller.smoothSegments(
                maxNumSmoothingIterations, numSmoothingIterationsPerFit,
                smoothingCredibleIntervalThresholdCopyRatio, smoothingCredibleIntervalThresholdAlleleFraction);

        //write final segments and parameters to file
        writeModeledSegmentsAndParameterFiles(modeller, FINAL_FIT_FILE_TAG);

        //write final segments for copy-ratio caller (TODO remove this and MEAN_LOG2_COPY_RATIO column when new caller is available)
        final OverlapDetector<CopyRatio> copyRatioMidpointOverlapDetector = denoisedCopyRatios.getMidpointOverlapDetector();
        final CopyRatioSegmentCollection copyRatioSegmentsFinal = new CopyRatioSegmentCollection(
                modeller.getModeledSegments().getMetadata(),
                modeller.getModeledSegments().getIntervals().stream()
                        .map(s -> new CopyRatioSegment(s, new ArrayList<>(copyRatioMidpointOverlapDetector.getOverlaps(s))))
                        .collect(Collectors.toList()));
        writeSegments(copyRatioSegmentsFinal, COPY_RATIO_SEGMENTS_FOR_CALLER_FILE_SUFFIX);

        //write IGV-compatible files
        final LegacySegmentCollection copyRatioLegacySegments = new LegacySegmentCollection(
                metadata,
                modeller.getModeledSegments().getRecords().stream()
                        .map(s -> new LegacySegment(
                                metadata.getSampleName(),
                                s.getInterval(),
                                s.getNumPointsCopyRatio(),
                                s.getLog2CopyRatioSimplePosteriorSummary().getDecile50()))
                        .collect(Collectors.toList()));
        final LegacySegmentCollection alleleFractionLegacySegments = new LegacySegmentCollection(
                metadata,
                modeller.getModeledSegments().getRecords().stream()
                        .map(s -> new LegacySegment(
                                metadata.getSampleName(),
                                s.getInterval(),
                                s.getNumPointsAlleleFraction(),
                                s.getMinorAlleleFractionSimplePosteriorSummary().getDecile50()))
                        .collect(Collectors.toList()));
        writeSegments(copyRatioLegacySegments, COPY_RATIO_LEGACY_SEGMENTS_FILE_SUFFIX);
        writeSegments(alleleFractionLegacySegments, ALLELE_FRACTION_LEGACY_SEGMENTS_FILE_SUFFIX);

        logger.info(String.format("%s complete.", getClass().getSimpleName()));

        return null;
    }

    private void validateArguments() {
        Utils.validateArg(!(inputDenoisedCopyRatiosFile == null && inputAllelicCountsFile == null),
                "Must provide at least a denoised-copy-ratios file or an allelic-counts file.");
        Utils.validateArg(!(inputAllelicCountsFile == null && inputNormalAllelicCountsFile != null),
                "Must provide an allelic-counts file for the case sample to run in matched-normal mode.");

        CopyNumberArgumentValidationUtils.validateInputs(
                inputDenoisedCopyRatiosFile,
                inputAllelicCountsFile,
                inputNormalAllelicCountsFile);
        Utils.nonEmpty(outputPrefix);
        CopyNumberArgumentValidationUtils.validateAndPrepareOutputDirectories(outputDir);

        modelingArguments.validateArguments();
    }

    private <T> T readOptionalFileOrNull(final File file,
                                         final Function<File, T> read) {
        if (file == null) {
            return null;
        }
        logger.info(String.format("Reading file (%s)...", file));
        return read.apply(file);
    }

    private CopyRatioSegmentCollection performCopyRatioSegmentation(final CopyRatioCollection denoisedCopyRatios) {
        logger.info("Starting segmentation of denoised copy ratios...");
        return new CopyRatioKernelSegmenter(denoisedCopyRatios)
                .findSegmentation(maxNumSegmentsPerChromosome, kernelVarianceCopyRatio, kernelApproximationDimension,
                        ImmutableSet.copyOf(windowSizes).asList(),
                        numChangepointsPenaltyFactor, numChangepointsPenaltyFactor);
    }

    private AlleleFractionSegmentCollection performAlleleFractionSegmentation(final AllelicCountCollection hetAllelicCounts) {
        logger.info("Starting segmentation of heterozygous allelic counts...");
        return new AlleleFractionKernelSegmenter(hetAllelicCounts)
                .findSegmentation(maxNumSegmentsPerChromosome, kernelVarianceAlleleFraction, kernelApproximationDimension,
                        ImmutableSet.copyOf(windowSizes).asList(),
                        numChangepointsPenaltyFactor, numChangepointsPenaltyFactor);
    }

    private void writeModeledSegmentsAndParameterFiles(final MultidimensionalModeller modeller,
                                                       final String fileTag) {
        final ModeledSegmentCollection modeledSegments = modeller.getModeledSegments();
        writeSegments(modeledSegments, fileTag + SEGMENTS_FILE_SUFFIX);
        final File copyRatioParameterFile = new File(outputDir, outputPrefix + fileTag + COPY_RATIO_MODEL_PARAMETER_FILE_SUFFIX);
        final File alleleFractionParameterFile = new File(outputDir, outputPrefix + fileTag + ALLELE_FRACTION_MODEL_PARAMETER_FILE_SUFFIX);
        modeller.writeModelParameterFiles(copyRatioParameterFile, alleleFractionParameterFile);
    }

    private void writeSegments(final AbstractRecordCollection<?, ?> segments,
                               final String fileSuffix) {
        final File segmentsFile = new File(outputDir, outputPrefix + fileSuffix);
        logger.info(String.format("Writing segments to %s...", segmentsFile.getAbsolutePath()));
        segments.write(segmentsFile);
    }
}
