package org.broadinstitute.hellbender.tools.copynumber.utils.genotyping;

import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.arguments.ModelSegmentsGenotypingArgumentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AllelicCount;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.stream.Collectors;

/**
 * Naive methods for binomial genotyping of heterozygous sites from pileup allele counts.
 * Filters for total count and overlap with copy-ratio intervals are also implemented.
 */
public final class NaiveHeterozygousPileupGenotypingUtils {
    private static final Logger logger = LogManager.getLogger(NaiveHeterozygousPileupGenotypingUtils.class);

    public static class NaiveHeterozygousPileupGenotypingResult {
        private final AllelicCountCollection hetAllelicCounts;
        private final AllelicCountCollection hetNormalAllelicCounts;


        public NaiveHeterozygousPileupGenotypingResult(final AllelicCountCollection hetAllelicCounts,
                                                       final AllelicCountCollection hetNormalAllelicCounts) {
            this.hetAllelicCounts = Utils.nonNull(hetAllelicCounts);
            this.hetNormalAllelicCounts = hetNormalAllelicCounts;
        }

        public AllelicCountCollection getHetAllelicCounts() {
            return hetAllelicCounts;
        }

        public AllelicCountCollection getHetNormalAllelicCounts() {
            return hetNormalAllelicCounts;
        }
    }

    /**
     * This method can be called either in matched-normal or case-only mode.
     * Optional inputs should be set to {@code null} as appropriate.
     * Validation of sequence dictionaries and sample names will be performed.
     * @param denoisedCopyRatios    if not {@code null}, sites not overlapping with copy-ratio intervals will be filtered out
     * @param allelicCounts         never {@code null}
     * @param normalAllelicCounts   if not {@code null}, sites that are homozygous in the normal will be filtered out
     * @return {@link NaiveHeterozygousPileupGenotypingResult}, with hetNormalAllelicCounts set to {@code null} for case-only mode
     */
    public static NaiveHeterozygousPileupGenotypingResult genotypeHets(final CopyRatioCollection denoisedCopyRatios,
                                                                       final AllelicCountCollection allelicCounts,
                                                                       final AllelicCountCollection normalAllelicCounts,
                                                                       final ModelSegmentsGenotypingArgumentCollection genotypingArguments) {
        Utils.nonNull(allelicCounts);
        Utils.nonNull(genotypingArguments);
        final SampleLocatableMetadata metadata = CopyNumberArgumentValidationUtils.getValidatedMetadata(
                denoisedCopyRatios, allelicCounts);
        final int minTotalAlleleCountCase = genotypingArguments.minTotalAlleleCountCase;
        final int minTotalAlleleCountNormal = genotypingArguments.minTotalAlleleCountNormal;
        final double genotypingHomozygousLogRatioThreshold = genotypingArguments.genotypingHomozygousLogRatioThreshold;
        final double genotypingBaseErrorRate = genotypingArguments.genotypingBaseErrorRate;

        logger.info("Genotyping heterozygous sites from available allelic counts...");

        AllelicCountCollection filteredAllelicCounts = allelicCounts;

        //filter on total count in case sample
        logger.info(String.format("Filtering allelic counts with total count less than %d...", minTotalAlleleCountCase));
        filteredAllelicCounts = new AllelicCountCollection(
                metadata,
                filteredAllelicCounts.getRecords().stream()
                        .filter(ac -> ac.getTotalReadCount() >= minTotalAlleleCountCase)
                        .collect(Collectors.toList()));
        logger.info(String.format("Retained %d / %d sites after filtering on total count...",
                filteredAllelicCounts.size(), allelicCounts.size()));

        //filter on overlap with copy-ratio intervals, if available
        if (denoisedCopyRatios != null) {
            logger.info("Filtering allelic-count sites not overlapping with copy-ratio intervals...");
            filteredAllelicCounts = new AllelicCountCollection(
                    metadata,
                    filteredAllelicCounts.getRecords().stream()
                            .filter(ac -> denoisedCopyRatios.getOverlapDetector().overlapsAny(ac))
                            .collect(Collectors.toList()));
            logger.info(String.format("Retained %d / %d sites after filtering on overlap with copy-ratio intervals...",
                    filteredAllelicCounts.size(), allelicCounts.size()));
        }

        final AllelicCountCollection hetAllelicCounts;
        final AllelicCountCollection hetNormalAllelicCounts;
        if (normalAllelicCounts == null) {
            //filter on homozygosity in case sample
            logger.info("No matched normal was provided, not running in matched-normal mode...");

            logger.info("Performing binomial testing and filtering homozygous allelic counts...");
            hetAllelicCounts = new AllelicCountCollection(
                    metadata,
                    filteredAllelicCounts.getRecords().stream()
                            .filter(ac -> calculateHomozygousLogRatio(ac, genotypingBaseErrorRate) < genotypingHomozygousLogRatioThreshold)
                            .collect(Collectors.toList()));
            logger.info(String.format("Retained %d / %d sites after testing for heterozygosity...",
                    hetAllelicCounts.size(), allelicCounts.size()));
            hetNormalAllelicCounts = null;
        } else {
            //use matched normal
            logger.info("Matched normal was provided, running in matched-normal mode...");
            logger.info("Performing binomial testing and filtering homozygous allelic counts in matched normal...");
            if (!normalAllelicCounts.getIntervals().equals(allelicCounts.getIntervals())) {
                throw new UserException.BadInput("Allelic-count sites in case sample and matched normal do not match. " +
                        "Run CollectAllelicCounts using the same interval list of sites for both samples.");
            }
            final SampleLocatableMetadata normalMetadata = normalAllelicCounts.getMetadata();
            if (!CopyNumberArgumentValidationUtils.isSameDictionary(
                    normalMetadata.getSequenceDictionary(),
                    metadata.getSequenceDictionary())) {
                logger.warn("Sequence dictionaries in allelic-count files do not match.");
            }

            //filter on total count in matched normal
            logger.info(String.format("Filtering allelic counts in matched normal with total count less than %d...", minTotalAlleleCountNormal));
            AllelicCountCollection filteredNormalAllelicCounts = new AllelicCountCollection(
                    normalMetadata,
                    normalAllelicCounts.getRecords().stream()
                            .filter(ac -> ac.getTotalReadCount() >= minTotalAlleleCountNormal)
                            .collect(Collectors.toList()));
            logger.info(String.format("Retained %d / %d sites in matched normal after filtering on total count...",
                    filteredNormalAllelicCounts.size(), normalAllelicCounts.size()));

            //filter matched normal on overlap with copy-ratio intervals, if available
            if (denoisedCopyRatios != null) {
                logger.info("Filtering allelic-count sites in matched normal not overlapping with copy-ratio intervals...");
                filteredNormalAllelicCounts = new AllelicCountCollection(
                        normalMetadata,
                        filteredNormalAllelicCounts.getRecords().stream()
                                .filter(ac -> denoisedCopyRatios.getOverlapDetector().overlapsAny(ac))
                                .collect(Collectors.toList()));
                logger.info(String.format("Retained %d / %d sites in matched normal after filtering on overlap with copy-ratio intervals...",
                        filteredNormalAllelicCounts.size(), normalAllelicCounts.size()));
            }

            //filter on homozygosity in matched normal
            hetNormalAllelicCounts = new AllelicCountCollection(
                    normalMetadata,
                    filteredNormalAllelicCounts.getRecords().stream()
                            .filter(ac -> calculateHomozygousLogRatio(ac, genotypingBaseErrorRate) < genotypingHomozygousLogRatioThreshold)
                            .collect(Collectors.toList()));
            logger.info(String.format("Retained %d / %d sites in matched normal after testing for heterozygosity...",
                    hetNormalAllelicCounts.size(), normalAllelicCounts.size()));

            //retrieve sites in case sample
            logger.info("Retrieving allelic counts at these sites in case sample...");
            hetAllelicCounts = new AllelicCountCollection(
                    metadata,
                    filteredAllelicCounts.getRecords().stream()
                            .filter(ac -> hetNormalAllelicCounts.getOverlapDetector().overlapsAny(ac))
                            .collect(Collectors.toList()));
        }
        return new NaiveHeterozygousPileupGenotypingResult(hetAllelicCounts, hetNormalAllelicCounts);
    }

    private static double calculateHomozygousLogRatio(final AllelicCount allelicCount,
                                                      final double genotypingBaseErrorRate) {
        final int r = allelicCount.getRefReadCount();
        final int n = allelicCount.getTotalReadCount();
        final double betaAll = Beta.regularizedBeta(1, r + 1, n - r + 1);
        final double betaError = Beta.regularizedBeta(genotypingBaseErrorRate, r + 1, n - r + 1);
        final double betaOneMinusError = Beta.regularizedBeta(1 - genotypingBaseErrorRate, r + 1, n - r + 1);
        final double betaHom = betaError + betaAll - betaOneMinusError;
        final double betaHet = betaOneMinusError - betaError;
        return FastMath.log(betaHom) - FastMath.log(betaHet);
    }
}
