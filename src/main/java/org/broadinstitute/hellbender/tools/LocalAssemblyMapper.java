package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SvCigarUtils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.minimap2.MiniMap2Aligner;
import org.broadinstitute.hellbender.utils.minimap2.MiniMap2Alignment;
import org.broadinstitute.hellbender.utils.minimap2.MiniMap2Index;
import org.broadinstitute.hellbender.utils.read.CigarUtils;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

@DocumentedFeature
@CommandLineProgramProperties(
        summary = "experiment",
        oneLineSummary = "experiment",
        usageExample = "gatk LocalAssemblyMapper",
        programGroup = CoverageAnalysisProgramGroup.class
)
@BetaFeature
public class LocalAssemblyMapper extends CommandLineProgram {

    @Argument(fullName = "fasta-list", doc = "File containing list of fasta-files to align")
    private String inputList;

    @Argument(fullName = "ref-index", doc = "The MiniMap2 index for the reference")
    private String refIndex;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Write output to this file")
    public String outputName;

    @Override
    protected Object doWork() {
        final List<String> refNames;
        try ( final MiniMap2Index mm2Index = new MiniMap2Index(refIndex) ) {
            logger.info(refIndex + " loaded");
            refNames = mm2Index.getRefNames();
            try ( final BufferedReader listRdr = new BufferedReader(new FileReader(inputList)) ) {
                String fileName;
                long tigId = 0;
                try ( final BufferedWriter writer = new BufferedWriter(new FileWriter(outputName)) ) {
                    while ( (fileName = listRdr.readLine()) != null ) {
                        logger.info(fileName + " being processed");
                        final String eventName = fileName.replaceFirst(".fa$", ".");
                        try ( final BufferedReader fastaRdr =
                                      new BufferedReader(new InputStreamReader(BucketUtils.openFile(fileName))) ) {
                            final List<String> contigCalls = new ArrayList<>();
                            String fastaLine = fastaRdr.readLine();
                            if ( fastaLine == null ) {
                                logger.warn(fileName + " was empty.");
                                continue;
                            }
                            int line = 1;
                            while ( fastaLine != null ) {
                                if ( fastaLine.charAt(0) != '>' ) {
                                    throw new UserException("fasta file header does not start with '>' in line " +
                                            line + " of " + fileName);
                                }
                                final StringBuilder sb = new StringBuilder();
                                while ( (fastaLine = fastaRdr.readLine()) != null ) {
                                    line += 1;
                                    if ( fastaLine.charAt(0) == '>' ) break;
                                    sb.append(fastaLine);
                                }
                                contigCalls.add(sb.toString());
                            }
                            final MiniMap2Aligner aligner = new MiniMap2Aligner(mm2Index, MiniMap2Aligner.Preset.ASM20);
                            final List<List<MiniMap2Alignment>> allAlignments =
                                    aligner.alignSeqs(contigCalls, String::getBytes);
                            final int nContigs = contigCalls.size();
                            for ( int idx = 0; idx != nContigs; ++idx ) {
                                final List<MiniMap2Alignment> alignments = allAlignments.get(idx);
                                final int nAlignments = alignments.size();
                                for ( int alignmentIdx = 0; alignmentIdx != nAlignments; ++alignmentIdx ) {
                                    final MiniMap2Alignment alignment = alignments.get(alignmentIdx);
                                    writer.write(eventName + idx);
                                    writer.write('\t');
                                    final int refId = alignment.getRefId();

                                    final int samFlag =
                                            (alignment.isRevStrand() ? SAMFlag.READ_REVERSE_STRAND.intValue() : 0) +
                                            (alignmentIdx == 0 ? 0 : SAMFlag.SUPPLEMENTARY_ALIGNMENT.intValue()) +
                                            (refId >= 0 ? 0 : SAMFlag.READ_UNMAPPED.intValue());
                                    writer.write(Integer.toString(samFlag));
                                    writer.write('\t');
                                    writer.write(refId >= 0 ? refNames.get(refId) : "*");
                                    writer.write('\t');
                                    writer.write(Integer.toString(alignment.getRefStart() + 1));
                                    writer.write('\t');
                                    writer.write(Integer.toString(alignment.getMapQ()));
                                    writer.write('\t');
                                    writer.write(alignment.getCigar());
                                    writer.write("\t*\t0\t0\t");
                                    if ( alignmentIdx != 0 ) {
                                        writer.write('*');
                                    } else if ( alignment.isRevStrand() ) {
                                        writer.write(SequenceUtil.reverseComplement(contigCalls.get(idx)));
                                    } else {
                                        writer.write(contigCalls.get(idx));
                                    }
                                    writer.write("\t*");
                                    writer.newLine();
                                }
                            }
                        } catch ( final IOException ioe ) {
                            throw new UserException("Can't find fasta input file " + inputList, ioe);
                        }
                    }
                } catch ( final IOException ioe ) {
                    throw new UserException("Can't write output file " + outputName);
                }
            } catch ( final IOException ioe ) {
                throw new UserException("Can't find fasta-list input file " + inputList, ioe);
            }
        }

        return null;
    }
}
