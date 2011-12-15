package tdg.results;

import com.beust.jcommander.JCommander;
import com.google.common.base.Charsets;
import com.google.common.io.Files;
import tdg.Constants;
import tdg.utils.GeneticCode;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;

/**
 * Given the model files produced by ModelWriter (S.txt etc.), this parses those files and produced the distribution
 * of selection coefficients for mutations and substitutions (all & non-synonymous only for both).
 *
 * TODO: PB2 plots significant vs. other sites (i.e. the black + red histogram in paper)
 *
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 * @see ModelWriter
 * @see FitnessExtractor
 */
public class DistributionWriter {
    Options options;
    private static final String ROOT_DIR = "./";

    private static final double HI = 21D;
    private static final double LOW = -21D;
    private static final int MUTS_BINS = 168; // Means each bin is (21 - -21) / 168 = 0.25
    private static final int SUBS_BINS = 167; // So that the bins are symmetric around zero for substitutions (reversible)
    private static final int S_LIMIT = 10;

    public static void main(String[] args) throws Exception {
        Options o = new Options();
        new JCommander(o, args);
        
        DistributionWriter dw = new DistributionWriter();
        dw.options = o;
        dw.run();
    }
    
    private void run() throws Exception {
        // We need to know the genetic code to determine non-synonymous changes
        char[] aaCode = GeneticCode.CURRENT_CODE.toCharArray();

        double[] mutsAll, mutsNonSyn, subsAll, subsNonSyn;
        mutsAll = mutsNonSyn = new double[MUTS_BINS];
        subsAll = subsNonSyn = new double[SUBS_BINS];

        // Get the neutral mutation matrix, Q0
        BufferedReader Q0Reader = new BufferedReader(new FileReader(ROOT_DIR + Constants.Q0_FILENAME));
        String[] Q0Parts = Q0Reader.readLine().split("\\s");
        double[] Q0 = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];
        for (int i = 0; i < (GeneticCode.CODON_STATES * GeneticCode.CODON_STATES); i++) {
            Q0[i] = Double.parseDouble(Q0Parts[i]);
        }
        Q0Reader.close();

        BufferedReader SReader = new BufferedReader(new FileReader(ROOT_DIR + Constants.S_FILENAME));
        BufferedReader QSReader = new BufferedReader(new FileReader(ROOT_DIR + Constants.QS_FILENAME));
        BufferedReader PiReader = new BufferedReader(new FileReader(ROOT_DIR + Constants.PI_FILENAME));

        double mutsAllDenom = 0, mutsNonSynDenom = 0, subsAllDenom = 0, subsNonSynDenom = 0;
        
        String SLine;
        while ((SLine = SReader.readLine()) != null) {
            String QSLine = QSReader.readLine();
            String PiLine = PiReader.readLine();

            String[] SParts = SLine.split("\\s"); // 64 * 64 = 4096 entries
            String[] QSParts = QSLine.split("\\s"); // 64 * 64 = 4096 entries
            String[] PiParts = PiLine.split("\\s"); // 64 entries

            double mutsAllSiteSum = 0, mutsNonSynSiteSum = 0, subsAllSiteSum = 0, subsNonSynSiteSum = 0;
            
            for (int i = 0; i < (GeneticCode.CODON_STATES * GeneticCode.CODON_STATES); i++) {
                
                int codon = (i - i % GeneticCode.CODON_STATES) / GeneticCode.CODON_STATES;
                char fromResidue = aaCode[codon];
                char toResidue = aaCode[i % GeneticCode.CODON_STATES];
                boolean isCodonChange = codon != (i % GeneticCode.CODON_STATES);
                boolean isNonSynChange = fromResidue != toResidue;

                double piValue = Double.parseDouble(PiParts[codon]);
                double qValue = Double.parseDouble(QSParts[i]);
                
                double deltaS = Double.parseDouble(SParts[i]);
                deltaS = Math.max(deltaS, -S_LIMIT);
                deltaS = Math.min(deltaS, S_LIMIT);

                double val = deltaS - LOW;
                int mutsBin = (int) (MUTS_BINS * (val / (HI - LOW)));
                int subsBin = (int) (SUBS_BINS * (val / (HI - LOW)));

                if (isNonSynChange) {
                    mutsNonSynSiteSum += piValue * Q0[i];
                    subsNonSynSiteSum += piValue * qValue;

                    mutsNonSyn[mutsBin] += piValue * Q0[i];
                    subsNonSyn[subsBin] += piValue * qValue;
                } 
                
                if (isCodonChange) {
                    mutsAllSiteSum += piValue * Q0[i];
                    subsAllSiteSum += piValue * qValue;
                    
                    mutsAll[mutsBin] += piValue * Q0[i];
                    subsAll[subsBin] += piValue * qValue;
                }
            }

            mutsAllDenom += mutsAllSiteSum;
            mutsNonSynDenom += mutsNonSynSiteSum;
            subsAllDenom += subsAllSiteSum;
            subsNonSynDenom += subsNonSynSiteSum;
        }

        SReader.close();
        QSReader.close();
        PiReader.close();

        BufferedWriter mutsWriter = Files.newWriter(new File("distribution.mutations.txt"), Charsets.US_ASCII);
        for (int i = 0 ; i < MUTS_BINS; i++) {
            mutsAll[i] /= mutsAllDenom;
            mutsNonSyn[i] /= mutsNonSynDenom;
            mutsWriter.write(String.format("%s\t%s\n", mutsAll[i], mutsNonSyn[i]));
        }
        mutsWriter.close();

        BufferedWriter subsWriter = Files.newWriter(new File("distribution.substitutions.txt"), Charsets.US_ASCII);
        for (int i = 0 ; i < SUBS_BINS; i++) {
            subsAll[i] /= subsAllDenom;
            subsNonSyn[i] /= subsNonSynDenom;
            subsWriter.write(String.format("%s\t%s\n", subsAll[i], subsNonSyn[i]));
        }
        subsWriter.close();
    }

}
