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

public class SDistributionWriter {
    Options options;
    private static final String ROOT_DIR = "./";
    
    private static final int MUTS_BINS = 168;
    private static final int SUBS_BINS = 167;
    private static final int S_LIMIT = 10;

    public static void main(String[] args) throws Exception {
        Options o = new Options();
        new JCommander(o, args);
        
        SDistributionWriter dw = new SDistributionWriter();
        dw.options = o;
        dw.run();
    }
    
    private void run() throws Exception {
        char[] aaCode = GeneticCode.CURRENT_CODE.toCharArray();

        double hi = 21D;
        double low = -21D;

        double[] mutsAll = new double[MUTS_BINS];
        double[] mutsNonSyn = new double[MUTS_BINS];
        double[] subsAll = new double[SUBS_BINS];
        double[] subsNonSyn = new double[SUBS_BINS];
        
        BufferedReader buffS = new BufferedReader(new FileReader(ROOT_DIR + Constants.S_FILENAME));
		BufferedReader buffQ = new BufferedReader(new FileReader(ROOT_DIR + Constants.QS_FILENAME));
		BufferedReader buffPi = new BufferedReader(new FileReader(ROOT_DIR + Constants.PI_FILENAME));
        BufferedReader buffQ0 = new BufferedReader(new FileReader(ROOT_DIR + Constants.Q0_FILENAME));

        String[] stringQ0 = buffQ0.readLine().split("\\s");
        double[] q0 = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];
        for (int iPair = 0; iPair < (GeneticCode.CODON_STATES * GeneticCode.CODON_STATES); iPair++) {
            q0[iPair] = Double.parseDouble(stringQ0[iPair]);
        }
        buffQ0.close();

        double mutsAllDenom = 0, mutsNonSynDenom = 0, subsAllDenom = 0, subsNonSynDenom = 0;
        
        String lineQ;
        
        while ((lineQ = buffQ.readLine()) != null) {
        
            String linePi = buffPi.readLine();
            String[] stringQ = lineQ.split("\\s");
            String[] stringPi = linePi.split("\\s");

            double mutsAllSiteSum = 0, mutsNonSynSiteSum = 0, subsAllSiteSum = 0, subsNonSynSiteSum = 0;
            
            for (int iPair = 0; iPair < (GeneticCode.CODON_STATES * GeneticCode.CODON_STATES); iPair++) {
                int iCodon = (iPair - iPair % GeneticCode.CODON_STATES) / GeneticCode.CODON_STATES;
                char iRes = aaCode[iCodon];
                char jRes = aaCode[iPair % GeneticCode.CODON_STATES];
                double piValue = Double.parseDouble(stringPi[iCodon]);
                double qValue = Double.parseDouble(stringQ[iPair]);

                boolean isCodonChange = iCodon != (iPair % GeneticCode.CODON_STATES);
                boolean isNonSynChange = iRes != jRes;

                if (isNonSynChange) {
                    mutsNonSynSiteSum += piValue * q0[iPair];
                    subsNonSynSiteSum += piValue * qValue;
                    mutsAllSiteSum += piValue * q0[iPair];
                    subsAllSiteSum += piValue * qValue;
                } else if (isCodonChange) {
                    mutsAllSiteSum += piValue * q0[iPair];
                    subsAllSiteSum += piValue * qValue;
                }
            }

            mutsAllDenom += mutsAllSiteSum;
            mutsNonSynDenom += mutsNonSynSiteSum;
            subsAllDenom += subsAllSiteSum;
            subsNonSynDenom += subsNonSynSiteSum;
            
        }
        buffQ.close();
        buffPi.close();

        buffQ = new BufferedReader(new FileReader(ROOT_DIR + Constants.QS_FILENAME));
		buffPi = new BufferedReader(new FileReader(ROOT_DIR + Constants.PI_FILENAME));

        double propsig = 0D;
        double propnonsig = 0D;
        int lines = 0;
        String lineS;
        while ((lineS = buffS.readLine()) != null) {
            lines++;
            lineQ = buffQ.readLine();
            String linePi = buffPi.readLine();
            String[] stringS = lineS.split("\\s");
            String[] stringQ = lineQ.split("\\s");
            String[] stringPi = linePi.split("\\s");

            for (int iPair = 0; iPair < (GeneticCode.CODON_STATES * GeneticCode.CODON_STATES); iPair++) {
                int iCodon = (iPair - iPair % GeneticCode.CODON_STATES) / GeneticCode.CODON_STATES;
                char iRes = aaCode[iCodon];
                char jRes = aaCode[iPair % GeneticCode.CODON_STATES];
                double piValue = Double.parseDouble(stringPi[iCodon]);
                double qValue = Double.parseDouble(stringQ[iPair]);
                boolean isCodonChange = iCodon != (iPair % GeneticCode.CODON_STATES);
                boolean isNonSynChange = iRes != jRes;

                double deltaS = Double.parseDouble(stringS[iPair]);
                deltaS = Math.max(deltaS, -S_LIMIT);
                deltaS = Math.min(deltaS, S_LIMIT);
                double val = deltaS - low;

                int mutsBin = (int) (MUTS_BINS * (val / (hi - low)));
                int subsBin = (int) (SUBS_BINS * (val / (hi - low)));

                if (isNonSynChange) {
                    mutsNonSyn[mutsBin] += (piValue * q0[iPair]) / mutsNonSynDenom;
                    subsNonSyn[subsBin] += (piValue * qValue) / subsNonSynDenom;
                    mutsAll[mutsBin] += (piValue * q0[iPair]) / mutsAllDenom;
                    subsAll[subsBin] += (piValue * qValue) / subsAllDenom;
                } else if (isCodonChange) {
                    mutsAll[mutsBin] += (piValue * q0[iPair]) / mutsAllDenom;
                    subsAll[subsBin] += (piValue * qValue) / subsAllDenom;
                }
            }
        }

        BufferedWriter w1 = Files.newWriter(new File("muts.all.txt"), Charsets.US_ASCII);
        BufferedWriter w2 = Files.newWriter(new File("muts.nonsyn.txt"), Charsets.US_ASCII);        
        for (int i = 0 ; i < mutsAll.length; i++) {
            w1.write(String.format("%s\n", mutsAll[i]));
            w2.write(String.format("%s\n", mutsNonSyn[i]));
            
        }
        w1.close();
        w2.close();
        
        BufferedWriter w3 = Files.newWriter(new File("subs.all.txt"), Charsets.US_ASCII);
        BufferedWriter w4 = Files.newWriter(new File("subs.nonsyn.txt"), Charsets.US_ASCII);
        for (int i = 0 ; i < subsAll.length; i++) {
            w4.write(String.format("%s\n", subsNonSyn[i]));
            w3.write(String.format("%s\n", subsAll[i]));
        }
        w3.close();
        w4.close();
    }

}
