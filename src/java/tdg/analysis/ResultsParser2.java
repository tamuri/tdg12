package tdg.analysis;

import com.google.common.base.Charsets;
import com.google.common.base.Function;
import com.google.common.base.Functions;
import com.google.common.collect.Lists;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import com.google.common.primitives.Doubles;
import tdg.models.TDGCodonModel;
import tdg.models.TDGGlobals;
import tdg.models.parameters.Fitness;
import tdg.utils.GeneticCode;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;

public class ResultsParser2 {
    // Mit
/*    double tau = 1.25000000000000e-02;
    double kappa = 7.06084250000000e+00;
    double[] pi = new double[]{1.85704375000000e-01, 2.73384375000000e-01, 4.78586875000000e-01, 6.23243750000000e-02};
    double mu = 2.41435441700000e+00;
    double gamma = 0;
*/

    // PB2
    double tau = 1.25010000000000e-02;
    double kappa = 7.86404250000000e+00;
    double[] pi = new double[]{2.36114375000000e-01, 1.95774375000000e-01, 3.65944375000000e-01, 2.02166875000000e-01};
    double mu = 3.13262666700000e+00;
    double gamma = 0;

    final TDGGlobals tdgGlobals = new TDGGlobals(tau, kappa, pi, mu, gamma);
    final FileWriter outS;
    final FileWriter outPiS;
    final FileWriter outQS;
    static String prefix = "/Users/atamuri/Documents/work/mitochondria/110329_TdG_PB2_FullMaxEnt/results.homog/";


    public static void main(String[] args) throws Exception {
//        GeneticCode.initialise(GeneticCode.VERTEBRATE_MITOCHONDRIAL_CODE);
        GeneticCode.initialise(GeneticCode.STANDARD_CODE);
        ResultsParser2 rp = new ResultsParser2();
//        System.setOut(new PrintStream(new File("/Users/atamuri/Documents/work/mitochondria/110318_TdG_Mit_ApproxVsFull/out.data")));
 //       rp.run("/Users/atamuri/Documents/work/mitochondria/110318_TdG_Mit_ApproxVsFull/full.fitness.ordered.txt");
//        System.setOut(new PrintStream(new File("/Users/atamuri/Documents/work/mitochondria/110328_Mit_PenalisedLnL_MaxEnt/out.data")));
  //      rp.run("/Users/atamuri/Documents/work/mitochondria/110328_Mit_PenalisedLnL_MaxEnt/fitness.sorted.txt");

        rp.run(prefix + "fitness.sorted.txt");
        rp.close();
       // System.setOut(new PrintStream(new File("/Users/atamuri/Documents/work/mitochondria/110329_TdG_PB2_FullMaxEnt/results.homog/out.data")));
        //rp.run("/Users/atamuri/Documents/work/mitochondria/110329_TdG_PB2_FullMaxEnt/results.homog/fitness.sorted.txt");
   }

    private void close() throws Exception {
        outS.close();
        outPiS.close();
        outQS.close();
    }

    public ResultsParser2() throws Exception {
        outS = new FileWriter(new File(prefix + "S.txt"));
        outPiS = new FileWriter(new File(prefix + "PiS.txt"));
        outQS = new FileWriter(new File(prefix + "QS.txt"));
    }

    private void run(String fitnessFilePath) throws Exception {
        Files.readLines(new File(fitnessFilePath), Charsets.UTF_8, new SWriter());
    }

    class SWriter implements LineProcessor<Object> {
        @Override
        public boolean processLine(String line) throws IOException {
            List<Double> fitnesses = Lists.transform(Arrays.asList(line.split(" ")), new Function<String,Double>() {
                @Override
                public Double apply(String s) {
                    return Double.parseDouble(s);
                }
            });

            List<Integer> aminoAcids = Lists.newArrayList(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19);
            TDGCodonModel tdg = new TDGCodonModel(tdgGlobals, new Fitness(Doubles.toArray(fitnesses), false), aminoAcids);
            tdg.updateModel();

            outS.write(String.format("%s\n", Doubles.join(" ", tdg.getS())));
            outPiS.write(String.format("%s\n", Doubles.join(" ", tdg.getCodonFrequencies())));
            outQS.write(String.format("%s\n", Doubles.join(" ", tdg.getFullQ())));
//            System.out.printf("%s\n", Doubles.join(" ", tdg.getS()));
  //          System.out.printf("%s\n", Doubles.join(" ", tdg.getCodonFrequencies()));
   //         System.out.printf("%s\n", Doubles.join(" ", tdg.getFullQ()));

            return true;
        }

        @Override
        public Object getResult() {
            return null;
        }
    }
}
