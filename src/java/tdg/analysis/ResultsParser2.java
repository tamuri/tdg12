package tdg.analysis;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.google.common.base.Charsets;
import com.google.common.base.Function;
import com.google.common.base.Functions;
import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import com.google.common.primitives.Doubles;
import tdg.cli.DoubleArrayConverter;
import tdg.cli.DoubleConverter;
import tdg.cli.GeneticCodeConverter;
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
   /* double tau = 1.25000000000000e-02;
    double kappa = 7.06084250000000e+00;
    double[] pi = new double[]{1.85704375000000e-01, 2.73384375000000e-01, 4.78586875000000e-01, 6.23243750000000e-02};
    double mu = 2.41435441700000e+00;
    double gamma = 0;*/


    // PB2
    /*double tau = 1.25010000000000e-02;
    double kappa = 7.86404250000000e+00;
    double[] pi = new double[]{2.36114375000000e-01, 1.95774375000000e-01, 3.65944375000000e-01, 2.02166875000000e-01};
    double mu = 3.13262666700000e+00;
    double gamma = 0;

    final TDGGlobals tdgGlobals = new TDGGlobals(tau, kappa, pi, mu, gamma);
    final FileWriter outS;
    final FileWriter outPiS;
    final FileWriter outQS;
    static String prefix = "/Users/atamuri/Documents/work/mitochondria/110329_TdG_PB2_FullMaxEnt/results.homog/";*/

    // Mit Sim 2
/*

    double tau = 4.37500000000000e-02;
    double kappa = 5.32475500000000e+00;
    double[] pi = new double[]{0.18383875,0.27358875,0.47528375,0.06728875};
    double mu = 2.68226099300000e+00;
    double gamma = 0;

*/

    TDGGlobals tdgGlobals;
    FileWriter outS;
    FileWriter outPiS;
    FileWriter outQS;
    //static String prefix = "/Users/atamuri/Documents/work/mitochondria/110506_MitSim/1/";
    //static String prefix = "/Users/atamuri/Documents/work/mitochondria/110328_Mit_PenalisedLnL_MaxEnt/";
    String prefix = "./";


    public static void main(String[] args) throws Exception {
       //GeneticCode.initialise(GeneticCode.VERTEBRATE_MITOCHONDRIAL_CODE);
       // GeneticCode.initialise(GeneticCode.STANDARD_CODE);
        ROptions o = new ROptions();
        new JCommander(o, args);

        ResultsParser2 rp = new ResultsParser2();
        rp.tdgGlobals = new TDGGlobals(o.tau, o.kappa, o.pi, o.mu, o.gamma);
//        System.setOut(new PrintStream(new File("/Users/atamuri/Documents/work/mitochondria/110318_TdG_Mit_ApproxVsFull/out.data")));
 //       rp.run("/Users/atamuri/Documents/work/mitochondria/110318_TdG_Mit_ApproxVsFull/full.fitness.ordered.txt");
//        System.setOut(new PrintStream(new File("/Users/atamuri/Documents/work/mitochondria/110328_Mit_PenalisedLnL_MaxEnt/out.data")));
  //      rp.run("/Users/atamuri/Documents/work/mitochondria/110328_Mit_PenalisedLnL_MaxEnt/fitness.sorted.txt");

        rp.run(rp.prefix + "fitness.sorted.txt");
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
        List<Double> Q0 = Lists.newArrayList();
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                int aa_from = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(i);
                int aa_to = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(j);

                if (aa_from < 0 || aa_to < 0) continue;


                if (i == j) Q0.add(-1.0);
                if (i != j) Q0.add(tdgGlobals.getNeutralMutationRate(i, j) * tdgGlobals.getNu());
            }
        }
        FileWriter outQ0 = new FileWriter(new File(prefix + "Q0.txt"));
        outQ0.write(Joiner.on(' ').join(Q0));
        outQ0.close();


        // System.exit(0);
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

            // exact method
            List<Integer> aminoAcids = Lists.newArrayList(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19);
            TDGCodonModel tdg = new TDGCodonModel(tdgGlobals, new Fitness(Doubles.toArray(fitnesses), false), aminoAcids);

            // approx method
            /*List<Integer> aa = Lists.newArrayList();
            List<Double> ff = Lists.newArrayList();
            int pos = 0;
            for (double f : fitnesses) {
                if (f != -21) {
                    // this amino acid is observed
                    aa.add(pos);
                    ff.add(f);
                    pos++;
                }
            }
            TDGCodonModel tdg = new TDGCodonModel(tdgGlobals, new Fitness(Doubles.toArray(ff), false), aa);
*/


            tdg.updateModel();

            double[] S = tdg.getS();
            double[] QS = tdg.getFullQ();
            double[] PiS = tdg.getCodonFrequencies();

            // without STOP codons:
            /*List<Double> S_nostop = Lists.newArrayList();
            List<Double> QS_nostop = Lists.newArrayList();
            List<Double> PiS_nostop = Lists.newArrayList();

            for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
                if (GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(i) >= 0) {
                    PiS_nostop.add(PiS[i]);
                    for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                        if (GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(j) >= 0) {
                            S_nostop.add(S[i * GeneticCode.CODON_STATES + j]);
                            QS_nostop.add(QS[i * GeneticCode.CODON_STATES + j]);
                        }
                    }
                }
            }
            outS.write(String.format("%s\n", Joiner.on(' ').join(S_nostop)));
            outPiS.write(String.format("%s\n", Joiner.on(' ').join(PiS_nostop)));
            outQS.write(String.format("%s\n", Joiner.on(' ').join(QS_nostop)));
*/
            outS.write(String.format("%s\n", Doubles.join(" ", S)));
            outPiS.write(String.format("%s\n", Doubles.join(" ", PiS)));
            outQS.write(String.format("%s\n", Doubles.join(" ", QS)));




            return true;
        }

        @Override
        public Object getResult() {
            return null;
        }
    }
}
    class ROptions {
    @Parameter(names = "-tau", description = "Rate of multiple substitutions.", converter = DoubleConverter.class, required = true)
    public double tau;

    @Parameter(names = "-kappa", description = "Transition/transversion bias.", converter = DoubleConverter.class, required = true)
    public double kappa;

    @Parameter(names = "-pi", description = "Comma-separated base nucleotide frequencies (T,C,A,G).", converter = DoubleArrayConverter.class, required = true)
    public double[] pi;

    @Parameter(names = "-mu", description = "Branch/rate scaling factor.", converter = DoubleConverter.class, required = true)
    public double mu;

    double gamma = 0;

    @Parameter(names = "-gc", description = "The genetic code translation to use (standard or vertebrate_mit).", required = true, converter = GeneticCodeConverter.class)
    public GeneticCode geneticCode;
    }
