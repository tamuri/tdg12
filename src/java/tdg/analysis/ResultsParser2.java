package tdg.analysis;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.google.common.base.Charsets;
import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableList;
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
import tdg.utils.Functions;
import tdg.utils.GeneticCode;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

public class ResultsParser2 {
    TDGGlobals tdgGlobals;
    String path;
    ROptions o;

    FileWriter outS;
    FileWriter outPiS;
    FileWriter outQS, outPiAA;
    List<Integer> aminoAcids = ImmutableList.copyOf(Lists.<Integer>newArrayList(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)) ;

    public static void main(String[] args) throws Exception {
        ROptions o = new ROptions();
        new JCommander(o, args);

        ResultsParser2 rp = new ResultsParser2(o, new TDGGlobals(o.tau, o.kappa, o.pi, o.mu, o.gamma), "fitness.sorted.txt");
        rp.run();
   }

    public ResultsParser2(ROptions o, TDGGlobals globals, String filePath) throws Exception {
        this.o = o;
        this.tdgGlobals = globals;
        this.path = filePath;
    }

    private void run() throws Exception {
        // Neutral Q
        List<Double> Q0 = Lists.newArrayListWithCapacity(GeneticCode.CODON_STATES * GeneticCode.CODON_STATES);
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                if (i == j) Q0.add(-1.0);
                else Q0.add(tdgGlobals.getNeutralMutationRate(i, j) * tdgGlobals.getNu());
            }
        }

        FileWriter outQ0 = new FileWriter(new File("Q0.txt"));
        outQ0.write(Joiner.on(' ').join(Q0));
        outQ0.close();

        // S, Q with selection and codon pi files
        outS = new FileWriter(new File("S.txt"));
        outPiS = new FileWriter(new File("PiS.txt"));
        outQS = new FileWriter(new File("QS.txt"));
        outPiAA = new FileWriter(new File("PiAA.txt"));

        Files.readLines(new File(this.path), Charsets.UTF_8, new SWriter());

        outS.close();
        outPiS.close();
        outQS.close();
        outPiAA.close();
        
    }

    class SWriter implements LineProcessor<Object> {
        @Override
        public boolean processLine(String line) throws IOException {

            List<Double> fitnesses = Lists.transform(Arrays.asList(line.split(" ")), Functions.stringToDouble());

            TDGCodonModel tdg;

            if (o.approx) {
                List<Integer> aa = Lists.newArrayList();
                List<Double> ff = Lists.newArrayList();
                int pos = 0;
                for (double f : fitnesses) {
                    if (!Double.isInfinite(f)) {
                        // this amino acid is observed
                        aa.add(pos);
                        ff.add(f);
                        pos++;
                    }
                }
                tdg = new TDGCodonModel(tdgGlobals, new Fitness(Doubles.toArray(ff), false), aa);
            } else {
                tdg = new TDGCodonModel(tdgGlobals, new Fitness(Doubles.toArray(fitnesses), false), aminoAcids);
            }

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
            outS.write(String.format("%s\n", Doubles.join(" ", S)));
            outPiS.write(String.format("%s\n", Doubles.join(" ", PiS)));
            outQS.write(String.format("%s\n", Doubles.join(" ", QS)));
            outPiAA.write(String.format("%s\n", Doubles.join(" ", tdg.getAminoAcidFrequencies())));
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

    @Parameter(names = "-approx", description = "Use the approximate method to optimise the likelihood")
    public boolean approx = false;

}
    