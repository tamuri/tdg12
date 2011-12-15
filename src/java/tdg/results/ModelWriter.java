package tdg.results;

import com.beust.jcommander.JCommander;
import com.google.common.base.Charsets;
import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import com.google.common.primitives.Doubles;
import tdg.Constants;
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

/**
 * Takes a list of fitnesses (usually made by running FitnessExtractor) and swMutSel0 global parameters and
 * writes the following files:
 *
 * 1. Q0.txt - the neutral substitution rate matrix for each site
 * 2. QS.txt - the substitution rate matrix, with selection
 * 3. S.txt - the selection coefficient (S_ij) matrix
 * 4. PiS.txt - the codon frequencies
 * 5. PiAA.txt - the amino acid frequencies
 *
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 * @see FitnessExtractor
 */
public class ModelWriter {
    TDGGlobals tdgGlobals;
    String path;
    Options o;

    FileWriter outS;
    FileWriter outPiS;
    FileWriter outQS, outPiAA;
    List<Integer> aminoAcids = ImmutableList.copyOf(Lists.<Integer>newArrayList(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)) ;

    public static void main(String[] args) throws Exception {
        Options o = new Options();
        new JCommander(o, args);
        
        ModelWriter rp = new ModelWriter(o, new TDGGlobals(o.tau, o.kappa, o.pi, o.mu, o.gamma), Constants.F_FILENAME);
        rp.run();
   }

    public ModelWriter(Options o, TDGGlobals globals, String filePath) throws Exception {
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

        FileWriter outQ0 = new FileWriter(new File(Constants.Q0_FILENAME));
        outQ0.write(Joiner.on(' ').join(Q0));
        outQ0.close();

        // S, Q with selection and codon pi files
        outS = new FileWriter(new File(Constants.S_FILENAME));
        outPiS = new FileWriter(new File(Constants.PI_FILENAME));
        outQS = new FileWriter(new File(Constants.QS_FILENAME));
        outPiAA = new FileWriter(new File(Constants.PIAA_FILENAME));

        Files.readLines(new File(this.path), Charsets.UTF_8, new FitnessProcessor());

        outS.close();
        outPiS.close();
        outQS.close();
        outPiAA.close();

    }

    class FitnessProcessor implements LineProcessor<Object> {
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

            tdg.updateModel();
            double[] S = tdg.getS();
            double[] QS = tdg.getFullQ();
            double[] PiS = tdg.getCodonFrequencies();
            outS.write(String.format("%s\n", Doubles.join(" ", S)));
            outPiS.write(String.format("%s\n", Doubles.join(" ", PiS)));
            outQS.write(String.format("%s\n", Doubles.join(" ", QS)));
            outPiAA.write(String.format("%s\n", Doubles.join(" ", tdg.getAminoAcidFrequencies())));

            return true;
        }

        @Override
        public Object getResult() {
            return null;
        }
    }

}
