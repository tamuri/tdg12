package tdg;

import com.google.common.collect.Lists;
import com.google.common.io.LineProcessor;
import pal.alignment.Alignment;
import pal.tree.ReadTree;
import pal.tree.Tree;
import pal.tree.TreeParseException;
import tdg.model.Fitness;
import tdg.model.TDGGlobals;
import tdg.utils.PhyloUtils;
import tdg.utils.Triple;

import java.io.IOException;
import java.io.PushbackReader;
import java.io.StringReader;
import java.sql.Timestamp;
import java.util.List;

/**
 * Author: Asif Tamuri (atamuri@ebi.ac.uk)
 * Date: 21/03/2013 22:13
 */
public class CheckpointLineProcessor implements LineProcessor<Triple<TDGGlobals, FitnessStore, Tree>> {
    private TDGGlobals globals;
    private Tree tree;
    private List<Fitness> fitnessList = Lists.newArrayList();
    private Alignment alignment;

    public CheckpointLineProcessor(Alignment alignment) {
        this.alignment = alignment;
    }

    @Override
    public boolean processLine(String line) throws IOException {

        String[] parts;

        if (line.startsWith("TDGGlobals")) {
            // TDGGlobals{ -tau 0.1056066601361324 -kappa 3.5149103405053466 -pi 0.2702656644776448,0.2524978838686726,0.3926835680493317,0.08455288360435098 -mu 3.6173836635167165 (1/nu=0.24262321603322343)}
            parts = line.split(" ");
            double tau = Double.parseDouble(parts[2]);
            double kappa = Double.parseDouble(parts[4]);

            String[] piParts = parts[6].split(",");
            double[] pi = new double[4];
            for (int i = 0; i < pi.length; i++) {
                pi[i] = Double.parseDouble(piParts[i]);
            }

            double mu = Double.parseDouble(parts[8]);

            this.globals = new TDGGlobals(tau, kappa, pi, mu, 0);
            System.out.printf("%s - tdg.Estimator loaded checkpoint (%s)\n", new Timestamp(System.currentTimeMillis()), this.globals.toString());
        } else if (line.startsWith("(")) {

            try {
                this.tree = new ReadTree(new PushbackReader(new StringReader(line)));
            } catch (TreeParseException e) {
                e.printStackTrace();
                throw new IOException("ERROR: Could not read checkpoint tree.");
            }

            // check the tree will work with the alignment
            if (!PhyloUtils.isTreeAndAlignmentValid(this.tree, alignment)) {
                throw new IOException("ERROR: Checkpoint tree does not match alignment.");
            }

            System.out.printf("%s - tdg.Estimator loaded checkpoint (tree length %s)\n", new Timestamp(System.currentTimeMillis()), PhyloUtils.getTotalTreeLength(this.tree));
        } else if (globals != null && tree != null && (parts = line.split(",")).length == 20) {
            // If we've successfully loaded mutational parameters and tree, and this line contains 20 comma-separated values

            // Re-order the fitnesses when we read them in
            int site = fitnessList.size() + 1;
            List<Integer> residues = PhyloUtils.getDistinctAminoAcids(PhyloUtils.getCleanedCodons(alignment, site).values());

            double[] f = new double[20];
            for (int i = 0; i < f.length; i++) {
                f[i] = Double.parseDouble(parts[residues.get(i)]);
            }

            fitnessList.add(new Fitness(f, true));
        }

        return true;
    }

    @Override
    public Triple<TDGGlobals, FitnessStore, Tree> getResult() {
        int site = 1;
        FitnessStore fitnessStore = new FitnessStore();
        for (Fitness f : fitnessList) fitnessStore.setFitness(site++, f);
        System.out.printf("%s - tdg.Estimator loaded checkpoint (%s fitness vectors)\n", new Timestamp(System.currentTimeMillis()), fitnessList.size());

        return Triple.of(globals, fitnessStore, tree);
    }
}
