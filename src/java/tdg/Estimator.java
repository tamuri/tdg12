package tdg;

import com.beust.jcommander.JCommander;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;
import org.apache.commons.math.optimization.RealConvergenceChecker;
import org.apache.commons.math.optimization.RealPointValuePair;
import org.apache.commons.math.optimization.SimpleScalarValueChecker;
import pal.alignment.Alignment;
import pal.tree.Tree;
import tdg.model.Fitness;
import tdg.model.TDGGlobals;
import tdg.utils.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.sql.Timestamp;
import java.util.List;

/**
 * The main class for model parameter estimation. This is the class called by end-users.
 */
public class Estimator {
    EstimatorOptions options = new EstimatorOptions();

    public static void main(String[] args) {
        Estimator e = new Estimator(args);
        e.run();
    }

    public Estimator(String... args) {
        JCommander jc = new JCommander(options);
        try {
            jc.parse(args);
        } catch (Exception e) {
            System.out.printf("Error: %s\n", e.getMessage());
            jc.setProgramName("java -cp tdg12.jar tdg.Estimator");
            jc.usage();
        }
    }

    private void run() {

        CoreUtils.msg("tdg.Estimator started.\n");

        Tree tree = PhyloUtils.readTree(options.tree);
        Alignment alignment = PhyloUtils.readAlignment(options.alignment);

        if (!PhyloUtils.isTreeAndAlignmentValid(tree, alignment)) {
            throw new RuntimeException("ERROR: tree and alignment do not have the same taxa.");
        }

        // TODO: is there a better place to put this?
        MatrixArrayPool.treeSize = tree.getInternalNodeCount();

        Runner runner = getRunner(alignment);

        // Step 0 - Set the initial parameters
        TDGGlobals globals;
        FitnessStore fitnessStore;
        Triple<TDGGlobals, FitnessStore, Tree> checkpoint = null;

        // If a checkpoint file has been supplied, try to load it
        if (options.checkpointFile != null) checkpoint = loadCheckpoint(options.checkpointFile, alignment);

        // Either we're starting estimation from scratch *or* checkpoint did not load successfully
        if (checkpoint == null) {
            globals = new TDGGlobals(); // Default initial parameters

            // Use the mutational matrix only for the first iteration (all 20 amino acids have F = 0)
            Fitness initialFitness = Fitness.getMutationOnlyFitness();
            fitnessStore = new FitnessStore(alignment.getSiteCount() / 3);
            for (int i = 1; i <= alignment.getSiteCount() / 3; i++) {
                fitnessStore.setFitness(i, initialFitness);
            }

            // Get rid of the current branch lengths and set to a sensible initial value
            PhyloUtils.setAllBranchLengths(tree, Constants.INITIAL_BRANCH_LENGTH);
        } else {
            globals = checkpoint.first;
            fitnessStore = checkpoint.second;
            tree = checkpoint.third;
        }

        RealConvergenceChecker convergenceChecker = new SimpleScalarValueChecker(-1, Constants.CONVERGENCE_TOL);
        RealPointValuePair lastOptima;
        RealPointValuePair thisOptima = new RealPointValuePair(new double[]{}, Double.NEGATIVE_INFINITY);

        int iteration = 0;
        boolean converged = false;

        while (!converged) {
            long start = System.currentTimeMillis();

            iteration++;

            // Step 1 - optimise the site-invariant mutational parameters, get new TDGGlobals
            Pair<Double, TDGGlobals> globalsOptResult = runner.optimiseMutationModel(tree, globals, fitnessStore, options.prior);
            CoreUtils.msg("%s - Mutation matrix optima: %s ( %s )\n", iteration, globalsOptResult.first, globalsOptResult.second.toString());
            globals = globalsOptResult.second;

            // Step 2 - optimise the branch lengths, get updated tree
            Pair<Double, Tree> treeOptResult = runner.optimiseBranchLengths(tree, globals, fitnessStore, options.prior);
            CoreUtils.msg("%s - Branch length optima: %s ( Total tree length: %s )\n", iteration, treeOptResult.first, PhyloUtils.getTotalTreeLength(treeOptResult.second));
            tree = treeOptResult.second;

            // Step 3 - optimise the fitness parameters
            double iterationLnL = runner.optimiseFitness(tree, globals, fitnessStore, options.prior);
            CoreUtils.msg("%s - Fitness optima: %s\n", iteration, iterationLnL);

            if (iteration == 1) {
                thisOptima = new RealPointValuePair(new double[]{}, iterationLnL);
            } else {
                lastOptima = thisOptima;
                thisOptima = new RealPointValuePair(new double[]{}, iterationLnL);
                converged = convergenceChecker.converged(iteration, lastOptima, thisOptima);
            }

            long end = System.currentTimeMillis();

            writeResults(start, end, iteration, alignment,
                    globals, globalsOptResult.first,
                    tree, treeOptResult.first,
                    fitnessStore, iterationLnL);
        }

        System.out.println();

        CoreUtils.msg("tdg.Estimator converged.\n");

        runner.close();
    }

    private Runner getRunner(Alignment alignment) {
        Runner runner;
        if (options.distributed) {
            // Distributed RPC alternatives: hessian, Protocolbuffer-rpc-pro, finagle,
            // http://code.google.com/p/missian/
            List<String> slaves;
            try {
                slaves = Files.readLines(new File(options.hostsFile), Charset.defaultCharset());
            } catch (IOException e) {
                throw new RuntimeException(e);
            }

            for (int i = 0; i < slaves.size(); i++) slaves.set(i, "http://" + slaves.get(i) + "/service");

            runner = new DistributedRunner(alignment, slaves);

        } else {
            // Object pooling: http://code.google.com/p/furious-objectpool/ or http://commons.apache.org/proper/commons-pool/
            runner = new MultiThreadedRunner(alignment, options.threads);
        }
        return runner;
    }

    private Triple<TDGGlobals, FitnessStore, Tree> loadCheckpoint(String checkpointFile, final Alignment alignment) {
        try {
            return Files.readLines(new File(checkpointFile), Charset.defaultCharset(), new CheckpointLineProcessor(alignment));
        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }
    }

    private void writeResults(long start, long end, int iteration,
                              Alignment alignment, TDGGlobals globals, double globalsOpt,
                              Tree tree, double treeOpt,
                              FitnessStore fitnessStore, double fitnessOpt) {

        String filename = String.format("%03d_optima.txt", iteration);

        try {
            BufferedWriter bw = Files.newWriter(new File(filename), Charset.defaultCharset());

            bw.write(new Timestamp(start).toString());
            bw.newLine();

            bw.write(globals.toString());
            bw.newLine();
            bw.write(tree.toString());
            bw.newLine();

            int site = 1;
            for (Fitness f : fitnessStore) {
                List<Integer> residues = PhyloUtils.getDistinctAminoAcids(PhyloUtils.getCleanedCodons(alignment, site).values());
                double[] fitnesses = new double[20];
                for (int i = 0; i < GeneticCode.AMINO_ACID_STATES; i++) fitnesses[i] = f.get()[residues.indexOf(i)];
                bw.write(Doubles.join(",", fitnesses));
                bw.newLine();
                site++;
            }

            bw.newLine();
            bw.write(String.format("lnL - %s %s %s\n", globalsOpt, treeOpt, fitnessOpt));
            bw.write(((end - start) / 1000.0) + " secs.");
            bw.newLine();

            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

}

