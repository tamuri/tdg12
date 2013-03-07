package tdg;

import com.beust.jcommander.JCommander;
import com.google.common.collect.Lists;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;
import org.apache.commons.math.optimization.RealConvergenceChecker;
import org.apache.commons.math.optimization.RealPointValuePair;
import org.apache.commons.math.optimization.SimpleScalarValueChecker;
import pal.alignment.Alignment;
import pal.tree.Tree;
import tdg.model.Fitness;
import tdg.model.TDGGlobals;
import tdg.utils.Pair;
import tdg.utils.PhyloUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.charset.Charset;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;

/**
 * The main class for model parameter estimation. This is the class called by end-users.
 *
 * TODO: Implement checkpointing, to read current globals, tree and fitness parameters
 */
public class Estimator {

    EstimatorOptions options = new EstimatorOptions();


    public static void main(String[] args) {
        Estimator e = new Estimator(args);
        e.run();
    }

    public Estimator(String... args) {
        JCommander jc = new JCommander(options);
        jc.parse(args);
    }

    private void run() {


        Tree tree = PhyloUtils.readTree(options.tree);
        Alignment alignment = PhyloUtils.readAlignment(options.alignment);

        if (!PhyloUtils.isTreeAndAlignmentValid(tree, alignment)) {
            throw new RuntimeException("ERROR: tree and alignment do not have the same taxa.");
        }

        MatrixArrayPool.treeSize = tree.getInternalNodeCount();

        Runner runner;

        if (options.distributed) {
            List<String> slaves = Lists.newArrayList();
            for (String s : options.hosts) slaves.add("http://localhost:" + s + "/service");
            runner = new DistributedRunner(alignment, slaves);
        } else {
            // Object pooling: http://code.google.com/p/furious-objectpool/ or http://commons.apache.org/proper/commons-pool/
            runner = new MultiThreadedRunner(alignment, options.threads);
        }


        // Step 0 - Set the initial parameters

        // Get rid of the current branch lengths and set to a sensible initial value
        PhyloUtils.setAllBranchLengths(tree, Constants.INITIAL_BRANCH_LENGTH);

        // Default constructor sets: -tau 0.01 -kappa 2.0 -pi 0.25,0.25,0.25 -mu 1.0
        TDGGlobals globals = new TDGGlobals();

        // Use the mutational matrix only for the first iteration (all 20 amino acids have F = 0)
        Fitness intialFitness = Fitness.getMutationOnlyFitness();
        FitnessStore fitnessStore = new FitnessStore(alignment.getSiteCount() / 3);
        for (int i = 1; i <= alignment.getSiteCount() / 3; i++) {
            fitnessStore.setFitness(i, intialFitness);
        }

        RealConvergenceChecker convergenceChecker = new SimpleScalarValueChecker(-1, 1e-4);
        RealPointValuePair previous;
        RealPointValuePair current = new RealPointValuePair(new double[]{}, Double.NEGATIVE_INFINITY);

        int iteration = 0;
        boolean converged = false;

        while (!converged) {
            long start = System.currentTimeMillis();
            Date starttime = new Date();

            iteration++;

            // Step 1 - optimise the site-invariant mutational parameters, get new TDGGlobals
            Pair<Double, TDGGlobals> globalsOptResult = runner.optimiseMutationModel(tree, globals, fitnessStore);
            System.out.printf("%s - Mutational parameters optimisation: %s ( %s )\n", iteration, globalsOptResult.first, globalsOptResult.second.toString());
            globals = globalsOptResult.second;

            // Step 2 - optimise the branch lengths, get updated tree
            Pair<Double, Tree> treeOptResult = runner.optimiseBranchLengths(tree, globals, fitnessStore);
            System.out.printf("%s - Branch length optimisation: %s ( Total tree length: %s )\n", iteration, treeOptResult.first, PhyloUtils.getTotalTreeLength(treeOptResult.second));
            tree = treeOptResult.second;

            // Step 3 - optimise the fitness parameters
            double fitnessOptResult = runner.optimiseFitness(tree, globals, fitnessStore);
            System.out.printf("%s - Fitness optimisation: %s\n", iteration, fitnessOptResult);

            if (iteration == 1) {
                current = new RealPointValuePair(new double[]{}, fitnessOptResult);
            } else {
                previous = current;
                current = new RealPointValuePair(new double[]{}, fitnessOptResult);
                converged = convergenceChecker.converged(iteration, previous, current);
            }

            long end = System.currentTimeMillis();

            writeResults(starttime, end - start, iteration, tree, globals, fitnessStore);
        }

        System.out.println();
        System.out.printf("%s\n", tree.toString());



        runner.close();

    }

    private void writeResults(Date date, long millis, int iteration, Tree tree, TDGGlobals globals, FitnessStore fitnessStore) {

        String filename = String.format("%03d_optima.txt", iteration);

        try {
            BufferedWriter bw = Files.newWriter(new File(filename), Charset.defaultCharset());

            bw.write(new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(date));
            bw.newLine();

            bw.write(globals.toString());
            bw.newLine();
            bw.write(tree.toString());
            bw.newLine();
            for (Fitness f : fitnessStore) {
                bw.write(Doubles.join(",", f.get()));
                bw.newLine();
            }
            bw.newLine();

            bw.write((millis / 1000.0) + " secs.");
            bw.newLine();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

}

