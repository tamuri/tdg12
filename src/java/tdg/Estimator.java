package tdg;

import com.beust.jcommander.JCommander;
import org.apache.commons.math.optimization.RealConvergenceChecker;
import org.apache.commons.math.optimization.RealPointValuePair;
import org.apache.commons.math.optimization.SimpleScalarValueChecker;
import pal.alignment.Alignment;
import pal.tree.Tree;
import tdg.model.Fitness;
import tdg.model.TDGGlobals;
import tdg.utils.Pair;
import tdg.utils.PhyloUtils;

/**
 * The main class for model parameter estimation. This is the class called by end-users.
 *
 * TODO: Implement checkpointing, to read current globals, tree and fitness parameters
 */
public class Estimator {
    private static final double INITIAL_BRANCH_LENGTH = 0.1;

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

        // Load a runner: SingleThreadRunner, MultiThreadRunner(threads), DistributedRunner(slaves)

        Runner runner = new MultiThreadedRunner(options.threads);

        if (options.runner.equals("distributed")) {
            // DistributedRunner(list of slaves)

            // Hessian seems like a good choice. See http://hessian.caucho.com/#HessianImplementationsDownload
        } else if (options.threads > 1) {
            // Object pooling: http://code.google.com/p/furious-objectpool/ or http://commons.apache.org/proper/commons-pool/
            runner = new MultiThreadedRunner(options.threads);
        } else {
            // SingleThreadedRunner
        }

        // Logic goes here, each time called "Runner" (we're not interested in how it runs)


        Tree tree = PhyloUtils.readTree(options.tree);


        Alignment alignment = PhyloUtils.readAlignment(options.alignment);

        if (!PhyloUtils.isTreeAndAlignmentValid(tree, alignment)) {
            throw new RuntimeException("ERROR: tree and alignment do not have the same taxa.");
        }

        // Object pooling: http://code.google.com/p/furious-objectpool/ or http://commons.apache.org/proper/commons-pool/

        /* MdR's routine for parameter estimation:
        (0) Set inital branch lengths to some reasonable value (say 0.1),
            set all F=0,
            and mutational parameters to some starting values as you suggest.
        (1) Optimise mutation model and mu (i.e. multiply all branches by mu)
        (2) Optimise branch lengths one by one.
        (3) Optimise fitnesses.
        (4) Repeat 1, 2 and 3 until convergence.
         */

        // Step 0
        PhyloUtils.setAllBranchLengths(tree, INITIAL_BRANCH_LENGTH);

        TDGGlobals globals = new TDGGlobals();

        FitnessStore fitnessStore = new FitnessStore(alignment.getSiteCount());
        for (int i = 1; i <= alignment.getSiteCount() / 3; i++) {
            fitnessStore.setFitness(
                    i,
                    new Fitness(
                            new double[PhyloUtils.getDistinctAminoAcids(
                                    PhyloUtils.getCleanedCodons(alignment, i).values()
                                    ).size()], true));
        }

        RealConvergenceChecker convergenceChecker = new SimpleScalarValueChecker(-1, 1e-3);
        RealPointValuePair previous;
        RealPointValuePair current = new RealPointValuePair(new double[]{}, Double.NEGATIVE_INFINITY);

        int iteration = 0;
        boolean converged = false;
        while (!converged) {
            iteration++;

            // Step 1
            Pair<Double, TDGGlobals> globalsOptResult = runner.optimiseMutationModel(tree, alignment, globals, fitnessStore);
            System.out.printf("%s - Mutational parameters optimisation: %s ( %s )\n", iteration, globalsOptResult.first, globalsOptResult.second.toString());
            globals = globalsOptResult.second;

            // Step 2
            Pair<Double, Tree> treeOptResult = runner.optimiseBranchLengths(tree, alignment, globals, fitnessStore);
            System.out.printf("%s - Branch length optimisation: %s ( Total tree length: %s )\n", iteration, treeOptResult.first, PhyloUtils.getTotalTreeLength(treeOptResult.second));
            tree = treeOptResult.second;

            // Step 3
            double fitnessOptResult = runner.optimiseFitness(tree, alignment, globals, fitnessStore);
            System.out.printf("%s - Fitness optimisation: %s\n", iteration, fitnessOptResult);

            if (iteration == 1) {
                current = new RealPointValuePair(new double[]{}, fitnessOptResult);
            } else {
                previous = current;
                current = new RealPointValuePair(new double[]{}, fitnessOptResult);

                converged = convergenceChecker.converged(iteration, previous, current);
            }
        }


        runner.close();



    }

}

