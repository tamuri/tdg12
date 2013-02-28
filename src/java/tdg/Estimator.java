package tdg;

import com.beust.jcommander.JCommander;
import com.google.common.collect.Lists;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;
import pal.alignment.Alignment;
import pal.tree.Tree;
import tdg.model.Fitness;
import tdg.model.TDGGlobals;
import tdg.utils.Functions;
import tdg.utils.Pair;
import tdg.utils.PhyloUtils;
import tdg.utils.m;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.List;

/**
 * The main class for model parameter estimation. This is the class called by end-users.
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

/*

        try {
            List<String> lines = Files.readLines(new File("/Users/atamuri/Documents/2013/tdg12/f1.txt"), Charset.defaultCharset());
            int site = 1;
            for (String l : lines) {
                String[] parts = l.split(",");
                List<Double> d = Lists.transform(Lists.newArrayList(parts), Functions.stringToDouble());
                d.remove(0);
                Fitness f = new Fitness(Doubles.toArray(d), true);
                System.out.printf("%s\n", f);
                fitnessStore.setFitness(site, f);
                site++;
                // if (site > 5) break;
            }

        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        globals = new TDGGlobals(0.18384062387520766, 4.128219745551803, new double[]{0.28119571640365953, 0.24370918125910898, 0.3741515588418039, 0.10094354349542756}, 2.6629836091051446,  0);
*/


        // System.exit(0);


        for (int j : m.seqi(0, 1)) {




            // Step 1
            Pair<Double, TDGGlobals> result1 = runner.optimiseMutationModel(tree, alignment, globals, fitnessStore);
            System.out.printf("Optima: %s = %s\n", result1.second.toString(), result1.first);
            globals = result1.second;



            // Step 2
            Pair<Double, Tree> result2 = runner.optimiseBranchLengths(tree, alignment, globals, fitnessStore);
            System.out.printf("Optima: %s\n%s\n", result2.second, result2.first);
            tree = result2.second;


            // Step 3


            double result3 = runner.optimiseFitness(tree, alignment, globals, fitnessStore);
            System.out.printf("Optima: %s\n", result3);

            /*for (int i = 1; i < alignment.getSiteCount() / 3; i++) {
                fitnessStore.setFitness(i, result3.second.get(i - 1));
            }
*/



            //  System.out.println("TIME: ");
            //System.out.printf("1: %ss\t 2: %ss\t3: %ss\ttotal: %ss\n\n", (start2 - start1) / 1000, (start3 - start2) / 1000, (end - start3) / 1000, (end - start1) / 1000);

        }


        // Do any cleanup that needs to be done
        runner.close();



    }

}

