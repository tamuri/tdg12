package tdg;

import com.beust.jcommander.JCommander;
import com.google.common.collect.Lists;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import com.google.common.primitives.Doubles;
import org.apache.commons.math.optimization.RealConvergenceChecker;
import org.apache.commons.math.optimization.RealPointValuePair;
import org.apache.commons.math.optimization.SimpleScalarValueChecker;
import pal.alignment.Alignment;
import pal.tree.ReadTree;
import pal.tree.Tree;
import pal.tree.TreeParseException;
import tdg.model.Fitness;
import tdg.model.TDGGlobals;
import tdg.utils.Pair;
import tdg.utils.PhyloUtils;
import tdg.utils.Triple;

import java.io.*;
import java.nio.charset.Charset;
import java.sql.Timestamp;
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
        System.out.printf("%s - tdg.Estimator started.\n", new Timestamp(System.currentTimeMillis()));

        Tree tree = PhyloUtils.readTree(options.tree);
        Alignment alignment = PhyloUtils.readAlignment(options.alignment);

        if (!PhyloUtils.isTreeAndAlignmentValid(tree, alignment)) {
            throw new RuntimeException("ERROR: tree and alignment do not have the same taxa.");
        }

        MatrixArrayPool.treeSize = tree.getInternalNodeCount();

        Runner runner;

        if (options.distributed) {

            List<String> slaves = null;
            try {
                slaves = Files.readLines(new File(options.hostsFile), Charset.defaultCharset());
            } catch (IOException e) {
                throw new RuntimeException(e);
            }

            // for (String s : options.hosts) slaves.add("http://localhost:" + s + "/service");
            for (int i = 0; i < slaves.size(); i++) {
                slaves.set(i, "http://" + slaves.get(i) + "/service");
            }
            runner = new DistributedRunner(alignment, slaves);
        } else {
            // Object pooling: http://code.google.com/p/furious-objectpool/ or http://commons.apache.org/proper/commons-pool/
            runner = new MultiThreadedRunner(alignment, options.threads);
        }


        // Step 0 - Set the initial parameters

        // TODO: implement a type of checkpointing....use provides a xxx_optima.txt file & we continue from there

        TDGGlobals globals;
        FitnessStore fitnessStore;
        Triple<TDGGlobals, FitnessStore, Tree> checkpoint = null;

        if (options.checkpointFile != null) {
            checkpoint = loadCheckpoint(options.checkpointFile, alignment);
        }

        if (checkpoint == null) {
            globals = new TDGGlobals(); // Default constructor sets: -tau 0.01 -kappa 2.0 -pi 0.25,0.25,0.25 -mu 1.0

            // Use the mutational matrix only for the first iteration (all 20 amino acids have F = 0)
            Fitness intialFitness = Fitness.getMutationOnlyFitness();
            fitnessStore = new FitnessStore(alignment.getSiteCount() / 3);
            for (int i = 1; i <= alignment.getSiteCount() / 3; i++) {
                fitnessStore.setFitness(i, intialFitness);
            }

            // Get rid of the current branch lengths and set to a sensible initial value
            PhyloUtils.setAllBranchLengths(tree, Constants.INITIAL_BRANCH_LENGTH);
        } else {
            globals = checkpoint.first;
            fitnessStore = checkpoint.second;
            tree = checkpoint.third;
        }

        RealConvergenceChecker convergenceChecker = new SimpleScalarValueChecker(-1, 1e-5);
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
            System.out.printf("%s - %s - Mutation matrix optima: %s ( %s )\n", new Timestamp(System.currentTimeMillis()), iteration, globalsOptResult.first, globalsOptResult.second.toString());
            globals = globalsOptResult.second;

            // Step 2 - optimise the branch lengths, get updated tree
            Pair<Double, Tree> treeOptResult = runner.optimiseBranchLengths(tree, globals, fitnessStore);
            System.out.printf("%s - %s - Branch length optima: %s ( Total tree length: %s )\n", new Timestamp(System.currentTimeMillis()), iteration, treeOptResult.first, PhyloUtils.getTotalTreeLength(treeOptResult.second));
            tree = treeOptResult.second;

            // Step 3 - optimise the fitness parameters
            double fitnessOptResult = runner.optimiseFitness(tree, globals, fitnessStore);
            System.out.printf("%s - %s - Fitness optima: %s\n", new Timestamp(System.currentTimeMillis()), iteration, fitnessOptResult);

            if (iteration == 1) {
                current = new RealPointValuePair(new double[]{}, fitnessOptResult);
            } else {
                previous = current;
                current = new RealPointValuePair(new double[]{}, fitnessOptResult);
                converged = convergenceChecker.converged(iteration, previous, current);
            }

            long end = System.currentTimeMillis();

            writeResults(starttime, end - start, iteration,
                    globals, globalsOptResult.first,
                    tree, treeOptResult.first,
                    fitnessStore, fitnessOptResult);
        }

        System.out.println();

        System.out.printf("%s - tdg.Estimator converged.\n", new Timestamp(System.currentTimeMillis()));




        runner.close();

    }

    private Triple<TDGGlobals, FitnessStore, Tree> loadCheckpoint(String checkpointFile, final Alignment alignment) {


        try {

            return Files.readLines(new File(checkpointFile), Charset.defaultCharset(),
                    new LineProcessor<Triple<TDGGlobals, FitnessStore, Tree>>() {
                        private TDGGlobals globals;
                        private Tree tree;
                        private List<Fitness> fitnessList = Lists.newArrayList();

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
                            } else if (line.startsWith("0.0") && (parts = line.split(",")).length == 20) {
                                double[] f = new double[20];
                                for (int i = 0; i < f.length; i++) {
                                    f[i] = Double.parseDouble(parts[i]);
                                }
                                fitnessList.add(new Fitness(f, true));
                            }

                            return true;
                        }

                        @Override
                        public Triple<TDGGlobals, FitnessStore, Tree> getResult() {
                            int site = 1;
                            FitnessStore fitnessStore = new FitnessStore(fitnessList.size());
                            for (Fitness f : fitnessList) fitnessStore.setFitness(site++, f);
                            System.out.printf("%s - tdg.Estimator loaded checkpoint (%s fitness vectors)\n", new Timestamp(System.currentTimeMillis()), fitnessList.size());

                            return Triple.of(globals, fitnessStore, tree);
                        }
                    });

        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }

    }

    private void writeResults(Date date, long millis, int iteration,
                              TDGGlobals globals, double globalsOpt,
                              Tree tree, double treeOpt,
                              FitnessStore fitnessStore, double fitnessOpt) {

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
            bw.write(String.format("lnL - %s %s %s\n", globalsOpt, treeOpt, fitnessOpt));
            bw.write((millis / 1000.0) + " secs.");
            bw.newLine();

            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

}

