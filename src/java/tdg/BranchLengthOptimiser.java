package tdg;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import com.google.common.primitives.Doubles;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.optimization.GoalType;
import org.apache.commons.math.optimization.UnivariateRealOptimizer;
import org.apache.commons.math.optimization.univariate.BrentOptimizer;
import pal.alignment.Alignment;
import pal.tree.Node;
import pal.tree.Tree;
import tdg.cli.GeneticCodeOption;
import tdg.cli.GlobalsOptions;
import tdg.model.Fitness;
import tdg.model.TDGCodonModel;
import tdg.model.TDGGlobals;
import tdg.trees.RerootedTreeIterator;
import tdg.utils.Functions;
import tdg.utils.GeneticCode;
import tdg.utils.PhyloUtils;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.*;

public class BranchLengthOptimiser {

    private static final double MIN_BRANCH_LENGTH = 1E-7;
    private static final double MAX_BRANCH_LENGTH = 10;

    private final List<SiteCalculator> siteCalculators = Lists.newArrayList();

    private final Tree originalTree;
    private final int threads;

    public static void main(String[] args) throws Exception {
        BranchLengthOptimiser blo = new BranchLengthOptimiser(args);
        blo.run();
    }

    public BranchLengthOptimiser(String... args) throws Exception{
        Options o = new Options();
        JCommander jc = new JCommander(o);
        jc.parse(args);

        threads = o.threads;
        originalTree = PhyloUtils.readTree(o.treeFile);

        Alignment alignment = PhyloUtils.readAlignment(o.alignmentFile);
        TDGGlobals globals = new TDGGlobals(o.globals.tau, o.globals.kappa, o.globals.pi, o.globals.mu, 0);

        List<List<Double>> fitnesses = Files.readLines(new File(o.fitnessFile), Charset.defaultCharset(), new LineProcessor<List<List<Double>>>() {
            List<List<Double>> f = Lists.newArrayList();
            @Override
            public boolean processLine(String s) throws IOException {
                f.add(Lists.transform(Lists.newArrayList(s.split(" ")), Functions.stringToDouble()));
                return true;
            }

            @Override
            public List<List<Double>> getResult() {
                return f;
            }
        });

        createSiteCalculators(fitnesses, alignment, globals);
    }

    private void createSiteCalculators(List<List<Double>> fitnesses, Alignment alignment, TDGGlobals globals) {
        for (int i = 0; i < fitnesses.size(); i++) {

            int site = i + 1;

            Map<String, Integer> p = getGaplessSitePattern(alignment, site);

            List<Integer> a = PhyloUtils.getDistinctAminoAcids(p.values());
            for (int j = 0; j < GeneticCode.AMINO_ACID_STATES; j++) if (!a.contains(j)) a.add(j);

            Fitness f = new Fitness(getOptimiserOrder(Doubles.toArray(fitnesses.get(i)), a), false);

            TDGCodonModel m = new TDGCodonModel(globals, f, a);

            siteCalculators.add(new SiteCalculator(p, f, m));
        }
    }

    private void run() throws Exception{
        System.out.printf("Original tree length: %s\n", PhyloUtils.getTotalTreeLength(originalTree));
        Tree optimisedTree = optimiseBranchLength(originalTree);
        System.out.printf("New tree length: %s\n", PhyloUtils.getTotalTreeLength(optimisedTree));
        System.out.printf("New tree:\n%s\n", PhyloUtils.prettyTreeString(optimisedTree));
    }

    private Tree optimiseBranchLength(final Tree originalTree) throws Exception {

        RerootedTreeIterator rti = new RerootedTreeIterator(originalTree);

        Set<Set<Node>> optimisedBranches = Sets.newHashSet();

        int count = 0;

        // loop through every rerooted tree
        for (Tree tree : rti) {

            System.out.printf("Rerooted tree %s\n", count++);

            // sanity check the current likelihood
            double previousLnL = updateSiteCalculatorTrees(tree);

            System.out.printf("Current lnL: %s\n", previousLnL);

            // get the root for this tree
            Node r = tree.getRoot();

            // loop through each child node of the root
            for (int i = 0; i < r.getChildCount(); i++) {

                final Node c = r.getChild(i);

                Set<Node> branch = Sets.newHashSet(r, c);

                // If we haven't optimised this branch yet
                if (!optimisedBranches.contains(branch)) {
                    optimisedBranches.add(branch);
                    double old = c.getBranchLength();

                    UnivariateRealOptimizer opt = new BrentOptimizer();
                    opt.optimize(new UnivariateRealFunction() {
                        @Override
                        public double value(double branchlength) throws FunctionEvaluationException {
                            return getLikelihoodSum(c, branchlength);
                        }
                    }, GoalType.MAXIMIZE, MIN_BRANCH_LENGTH, MAX_BRANCH_LENGTH);

                    System.out.printf("%s -> %s\t%s -> %s\n", previousLnL, opt.getFunctionValue(), old, opt.getResult());
                    previousLnL = opt.getFunctionValue();
                    c.setBranchLength(opt.getResult());

                }
            }

        }

        return rti.getOriginalRooting();
    }

    private double getLikelihoodSum(final Node node, final double newBranchLength) {
        final ExecutorService threadPool = Executors.newFixedThreadPool(threads);

        final List<Future<Double>> futures = Lists.newArrayList();

        for (final SiteCalculator s : siteCalculators) {
            Future<Double> future = threadPool.submit(new Callable<Double>() {
                @Override
                public Double call() throws Exception {
                    return s.getNodeLikelihood(node, newBranchLength);
                }
            });
            futures.add(future);
        }

        double totalLnL = 0;
        int count = 1;
        for (final Future<Double> f : futures) {
            try {
                totalLnL += f.get();
                System.out.printf("%s / %s                \r", count++, siteCalculators.size());
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }

        threadPool.shutdown();

        return totalLnL;
    }

    /**
     * Loops through each of the sets of fitnesses we've loaded and creates a likelihood calculator for each. This
     * means getting the site pattern for the site, ordering the fitnesses in the way that LikelihoodCalculator (and
     * the optimiser) expect. It is assumed that each list of fitnesses is ordered by site location.
     */
    private double updateSiteCalculatorTrees(final Tree tree) {
        final ExecutorService threadPool = Executors.newFixedThreadPool(threads);

        List<Future<Double>> futures = Lists.newArrayList();

        for (final SiteCalculator sc : siteCalculators) {
            Future<Double> future = threadPool.submit(new Callable<Double>() {
                @Override
                public Double call() throws Exception {
                    return sc.updateTree(tree);
                }
            });
            futures.add(future);
        }

        double sum = 0;
        int count = 1;

        for (Future<Double> f : futures) {
            try {
                sum += f.get();
                System.out.printf("%s / %s                \r", count++, siteCalculators.size());
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }
        System.out.println();

        threadPool.shutdown();
        return sum;
    }

    /**
     * Fitnesses are usually provided loaded in the canonical amino acid order. The likelihood calculator expects
     * fitnesses ordered by the PhyloUtils.getDistinctAminoAcids function (i.e. ordered by observed frequency, then
     * canonical amino acid order). The first fitness should be 0.
     */
    private double[] getOptimiserOrder(double[] fitness, List<Integer> order) {
        List<Double> orderedFitness = Lists.newArrayList();

        for (Integer aminoacid : order) {
            orderedFitness.add(fitness[aminoacid]);
        }

        for (int i = 0; i < GeneticCode.AMINO_ACID_STATES; i++) {
            if (!order.contains(i)) {
                orderedFitness.add(fitness[i]);
            }
        }

        return Doubles.toArray(orderedFitness);
    }

    /**
     * Given an alignment and site positions, returns a site pattern with bad or STOP codons replaced with unknown
     * states.
     */
    private Map<String, Integer> getGaplessSitePattern(Alignment a, int site) {
        Map<String, Integer> sitePattern = PhyloUtils.getCodonsAtSite(a, site);

        // Remove any stop codons and treat them as gaps
        for (Map.Entry<String, Integer> e : sitePattern.entrySet()) {
            if (GeneticCode.getInstance().isUnknownCodonState(GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(e.getValue()))) {
                sitePattern.put(e.getKey(), GeneticCode.UNKNOWN_STATE);
            }
        }

        return sitePattern;
    }

    class Options {
        @Parameter(names = "-tree", description = "Tree file in Newick format", required = true)
        public String treeFile;

        @Parameter(names = "-seq", description = "Codon alignment file in Phylip sequential format.", required = true)
        public String alignmentFile;

        @Parameter(names = "-fitnessfile", required = true)
        public String fitnessFile;

        @Parameter(names = "-threads", description = "The number of threads to use.", required = false)
        public int threads = 1;

        @ParametersDelegate
        public GlobalsOptions globals = new GlobalsOptions();

        @ParametersDelegate
        public GeneticCodeOption gc = new GeneticCodeOption();
    }
}
