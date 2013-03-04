package tdg;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.google.common.primitives.Ints;
import com.google.common.util.concurrent.AtomicDouble;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.analysis.MultivariateRealFunction;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.optimization.*;
import org.apache.commons.math.optimization.direct.DirectSearchOptimizer;
import org.apache.commons.math.optimization.direct.NelderMead;
import org.apache.commons.math.optimization.univariate.BrentOptimizer;
import pal.alignment.Alignment;
import pal.tree.Node;
import pal.tree.Tree;
import tdg.model.*;
import tdg.optim.LikelihoodFunctionWrapper;
import tdg.trees.RerootedTreeIterator;
import tdg.utils.*;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.*;


/**
 * REMEMBER: lots of cpus != lots of memory
 *  imagine that you're running with 100s of cpus but <1GB memory
 *  Think about a small tree with lots of sites - can't store all likelihood calculators
 */
public class MultiThreadedRunner extends AbstractRunner {
    private final ExecutorService threadPool;

    public MultiThreadedRunner(int threads) {
        this.threadPool =  Executors.newFixedThreadPool(threads);
    }

    /**
     * Optimises mutation model *only* (all fitnesses are set to 0)
     */
    @Override public Pair<Double, TDGGlobals> optimiseMutationModel(final Tree tree, final Alignment alignment, final TDGGlobals startGlobals, final FitnessStore fitnessStore) {
        DirectSearchOptimizer optimiser = new NelderMead();
        optimiser.setConvergenceChecker(new SimpleScalarValueChecker(-1, 1E-3));

        RealPointValuePair optima;

        // The mutation only fitness is used *only* at the start of the estimation procedure. We can use this to assume
        // that: i) All sites will have the same Fitness & TDGCodonModel and ii) All branch lengths are Constants.INITIAL_BRANCH_LENGTH
        final boolean mutationOnly = (fitnessStore.getFitness(1) == Fitness.getMutationOnlyFitness());

        try {

            optima = optimiser.optimize(new MultivariateRealFunction() {
                @Override
                public double value(double[] point) throws FunctionEvaluationException, IllegalArgumentException {

                    final TDGGlobals globals = new TDGGlobals(point);

                    if (globals.getTau() < 0) return Constants.VERY_BAD_LIKELIHOOD;
                    if (globals.getKappa() < 0) return Constants.VERY_BAD_LIKELIHOOD;
                    if (globals.getMu() < 0) return Constants.VERY_BAD_LIKELIHOOD;

                    List<Future<Double>> futures = Lists.newArrayList();

                    final TDGCodonModel mutationModel;

                    if (mutationOnly) {
                        mutationModel = new CachedTDGCodonModel(globals, fitnessStore.getFitness(1), Ints.asList(CoreUtils.seqi(0, 19)));
                        mutationModel.updateModel();
                        mutationModel.getProbabilityMatrix(new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES], Constants.INITIAL_BRANCH_LENGTH);
                    } else {
                        mutationModel = null;
                    }

                    for (int i = 1; i <= alignment.getSiteCount() / 3; i++) {
                        final int site = i;
                        Future<Double> future = threadPool.submit(new Callable<Double>() {
                            @Override
                            public Double call() throws Exception {
                                double result;

                                Map<String, Integer> states = PhyloUtils.getCleanedCodons(alignment, site);

                                LikelihoodCalculator calculator = new LikelihoodCalculator(tree, states, null);
                                calculator.getStorage();

                                if (mutationOnly) {
                                    calculator.addCladeModel("ALL", mutationModel);
                                    result = calculator.calculateLogLikelihood();
                                } else {
                                    Fitness f = fitnessStore.getFitness(site);
                                    TDGCodonModel model = new TDGCodonModel(globals, f, PhyloUtils.getDistinctAminoAcids(states.values()));
                                    calculator.setParameters(f);
                                    calculator.addCladeModel("ALL", model);
                                    result = calculator.function(calculator.getMinimisationParameters().getParameters());
                                }

                                calculator.releaseStorage();

                                return result;
                            }
                        });

                        futures.add(future);
                    }

                    double total = 0;
                    for (double x : getAllResults(futures)) total += x;

                    // System.out.printf("%s = %s\n", g.toString(), total);

                    return total;
                }
            }, GoalType.MAXIMIZE, startGlobals.getOptimiserParameters());
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }

        if (optima == null) throw new RuntimeException("MultithreadedRunner.optimiseMutationModel was unsuccessful.");

        return Pair.of(optima.getValue(), new TDGGlobals(optima.getPoint()));
    }


    @Override public Pair<Double, Tree> optimiseBranchLengths(final Tree originalTree, final Alignment alignment, final TDGGlobals globals, final FitnessStore fitnessStore) {

        RerootedTreeIterator rti = new RerootedTreeIterator(originalTree);

        // Save branches we've already optimised so we don't do them again
        Set<Set<Node>> optimisedBranches = Sets.newHashSet();

        double previousLnL = Double.NEGATIVE_INFINITY;

        final List<LikelihoodCalculator> allCalculators = Lists.newArrayList();

        // loop through every rerooted tree
        for (Tree tree : rti) {

            // System.out.printf("Rerooted tree %s\n", count++);

            // sanity check the current likelihood
            previousLnL = updateSiteCalculatorTrees(allCalculators, tree, alignment, globals, fitnessStore);
            // System.out.printf("Current lnL: %s\n", previousLnL);

            // get the root for this tree
            Node root = tree.getRoot();

            // loop through each child node of the root
            for (int i = 0; i < root.getChildCount(); i++) {

                final Node child = root.getChild(i);

                Set<Node> branch = Sets.newHashSet(root, child);

                // If we haven't optimised this branch yet
                if (!optimisedBranches.contains(branch)) {
                    optimisedBranches.add(branch);
                    double old = child.getBranchLength();

                    UnivariateRealOptimizer opt = new BrentOptimizer();

                    try {
                        opt.optimize(new UnivariateRealFunction() {
                            @Override
                            public double value(double branchlength) throws FunctionEvaluationException {
                                return getLikelihoodSum(allCalculators, child, branchlength);
                            }
                        }, GoalType.MAXIMIZE, MIN_BRANCH_LENGTH, MAX_BRANCH_LENGTH);

                        // System.out.printf("%s -> %s\t%s -> %s\n", previousLnL, opt.getFunctionValue(), old, opt.getResult());
                        previousLnL = opt.getFunctionValue();
                        child.setBranchLength(opt.getResult());

                    } catch (Exception e) {
                        throw new RuntimeException(e);
                    }

                }
            }
        }

        return Pair.of(previousLnL, rti.getOriginalRooting());
    }

    @Override public double optimiseFitness(final Tree tree, final Alignment alignment, final TDGGlobals globals, final FitnessStore fitnessStore) {
        final List<Future<Pair<Integer, Fitness>>> futures = Lists.newArrayList();

        final AtomicDouble total = new AtomicDouble(0);
        
        for (int i = 1; i <= alignment.getSiteCount() / 3; i++) {
            final int site = i;
            Future<Pair<Integer, Fitness>> future = threadPool.submit(new Callable<Pair<Integer, Fitness>>() {
                @Override
                public Pair<Integer, Fitness> call() throws Exception {
                    Map<String, Integer> states = PhyloUtils.getCleanedCodons(alignment, site);
                    List<Integer> residues = PhyloUtils.getDistinctAminoAcids(states.values());

                    LikelihoodCalculator calculator = new LikelihoodCalculator(tree, states, null);
                    calculator.getStorage();

                    Fitness fitness = new Fitness(new double[residues.size()], true);
                    TDGCodonModel model = new TDGCodonModel(globals, fitness, residues);

                    calculator.addCladeModel("ALL", model);
                    calculator.setParameters(fitness);

                    DirectSearchOptimizer optimiser = new NelderMead();
                    optimiser.setMaxEvaluations(Constants.MAX_EVALUATIONS);
                    RealConvergenceChecker convergenceChecker = new SimpleScalarValueChecker(-1, Constants.CONVERGENCE_TOL);
                    optimiser.setConvergenceChecker(convergenceChecker);

                    LikelihoodFunctionWrapper wrapper = new LikelihoodFunctionWrapper();
                    wrapper.setLc(calculator);

                    if (residues.size() == 1) {
                        double lnl = calculator.function(new double[]{});
                        calculator.releaseStorage();
                        //System.out.printf("Site %s - 0.0\n", final_i);
                        //System.out.printf("Site %s - %s\n", final_i, lnl);
                        total.getAndAdd(lnl);
                        return Pair.of(site, fitness);
                    } else {
                        RealPointValuePair optima;
                        try {
                            optima = optimiser.optimize(wrapper, GoalType.MAXIMIZE, calculator.getMinimisationParameters().getParameters());

                        } catch (Exception e) {
                            throw new RuntimeException(e);
                        }

                        if (optima == null) return null;

                        calculator.function(optima.getPoint());

                        calculator.releaseStorage();

                        //System.out.printf("Site %s - %s\n", final_i, Doubles.join(", ", fitness.get()));
                        //System.out.printf("Site %s - %s\n", final_i, optima.getValue());

                        total.getAndAdd(optima.getValue());
                        
                        return Pair.of(site, fitness);
                    }
                }
            });

            futures.add(future);
        }

        for (Pair<Integer, Fitness> p : getAllResults(futures))
            fitnessStore.setFitness(p.first, p.second);

        return total.get();

    }

    private double updateSiteCalculatorTrees(final List<LikelihoodCalculator> calculators, final Tree tree, final Alignment alignment, final TDGGlobals globals, final FitnessStore fitnessStore) {
        List<Future<Pair<Double, LikelihoodCalculator>>> futures = Lists.newArrayList();

        for (int i = 1; i <= alignment.getSiteCount() / 3; i++) {
            final int site = i;
            Future<Pair<Double, LikelihoodCalculator>> future = threadPool.submit(new Callable<Pair<Double, LikelihoodCalculator>>() {
                @Override
                public Pair<Double, LikelihoodCalculator> call() throws Exception {
                    Map<String, Integer> states = PhyloUtils.getCleanedCodons(alignment, site);
                    List<Integer> aminoAcids = PhyloUtils.getDistinctAminoAcids(states.values());

                    TDGCodonModel model = new TDGCodonModel(globals, fitnessStore.getFitness(site), aminoAcids);
                    model.updateModel();

                    LikelihoodCalculator calculator = new LikelihoodCalculator(tree, states, null);
                    calculator.getStorage();

                    calculator.setParameters(fitnessStore.getFitness(site));
                    calculator.addCladeModel("ALL", model);

                    double d = calculator.calculateLogLikelihood();
                    calculator.releaseStorage();

                    return Pair.of(d, calculator);
                }
            });

            futures.add(future);
        }

        double total = 0;
        calculators.clear();

        for (Pair<Double, LikelihoodCalculator> p : getAllResults(futures)) {
            total += p.first;
            calculators.add(p.second);
        }

        return total;
    }

    private double getLikelihoodSum(final List<LikelihoodCalculator> calculators, final Node node, final double newBranchLength) {

        final List<Future<Double>> futures = Lists.newArrayList();

        for (final LikelihoodCalculator s : calculators) {
            Future<Double> future = threadPool.submit(new Callable<Double>() {
                @Override
                public Double call() throws Exception {
                    return s.getNodeLikelihood(node, newBranchLength);
                }
            });
            futures.add(future);
        }

        double total = 0;
        for (double x : getAllResults(futures)) total += x;
        
        return total;
    }

    public void close() {
        threadPool.shutdown();
    }

    private <T> List<T> getAllResults(List<Future<T>> futures) {
        List<T> results = Lists.newArrayList();

        for (Future<T> f : futures) {
            try {
                results.add(f.get());
            } catch (Exception e) {
                e.printStackTrace();
                throw new RuntimeException(e);
            }
        }

        return results;
    }
}
