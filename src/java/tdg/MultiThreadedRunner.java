package tdg;

import com.google.common.base.Function;
import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;
import com.google.common.util.concurrent.AtomicDouble;
import org.apache.commons.math.ConvergenceException;
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
import tdg.model.Fitness;
import tdg.model.LikelihoodCalculator;
import tdg.model.TDGCodonModel;
import tdg.model.TDGGlobals;
import tdg.optim.LikelihoodFunctionWrapper;
import tdg.trees.RerootedTreeIterator;
import tdg.utils.GeneticCode;
import tdg.utils.Pair;
import tdg.utils.PhyloUtils;
import tdg.utils.m;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.*;


/**
 * REMEMBER: lots of cpus != lots of memory
 *  imagine that you're running with 100s of cpus but <1GB memory
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

        try {

            optima = optimiser.optimize(new MultivariateRealFunction() {
                @Override
                public double value(double[] point) throws FunctionEvaluationException, IllegalArgumentException {

                    final TDGGlobals g = new TDGGlobals(point);

                    if (g.getTau() < 0) return Double.NEGATIVE_INFINITY;
                    if (g.getKappa() < 0) return Double.NEGATIVE_INFINITY;
                    if (g.getMu() < 0) return Double.NEGATIVE_INFINITY;

                    List<Future<Double>> futures = Lists.newArrayList();

                    for (int i = 1; i <= alignment.getSiteCount() / 3; i++) {
                        final int final_i = i;
                        Future<Double> future = threadPool.submit(new Callable<Double>() {
                            @Override
                            public Double call() throws Exception {
                                Map<String, Integer> states = PhyloUtils.getCleanedCodons(alignment, final_i);
                                List<Integer> residues = PhyloUtils.getDistinctAminoAcids(states.values());
                                Fitness fit = fitnessStore.getFitness(final_i);
                                TDGCodonModel model = new TDGCodonModel(g, fit, residues);
                                LikelihoodCalculator calculator = new LikelihoodCalculator(tree, states, null);
                                calculator.setParameters(fit);
                                calculator.addCladeModel("ALL", model);
                                double d = calculator.function(calculator.getMinimisationParameters().getParameters());
                                return d;
                            }
                        });

                        futures.add(future);
                    }

                    double total = 0;
                    for (Future<Double> f : futures) {
                        try {
                            total += f.get();
                        } catch (Exception e) {
                            throw new RuntimeException(e);
                        }
                    }

                    System.out.printf("%s = %s\n", g.toString(), total);

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

        Set<Set<Node>> optimisedBranches = Sets.newHashSet();

        int count = 0;

        double previousLnL = Double.NEGATIVE_INFINITY;

        final List<LikelihoodCalculator> allCalculators = Lists.newArrayList();

        // loop through every rerooted tree
        for (Tree tree : rti) {

            System.out.printf("Rerooted tree %s\n", count++);

            // sanity check the current likelihood
            previousLnL = updateSiteCalculatorTrees(allCalculators, tree, alignment, globals, fitnessStore);
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

                    try {
                        opt.optimize(new UnivariateRealFunction() {
                            @Override
                            public double value(double branchlength) throws FunctionEvaluationException {
                                return getLikelihoodSum(allCalculators, c, branchlength);
                            }
                        }, GoalType.MAXIMIZE, MIN_BRANCH_LENGTH, MAX_BRANCH_LENGTH);

                        System.out.printf("%s -> %s\t%s -> %s\n", previousLnL, opt.getFunctionValue(), old, opt.getResult());
                        previousLnL = opt.getFunctionValue();
                        c.setBranchLength(opt.getResult());

                    } catch (Exception e) {
                        throw new RuntimeException(e);
                    }

                }
            }
        }

        // TODO: we want to return the sum log-likelihood
        return Pair.of(previousLnL, rti.getOriginalRooting());
    }

    @Override public double optimiseFitness(final Tree tree, final Alignment alignment, final TDGGlobals globals, final FitnessStore fitnessStore) {
        final List<Future<Pair<Integer, Fitness>>> futures = Lists.newArrayList();

        final AtomicDouble total = new AtomicDouble(0);
        
        for (int i = 1; i <= alignment.getSiteCount() / 3; i++) {
            final int final_i = i;
            Future<Pair<Integer, Fitness>> future = threadPool.submit(new Callable<Pair<Integer, Fitness>>() {
                @Override
                public Pair<Integer, Fitness> call() throws Exception {
                    Map<String, Integer> states = PhyloUtils.getCleanedCodons(alignment, final_i);
                    List<Integer> residues = PhyloUtils.getDistinctAminoAcids(states.values());
                    /*for (int i = 0; i < GeneticCode.AMINO_ACID_STATES; i++) {
                        if (!residues.contains(i)) residues.add(i);
                    }*/

                    LikelihoodCalculator calculator = new LikelihoodCalculator(tree, states, null);

                    Fitness fitness = new Fitness(new double[residues.size()], true);
                    TDGCodonModel model = new TDGCodonModel(globals, fitness, residues);

                    List<String> residueStrings = Lists.transform(residues, new Function<Integer, String>() {
                        int pos = 1;

                        public String apply(Integer input) {
                            return String.format("%s:(%s, %s)", pos++, input, GeneticCode.getInstance().getAminoAcidCharByIndex(input));
                        }
                    });
                    System.out.printf("Site %s - Residues: [%s/%s] { %s }\n", final_i, residues.size(), residues.size(), Joiner.on(", ").join(residueStrings));

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
                        System.out.printf("Site %s - 0.0\n", final_i);
                        System.out.printf("Site %s - %s\n", final_i, lnl);
                        total.getAndAdd(lnl);
                        return Pair.of(final_i, fitness);
                    } else {
                        RealPointValuePair optima;
                        try {
                            optima = optimiser.optimize(wrapper, GoalType.MAXIMIZE, calculator.getMinimisationParameters().getParameters());
                        } catch (Exception e) {
                            throw new RuntimeException(e);
                        }

                        if (optima == null) return null;

                        calculator.function(optima.getPoint());


                        System.out.printf("Site %s - %s\n", final_i, Doubles.join(", ", fitness.get()));
                        System.out.printf("Site %s - %s\n", final_i, optima.getValue());

                        total.getAndAdd(optima.getValue());
                        
                        return Pair.of(final_i, fitness);
                    }
                }
            });

            futures.add(future);
        }

        // List<Fitness> fitnesses = Lists.newArrayList();

        for (Future<Pair<Integer, Fitness>> f : futures) {
            try {
                Pair<Integer, Fitness> out = f.get();
                fitnessStore.setFitness(out.first, out.second);
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
        }

        return total.get();

    }

    private double updateSiteCalculatorTrees(final List<LikelihoodCalculator> calculators, final Tree tree, final Alignment alignment, final TDGGlobals globals, final FitnessStore fitnessStore) {
        List<Future<Pair<Double, LikelihoodCalculator>>> futures = Lists.newArrayList();

        for (int i = 1; i <= alignment.getSiteCount() / 3; i++) {
            final int final_i = i;
            Future<Pair<Double, LikelihoodCalculator>> future = threadPool.submit(new Callable<Pair<Double, LikelihoodCalculator>>() {
                @Override
                public Pair<Double, LikelihoodCalculator> call() throws Exception {
                    Map<String, Integer> states = PhyloUtils.getCleanedCodons(alignment, final_i);
                    List<Integer> aminoAcids = PhyloUtils.getDistinctAminoAcids(states.values());

                    TDGCodonModel model = new TDGCodonModel(globals, fitnessStore.getFitness(final_i), aminoAcids);
                    model.updateModel();
                    LikelihoodCalculator calculator = new LikelihoodCalculator(tree, states, null);
                    calculator.setParameters(fitnessStore.getFitness(final_i));
                    calculator.addCladeModel("ALL", model);

                    return Pair.of(calculator.calculateLogLikelihood(), calculator);
                }
            });

            futures.add(future);
        }

        double sum = 0;

        calculators.clear();

        for (Future<Pair<Double, LikelihoodCalculator>> f : futures) {
            try {
                Pair<Double, LikelihoodCalculator> out = f.get();
                sum += out.first;
                calculators.add(out.second);
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }
        System.out.println();

        return sum;
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

        double totalLnL = 0;
        for (final Future<Double> f : futures) {
            try {
                totalLnL += f.get();
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }

        return totalLnL;
    }

    public void close() {
        threadPool.shutdown();
    }
}
