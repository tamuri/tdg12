package tdg;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.google.common.primitives.Doubles;
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
    private final Calculator calculator;
    private int[] sites;

    public MultiThreadedRunner(Alignment alignment, int threads) {
        setAlignment(alignment);
        this.sites = CoreUtils.seqi(1, alignment.getSiteCount() / 3);
        this.threadPool =  Executors.newFixedThreadPool(threads);
        this.calculator = new Calculator();
    }

    @Override
    protected void runnerSetTree(Tree tree) {
        // don't need to do anything
    }

    @Override
    protected void runnerSetFitnessStore(FitnessStore fitnessStore) {
        // don't need to do anything
    }

    @Override
    protected double runnerGetLogLikelihood(final Tree tree, final FitnessStore fitnessStore, final TDGGlobals globals) {
        /* final TDGCodonModel mutationModel;

                    if (mutationOnly) {
                        mutationModel = new CachedTDGCodonModel(globals, fitnessStore.getFitness(1), Ints.asList(CoreUtils.seqi(0, 19)));
                        mutationModel.updateModel();
                        mutationModel.getProbabilityMatrix(new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES], Constants.INITIAL_BRANCH_LENGTH);
                    } else {
                        mutationModel = null;
                    }*/

        List<Future<Double>> futures = Lists.newArrayList();

        for (int i : sites) {
            final int site = i;
            Future<Double> future = threadPool.submit(new Callable<Double>() {
                @Override
                public Double call() throws Exception {
                    return calculator.getLogLikelihood(getAlignment(), tree, site, fitnessStore.getFitness(site), globals);

                                /*  Could simplify:          if (mutationOnly) {
                                    calculator.addCladeModel("ALL", mutationModel);
                                    result = calculator.calculateLogLikelihood();
                                } */
                }
            });

            futures.add(future);
        }

        double total = 0;
        for (double x : getAllResults(futures)) total += x;
        return total;
    }

    @Override public double optimiseFitness(final Tree tree, final TDGGlobals globals, final FitnessStore fitnessStore) {
        final List<Future<Pair<Integer, Fitness>>> futures = Lists.newArrayList();

        final AtomicDouble total = new AtomicDouble(0);
        
        for (int i : sites) {
            final int site = i;
            Future<Pair<Integer, Fitness>> future = threadPool.submit(new Callable<Pair<Integer, Fitness>>() {
                @Override
                public Pair<Integer, Fitness> call() throws Exception {
                    Map<String, Integer> states = PhyloUtils.getCleanedCodons(getAlignment(), site);
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

    private final List<LikelihoodCalculator> calculators = Lists.newArrayList();

    @Override public double updateSiteCalculatorTrees(final Tree tree, final TDGGlobals globals, final FitnessStore fitnessStore) {
        List<Future<Pair<Double, LikelihoodCalculator>>> futures = Lists.newArrayList();

        for (int i : sites) {
            final int site = i;
            Future<Pair<Double, LikelihoodCalculator>> future = threadPool.submit(new Callable<Pair<Double, LikelihoodCalculator>>() {
                @Override
                public Pair<Double, LikelihoodCalculator> call() throws Exception {
                    Map<String, Integer> states = PhyloUtils.getCleanedCodons(getAlignment(), site);
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

    @Override public double getLikelihoodSum(final Node node, final double newBranchLength) {

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

    @Override
    protected void updateBranchLength(Node child, double branchlength) {

        for (LikelihoodCalculator lc : calculators) {
            lc.setBranch(child, branchlength);
        }
    }

    public void close() {
        threadPool.shutdown();
    }


}
