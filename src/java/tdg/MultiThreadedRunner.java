package tdg;

import com.google.common.collect.Lists;
import com.google.common.util.concurrent.AtomicDouble;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.optimization.GoalType;
import org.apache.commons.math.optimization.RealConvergenceChecker;
import org.apache.commons.math.optimization.RealPointValuePair;
import org.apache.commons.math.optimization.SimpleScalarValueChecker;
import org.apache.commons.math.optimization.direct.DirectSearchOptimizer;
import org.apache.commons.math.optimization.direct.NelderMead;
import org.apache.commons.math.random.MersenneTwister;
import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;
import pal.alignment.Alignment;
import pal.tree.Tree;
import tdg.model.*;
import tdg.optim.LikelihoodFunctionWrapper;
import tdg.utils.CoreUtils;
import tdg.utils.Pair;
import tdg.utils.PhyloUtils;

import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;


/**
 * REMEMBER: lots of cpus != lots of memory
 * imagine that you're running with 100s of cpus but <1GB memory
 * Think about a small tree with lots of sites - can't store all likelihood calculators
 */
public class MultiThreadedRunner extends AbstractRunner {
    private final ExecutorService threadPool;
    private int[] sites;

    public MultiThreadedRunner(Alignment alignment, int threads) {
        setAlignment(alignment);
        this.sites = CoreUtils.seqi(1, alignment.getSiteCount() / 3);
        this.threadPool = Executors.newFixedThreadPool(threads);
    }

    public void setSites(int[] sites) {
        this.sites = sites;
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
    public double runnerGetLogLikelihood(final Tree tree, final FitnessStore fitnessStore, final TDGGlobals globals, final Prior prior) {
        List<Future<Double>> futures = Lists.newArrayList();

        for (final int site : sites) {
            Future<Double> future = threadPool.submit(new Callable<Double>() {
                @Override
                public Double call() throws Exception {
                    Fitness fitness = fitnessStore.getFitness(site);
                    Map<String, Integer> states = PhyloUtils.getCleanedCodons(getAlignment(), site);
                    TDGCodonModel model = new TDGCodonModel(globals, fitness, PhyloUtils.getDistinctAminoAcids(states.values()));

                    LikelihoodCalculator calculator = new LikelihoodCalculator(tree, states, prior);
                    calculator.setParameters(fitness);
                    calculator.addCladeModel("ALL", model);

                    calculator.getStorage();
                    double result = calculator.function(calculator.getMinimisationParameters().getParameters());
                    calculator.releaseStorage();

                    return result;
                }
            });

            futures.add(future);
        }

        double total = 0;
        for (double x : CoreUtils.getFutureResults(futures)) total += x;
        return total;
    }

    @Override
    public double optimiseFitness(final Tree tree, final TDGGlobals globals, final FitnessStore fitnessStore, final Prior prior) {
        final List<Future<Pair<Integer, Fitness>>> futures = Lists.newArrayList();

        final AtomicDouble total = new AtomicDouble(0);
        final RandomData randomData = new RandomDataImpl(new MersenneTwister(987654321));

        for (final int site : sites) {
            Future<Pair<Integer, Fitness>> future = threadPool.submit(new Callable<Pair<Integer, Fitness>>() {
                @Override
                public Pair<Integer, Fitness> call() throws Exception {
                    Map<String, Integer> states = PhyloUtils.getCleanedCodons(getAlignment(), site);
                    List<Integer> residues = PhyloUtils.getDistinctAminoAcids(states.values());

                    LikelihoodCalculator calculator = new LikelihoodCalculator(tree, states, prior);
                    calculator.getStorage();

                    // TODO: have multiple run with multiple initial starting parameters
                    // TODO: handle convergence problems!
                    // Fitness fitness = new Fitness(fitnessStore.getFitness(fitness_i).get(), true);

                    List<Pair<Double, double[]>> optimals = Lists.newArrayList();

                    for (int run = 0; run < 3; run++) {
                        double[] f;

                        if (run == 0) {
                            f = fitnessStore.getFitness(site).get();
                        } else if (run == 1) {
                            f = new double[residues.size()];
                        } else {
                            f = new double[residues.size()];
                            for (int j = 0; j < f.length; j++)
                                f[j] = randomData.nextUniform(-Constants.RANDOM_INITIAL_FITNESS_RANGE, Constants.RANDOM_INITIAL_FITNESS_RANGE);
                        }

                        Fitness fitness = new Fitness(f, true);
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
                            // total.getAndAdd(lnl);
                            optimals.add(Pair.of(lnl, fitness.get().clone()));
                            // return Pair.of(fitness_i, fitness);
                        } else {
                            RealPointValuePair optima;
                            try {
                                optima = optimiser.optimize(wrapper, GoalType.MAXIMIZE, calculator.getMinimisationParameters().getParameters());
                            } catch (FunctionEvaluationException e) {
                                throw new RuntimeException(e);
                            }

                            if (optima == null) return null;

                            calculator.function(optima.getPoint());

                            calculator.releaseStorage();

                            //System.out.printf("Site %s - %s\n", final_i, Doubles.join(", ", fitness.get()));
                            //System.out.printf("Site %s - %s\n", final_i, optima.getValue());

                            // total.getAndAdd(optima.getValue());
                            optimals.add(Pair.of(optima.getValue(), fitness.get().clone()));

                            // return Pair.of(fitness_i, fitness);
                        }

                    }

                    int best = 0;
                    for (int i = 1; i < optimals.size(); i++)
                        if (optimals.get(i).first > optimals.get(best).first) best = i;

                    total.getAndAdd(optimals.get(best).first);
                    System.out.printf("%s - %s\n", site, optimals.get(best).first);
                    return Pair.of(site, new Fitness(optimals.get(best).second, true));

                }
            });
            futures.add(future);
        }

        for (Pair<Integer, Fitness> p : CoreUtils.getFutureResults(futures))
            fitnessStore.setFitness(p.first, p.second);

        return total.get();

    }

    private final List<LikelihoodCalculator> calculators = Lists.newArrayList();

    @Override
    public double updateSiteCalculatorTrees(final Tree tree, final TDGGlobals globals, final FitnessStore fitnessStore, final Prior prior) {
        List<Future<Pair<Double, LikelihoodCalculator>>> futures = Lists.newArrayList();

        for (final int site : sites) {
            Future<Pair<Double, LikelihoodCalculator>> future = threadPool.submit(new Callable<Pair<Double, LikelihoodCalculator>>() {
                @Override
                public Pair<Double, LikelihoodCalculator> call() throws Exception {
                    Map<String, Integer> states = PhyloUtils.getCleanedCodons(getAlignment(), site);
                    List<Integer> aminoAcids = PhyloUtils.getDistinctAminoAcids(states.values());

                    TDGCodonModel model = new TDGCodonModel(globals, fitnessStore.getFitness(site), aminoAcids);
                    model.updateModel();

                    LikelihoodCalculator calculator = new LikelihoodCalculator(tree, states, prior);
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

        for (Pair<Double, LikelihoodCalculator> p : CoreUtils.getFutureResults(futures)) {
            total += p.first;
            calculators.add(p.second);
        }

        return total;
    }

    @Override
    public double getLikelihoodSum(final int node, final double newBranchLength) {

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
        for (double x : CoreUtils.getFutureResults(futures)) total += x;

        return total;
    }

    @Override
    public void updateBranchLength(int child, double branchlength) {

        for (LikelihoodCalculator lc : calculators) {
            lc.setBranch(child, branchlength);
        }
    }

    public void close() {
        threadPool.shutdown();
    }


}
