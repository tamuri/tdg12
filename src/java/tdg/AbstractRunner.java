package tdg;


import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.analysis.MultivariateRealFunction;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.optimization.GoalType;
import org.apache.commons.math.optimization.RealPointValuePair;
import org.apache.commons.math.optimization.SimpleScalarValueChecker;
import org.apache.commons.math.optimization.UnivariateRealOptimizer;
import org.apache.commons.math.optimization.direct.DirectSearchOptimizer;
import org.apache.commons.math.optimization.direct.NelderMead;
import org.apache.commons.math.optimization.univariate.BrentOptimizer;
import pal.alignment.Alignment;
import pal.tree.Node;
import pal.tree.Tree;
import tdg.model.Fitness;
import tdg.model.LikelihoodCalculator;
import tdg.model.TDGGlobals;
import tdg.trees.RerootedTreeIterator;
import tdg.utils.Pair;

import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.Future;

public abstract class AbstractRunner implements Runner {
    private Alignment alignment;

    /**
     * Optimises TDGGlobal (mutational) parameters
     */
    @Override public Pair<Double, TDGGlobals> optimiseMutationModel(final Tree tree, final TDGGlobals startGlobals, final FitnessStore fitnessStore) {
        DirectSearchOptimizer optimiser = new NelderMead();
        optimiser.setConvergenceChecker(new SimpleScalarValueChecker(-1, 1E-3));

        RealPointValuePair optima;

        // The mutation only fitness is used *only* at the start of the estimation procedure. We can use this to assume
        // that: i) All sites will have the same Fitness & TDGCodonModel and ii) All branch lengths are Constants.INITIAL_BRANCH_LENGTH
        final boolean mutationOnly = (fitnessStore.getFitness(1) == Fitness.getMutationOnlyFitness());

        // Tree and FitnessStore does not need to change, so allow runners to do something intelligent
        runnerSetTree(tree);
        runnerSetFitnessStore(fitnessStore);

        try {
            optima = optimiser.optimize(new MultivariateRealFunction() {
                @Override
                public double value(double[] point) throws FunctionEvaluationException, IllegalArgumentException {

                    final TDGGlobals globals = new TDGGlobals(point);

                    if (globals.getTau() < 0) return Constants.VERY_BAD_LIKELIHOOD;
                    if (globals.getKappa() < 0) return Constants.VERY_BAD_LIKELIHOOD;
                    if (globals.getMu() < 0) return Constants.VERY_BAD_LIKELIHOOD;

                    return runnerGetLogLikelihood(tree, fitnessStore, globals);
                }
            }, GoalType.MAXIMIZE, startGlobals.getOptimiserParameters());
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }

        if (optima == null) throw new RuntimeException("AbstractRunner.optimiseMutationModel was unsuccessful.");

        return Pair.of(optima.getValue(), new TDGGlobals(optima.getPoint()));
    }

    @Override
    public Pair<Double, Tree> optimiseBranchLengths(Tree originalTree, TDGGlobals globals, FitnessStore fitnessStore) {

        // TODO: fitness store does not change (so is tdgglobals but they are pretty lightweight anyway)
        runnerSetFitnessStore(fitnessStore);

        RerootedTreeIterator rti = new RerootedTreeIterator(originalTree);

        // Save branches we've already optimised so we don't do them again
        Set<Set<Node>> optimisedBranches = Sets.newHashSet();

        double previousLnL = Double.NEGATIVE_INFINITY;

        // loop through every rerooted tree
        for (Tree tree : rti) {

            // System.out.printf("Rerooted tree %s\n", count++);

            // sanity check the current likelihood
            previousLnL = updateSiteCalculatorTrees(tree, globals, fitnessStore);
            // System.out.printf("Current lnL: %s\n", previousLnL);

            // get the root for this tree
            Node root = tree.getRoot();

            // System.out.printf("%s\n", tree);


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
                                return getLikelihoodSum(child, branchlength);
                            }
                        }, GoalType.MAXIMIZE, MIN_BRANCH_LENGTH, MAX_BRANCH_LENGTH);

                        //System.out.printf("%s -> %s\t%s -> %s\n", previousLnL, opt.getFunctionValue(), old, opt.getResult());
                        previousLnL = opt.getFunctionValue();
                        child.setBranchLength(opt.getResult()); // TODO: this only updates the local tree
                        updateBranchLength(child, opt.getResult()); // TODO: this updates the saved partials
                    } catch (Exception e) {
                        throw new RuntimeException(e);
                    }

                }
            }
        }

        return Pair.of(previousLnL, rti.getOriginalRooting());
    }

    protected abstract double getLikelihoodSum(Node child, double branchlength);

    protected abstract double updateSiteCalculatorTrees(Tree tree, TDGGlobals globals, FitnessStore fitnessStore);

    protected abstract void updateBranchLength(Node child, double branchlength);

    protected abstract double runnerGetLogLikelihood(final Tree tree, final FitnessStore fitnessStore, final TDGGlobals globals);

    protected abstract void runnerSetTree(final Tree tree);

    protected abstract void runnerSetFitnessStore(final FitnessStore fitnessStore);

    protected <T> List<T> getAllResults(List<Future<T>> futures) {
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

    public Alignment getAlignment() {
        return alignment;
    }

    public void setAlignment(Alignment alignment) {
        this.alignment = alignment;
    }
}
