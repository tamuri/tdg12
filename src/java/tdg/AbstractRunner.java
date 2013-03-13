package tdg;


import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
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
import tdg.model.TDGGlobals;
import tdg.trees.RerootedTreeIterator;
import tdg.utils.Pair;

import java.sql.Timestamp;
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
        optimiser.setConvergenceChecker(new SimpleScalarValueChecker(-1, 1E-6));

        RealPointValuePair optima;

        // The mutation only fitness is used *only* at the start of the estimation procedure. We can use this to assume
        // that: i) All sites will have the same Fitness & TDGCodonModel and ii) All branch lengths are Constants.INITIAL_BRANCH_LENGTH
        final boolean mutationOnly = (fitnessStore.getFitness(1) == Fitness.getMutationOnlyFitness());

        // Tree and FitnessStore does not need to change, so allow runners to do something intelligent
        runnerSetTree(tree);
        runnerSetFitnessStore(fitnessStore);

        try {
            optima = optimiser.optimize(new MultivariateRealFunction() {
                private int evaluations = 0;
                @Override
                public double value(double[] point) throws FunctionEvaluationException, IllegalArgumentException {

                    final TDGGlobals globals = new TDGGlobals(point);

                    if (globals.getTau() < 0) return Constants.VERY_BAD_LIKELIHOOD;
                    if (globals.getKappa() < 0) return Constants.VERY_BAD_LIKELIHOOD;
                    if (globals.getMu() < 0) return Constants.VERY_BAD_LIKELIHOOD;

                    double d = runnerGetLogLikelihood(tree, fitnessStore, globals);

                    evaluations++;
                    if (evaluations % 25 == 0)
                        System.out.printf("%s - optimiseMutationModel evaluation %s: %s -> %s\n", new Timestamp(System.currentTimeMillis()), evaluations, globals.toString(), d);

                    return d;
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

        // System.out.printf("original tree lnl: %s\n", updateSiteCalculatorTrees(originalTree, globals, fitnessStore));

        double previousLikelihood = Double.NEGATIVE_INFINITY;
        double newLikelihood = Double.NEGATIVE_INFINITY;

        RealConvergenceChecker convergenceChecker = new SimpleScalarValueChecker(0, 1e-6);

        int iteration = 0;

        RerootedTreeIterator rti = new RerootedTreeIterator(originalTree);

        boolean converged = false;

        while (!converged) {

            // Save branches we've already optimised so we don't do them again
            Set<Set<Node>> optimisedBranches = Sets.newHashSet();

            iteration++;

            int count = 0;

            // loop through every rerooted tree
            for (Tree tree : rti) {

                //System.out.printf("%s - optimiseBranchLengths rerooted tree %s\n", new Timestamp(System.currentTimeMillis()), count++);

                // sanity check the current likelihood
                newLikelihood = updateSiteCalculatorTrees(tree, globals, fitnessStore);

                //System.out.printf("%s - optimiseBranchLengths current lnL: %s\n", new Timestamp(System.currentTimeMillis()), previousLnL);

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
                                    double d = getLikelihoodSum(child.getNumber(), branchlength);
                                    return d;
                                }
                            }, GoalType.MAXIMIZE, MIN_BRANCH_LENGTH, MAX_BRANCH_LENGTH, child.getBranchLength());

                            //System.out.printf("%s - %s -> %s\t%s -> %s\n", new Timestamp(System.currentTimeMillis()), previousLnL, opt.getFunctionValue(), old, opt.getResult());

                            newLikelihood = opt.getFunctionValue();
                            child.setBranchLength(opt.getResult()); // TODO: this only updates the local tree
                            updateBranchLength(child.getNumber(), opt.getResult()); // TODO: this updates the saved partials
                        } catch (Exception e) {
                            throw new RuntimeException(e);
                        }

                    }
                }
            }

            converged = convergenceChecker.converged(iteration,
                    new RealPointValuePair(new double[]{}, previousLikelihood),
                    new RealPointValuePair(new double[]{}, newLikelihood));

            previousLikelihood = newLikelihood;

            System.out.printf("%s - Branch estimation %s -> %s\n", new Timestamp(System.currentTimeMillis()), iteration, newLikelihood);

        }

        return Pair.of(newLikelihood, rti.getOriginalRooting());
    }

    protected abstract double getLikelihoodSum(int child, double branchlength);

    protected abstract double updateSiteCalculatorTrees(Tree tree, TDGGlobals globals, FitnessStore fitnessStore);

    protected abstract void updateBranchLength(int child, double branchlength);

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
