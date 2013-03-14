package tdg.distributed;

import com.google.common.primitives.Ints;
import pal.alignment.Alignment;
import pal.tree.Tree;
import tdg.FitnessStore;
import tdg.MatrixArrayPool;
import tdg.MultiThreadedRunner;
import tdg.model.Prior;
import tdg.model.TDGGlobals;

/**
 * User: atamuri
 * Date: 05/03/2013 15:50
 */
public class ServiceImpl implements ServiceAPI {
    private int threads;

    private Tree tree;
    private FitnessStore fitnessStore;

    private final MultiThreadedRunner runner;

    public ServiceImpl(Alignment alignment, int threads) {
        this.runner = new MultiThreadedRunner(alignment, threads);
    }

    public void setTree(Tree tree) {
        this.tree = tree;
        MatrixArrayPool.treeSize = tree.getInternalNodeCount();
    }

    public int getThreads() {
        return threads;
    }

    @Override
    public void setSites(int[] sites) {
        runner.setSites(sites);
        System.out.printf("This slave is responsible for sites %s\n", Ints.join(", ", sites));
    }

    @Override
    public void setFitnessStore(FitnessStore fitnessStore) {
        this.fitnessStore = fitnessStore;
    }

    @Override
    public double optimiseMutationModel(TDGGlobals globals, Prior prior) {
        //System.out.printf("%s - optimiseMutationModel in\n", new Timestamp(System.currentTimeMillis()));
        double d =  runner.runnerGetLogLikelihood(this.tree, this.fitnessStore, globals, prior);
        //System.out.printf("%s - optimiseMutationModel out\n", new Timestamp(System.currentTimeMillis()));
        return d;
    }

    @Override
    public double updateLikelihoodCalculators(TDGGlobals globals, Prior prior) {
        //System.out.printf("%s - updatelikelihoodcalculators in\n", new Timestamp(System.currentTimeMillis()));

        double d =  runner.updateSiteCalculatorTrees(this.tree, globals, this.fitnessStore, prior);
        //System.out.printf("%s - updatelikelihoodcalculators out\n", new Timestamp(System.currentTimeMillis()));

        return d;
    }

    @Override
    public double getNodeLikelihood(int node, double newBranchLength) {
        //System.out.printf("%s - getNodeLikelihood in\n", new Timestamp(System.currentTimeMillis()));

        double d =  runner.getLikelihoodSum(node, newBranchLength);
        //System.out.printf("%s - getNodeLikelihood out\n", new Timestamp(System.currentTimeMillis()));

        return d;
    }

    @Override
    public void setBranchLength(int node, double branchlength) {
        runner.updateBranchLength(node, branchlength);
    }

    @Override
    public double optimiseFitness(final TDGGlobals globals, final Prior prior) {
        return runner.optimiseFitness(this.tree, globals, fitnessStore, prior);
    }

    @Override
    public FitnessStore getFitnessStore() {
        return fitnessStore;
    }
}
