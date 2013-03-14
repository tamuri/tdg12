package tdg.distributed;

import pal.tree.Tree;
import tdg.FitnessStore;
import tdg.model.Prior;
import tdg.model.TDGGlobals;

/**
 * User: atamuri
 * Date: 05/03/2013 15:47
 */
public interface ServiceAPI {
    public void setSites(int[] sites);

    public double optimiseMutationModel(TDGGlobals globals, Prior prior);

    public int getThreads();

    public void setFitnessStore(FitnessStore fitnessStore);

    public void setTree(Tree tree);

    public double updateLikelihoodCalculators(TDGGlobals globals, Prior prior);

    public double getNodeLikelihood(int node, double newBranchLength);

    public void setBranchLength(int node, double branchlength);

    public double optimiseFitness(final TDGGlobals globals, final Prior prior);

    public FitnessStore getFitnessStore();
}
