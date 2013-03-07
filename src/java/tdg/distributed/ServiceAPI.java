package tdg.distributed;

import pal.tree.Node;
import pal.tree.Tree;
import tdg.FitnessStore;
import tdg.model.TDGGlobals;

/**
 * User: atamuri
 * Date: 05/03/2013 15:47
 */
public interface ServiceAPI {
    public void setSites(int[] sites);

    public double optimiseMutationModel(TDGGlobals globals);

    public int getThreads();

    public void setFitnessStore(FitnessStore fitnessStore);

    public void setTree(Tree tree);

    public double updateLikelihoodCalculators(TDGGlobals globals);

    public double getNodeLikelihood(Node node, double newBranchLength);

    public void setBranchLength(Node node, double branchlength);
}
