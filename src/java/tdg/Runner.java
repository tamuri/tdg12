package tdg;

import pal.tree.Tree;
import tdg.model.Prior;
import tdg.model.TDGGlobals;
import tdg.utils.Pair;

public interface Runner {
    public static final double MIN_BRANCH_LENGTH = 1e-9;
    public static final double MAX_BRANCH_LENGTH = 10;

    public Pair<Double, TDGGlobals> optimiseMutationModel(final Tree tree, final TDGGlobals globals, final FitnessStore fitnessStore, final Prior prior);

    public Pair<Double, Tree> optimiseBranchLengths(final Tree originalTree, final TDGGlobals globals, final FitnessStore fitnessStore, final Prior prior);

    public double optimiseFitness(final Tree tree, final TDGGlobals globals, final FitnessStore fitnessStore, final Prior prior);

    public void close();
}
