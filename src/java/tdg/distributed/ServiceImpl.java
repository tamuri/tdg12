package tdg.distributed;

import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;
import pal.alignment.Alignment;
import pal.tree.Node;
import pal.tree.Tree;
import tdg.FitnessStore;
import tdg.MatrixArrayPool;
import tdg.model.*;
import tdg.utils.CoreUtils;
import tdg.utils.PhyloUtils;

import java.util.List;
import java.util.Map;

/**
 * User: atamuri
 * Date: 05/03/2013 15:50
 */
public class ServiceImpl implements ServiceAPI {
    private final Alignment alignment;
    private int threads;
    private int[] sites;

    private Tree tree;
    private FitnessStore fitnessStore;
    private final Calculator calculator;

    public ServiceImpl(Alignment alignment, int threads) {
        this.alignment = alignment;
        this.threads = threads;
        this.calculator = new Calculator();
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
        this.sites = sites;
        System.out.printf("This slave is responsible for sites %s\n", Ints.join(", ", sites));
    }

    @Override
    public void setFitnessStore(FitnessStore fitnessStore) {
        this.fitnessStore = fitnessStore;
    }

    @Override
    public double optimiseMutationModel(TDGGlobals globals) {

        // TODO: make multi-threaded

        double sum = 0;
        for (int pos = 0; pos < sites.length; pos++) {
            sum += this.calculator.getLogLikelihood(alignment, tree, sites[pos], fitnessStore.getFitness(pos + 1), globals);
        }

        return sum;
    }

    private List<LikelihoodCalculator> calculators = Lists.newArrayList();

    @Override
    public double updateLikelihoodCalculators(TDGGlobals globals) {
        // TODO: make multi-threaded
        calculators.clear();

        double sum = 0;

        for (int pos = 0; pos < sites.length; pos++) {

            Map<String, Integer> states = PhyloUtils.getCleanedCodons(this.alignment, sites[pos]);


            TDGCodonModel model = new TDGCodonModel(globals, fitnessStore.getFitness(pos + 1), PhyloUtils.getDistinctAminoAcids(states.values()));
            model.updateModel();

            LikelihoodCalculator calculator = new LikelihoodCalculator(this.tree, states, null);
            calculator.getStorage();
            calculator.setParameters(fitnessStore.getFitness(pos + 1));
            calculator.addCladeModel("ALL", model);

            double result = calculator.calculateLogLikelihood();
            sum += result;

            calculator.releaseStorage();

            calculators.add(calculator);
        }

        return sum;
    }

    @Override
    public double getNodeLikelihood(Node node, double newBranchLength) {
        double sum = 0;
        for (LikelihoodCalculator c : calculators) {
            sum += c.getNodeLikelihood(node, newBranchLength);
        }

        return sum;
    }

    @Override
    public void setBranchLength(Node node, double branchlength) {
        for (LikelihoodCalculator c : calculators) {
            c.setBranch(node, branchlength);
        }
    }
}
