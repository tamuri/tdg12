package tdg;

import pal.tree.Node;
import pal.tree.Tree;
import tdg.model.Fitness;
import tdg.model.LikelihoodCalculator;
import tdg.model.TDGCodonModel;

import java.util.Map;

/**
 * Each SiteCalculator holds the LikelihoodCalculator and prepares a Callable that can be used to calculate.
 */
public class SiteCalculator {
    final private Map<String, Integer> sitePattern;
    final private Fitness fitness;
    final private TDGCodonModel codonModel;

    private LikelihoodCalculator likelihoodCalculator;

    SiteCalculator(Map<String, Integer> sitePattern, Fitness fitness, TDGCodonModel codonModel) {
        this.sitePattern = sitePattern;
        this.fitness = fitness;
        this.codonModel = codonModel;
    }

    public double updateTree(Tree tree) {
        // TODO: can we update lc parameters rather than creating a new object each time?
        this.likelihoodCalculator = new LikelihoodCalculator(tree, this.sitePattern, null);
        likelihoodCalculator.setParameters(this.fitness);
        likelihoodCalculator.addCladeModel("ALL", this.codonModel);
        return likelihoodCalculator.function(likelihoodCalculator.getMinimisationParameters().getParameters());
    }

    public double getNodeLikelihood(Node node, double branchLength) {
        return likelihoodCalculator.getNodeLikelihood(node, branchLength);
    }
}
