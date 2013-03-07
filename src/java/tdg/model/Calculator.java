package tdg.model;

import pal.alignment.Alignment;
import pal.tree.Tree;
import tdg.utils.PhyloUtils;

import java.util.Map;

/**
 * User: atamuri
 * Date: 07/03/2013 10:15
 */
public class Calculator {

    public double getLogLikelihood(Alignment alignment, Tree tree, int site, Fitness fitness, TDGGlobals globals) {
        Map<String, Integer> states = PhyloUtils.getCodonsAtSite(alignment, site);
        TDGCodonModel model = new TDGCodonModel(globals, fitness, PhyloUtils.getDistinctAminoAcids(states.values()));

        LikelihoodCalculator calculator = new LikelihoodCalculator(tree, states, null);
        calculator.setParameters(fitness);
        calculator.addCladeModel("ALL", model);

        calculator.getStorage();
        double result = calculator.function(calculator.getMinimisationParameters().getParameters());
        calculator.releaseStorage();

        return result;
    }


}
