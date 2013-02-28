package tdg;


import pal.alignment.Alignment;
import pal.tree.Tree;
import tdg.model.Fitness;
import tdg.model.LikelihoodCalculator;
import tdg.model.TDGCodonModel;
import tdg.model.TDGGlobals;
import tdg.utils.GeneticCode;
import tdg.utils.PhyloUtils;

import java.util.List;
import java.util.Map;

public abstract class AbstractRunner implements Runner {
    public LikelihoodCalculator getLikelihoodCalculatorForSite(Tree tree, Alignment alignment, TDGGlobals g, int site, boolean useapprox) {
        Map<String, Integer> states = PhyloUtils.getCleanedCodons(alignment, site);
        List<Integer> aminoAcids = PhyloUtils.getDistinctAminoAcids(states.values());

        if (!useapprox) {
            for (int i = 0; i < GeneticCode.AMINO_ACID_STATES; i++) {
                if (!aminoAcids.contains(i)) aminoAcids.add(i);
            }
        }

        Fitness f = new Fitness(new double[aminoAcids.size()], false);


        TDGCodonModel model = new TDGCodonModel(g, f, aminoAcids);
        model.updateModel();
        LikelihoodCalculator calculator = new LikelihoodCalculator(tree, states, null);
        calculator.setParameters(f);
        calculator.addCladeModel("ALL", model);
        return calculator;
    }
}
