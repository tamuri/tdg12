package tdg.model;

import tdg.utils.GeneticCode;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * User: atamuri
 * Date: 01/03/2013 12:50
 */
public class CachedTDGCodonModel extends TDGCodonModel {
    private static final int STORE_SIZE = 10;

    private Map<Double, double[]> store = new LinkedHashMap<Double, double[]>() {
        @Override
        protected boolean removeEldestEntry(Map.Entry<Double, double[]> eldest) {
            return this.size() > STORE_SIZE;
        }
    };

    public CachedTDGCodonModel(TDGGlobals globals, Fitness fitness, List<Integer> aminoAcids) {
        super(globals, fitness, aminoAcids);
    }

    @Override
    public void getProbabilityMatrix(double[] matrix, double branchLength) {
        if (store.containsKey(branchLength)) {
            double[] ref = store.get(branchLength);
            System.arraycopy(ref, 0, matrix, 0, matrix.length);
        } else {
            synchronized (this) {
                super.getProbabilityMatrix(matrix, branchLength);
                store.put(branchLength, Arrays.copyOf(matrix, GeneticCode.CODON_STATES * GeneticCode.CODON_STATES));
            }
        }
    }

    @Override
    public void updateModel() {
        super.updateModel();
        store.clear();
    }
}
