package tdg;

import com.google.common.collect.Lists;
import tdg.utils.GeneticCode;

import java.util.List;

/**
 * User: atamuri
 * Date: 01/03/2013 22:54
 */
// TODO: Move these pools to the classes where they are used: LikelihoodCalculator and TDGCodonModel!
public class MatrixArrayPool {

    public final static List<double[]> POOL = Lists.newLinkedList();
    public final static List<double[][]> POOL2 = Lists.newLinkedList();

    public static int treeSize;

    public static synchronized double[] pop() {
        if (POOL.isEmpty()) {
            return new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];
        } else {
            return POOL.remove(0);
        }
    }

    public static synchronized void push(double[] matrix) {
        POOL.add(matrix);
    }

    public static synchronized double[][] pop2() {
        if (POOL2.isEmpty()) {
            return new double[treeSize][GeneticCode.CODON_STATES];
        } else {
            return POOL2.remove(0);
        }
    }

    public static synchronized void push2(double[][] matrix) {

        POOL2.add(matrix);
    }
}

