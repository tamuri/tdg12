package tdg.utils;

import com.google.common.base.Function;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class Functions {
    public static Function<String, Double> stringToDouble() {
        return new Function<String, Double>() {
            public Double apply(String s) {
                return Double.parseDouble(s);
            }
        };
    }

    public static Function<String, Integer> stringToInt() {
        return new Function<String, Integer>() {
            public Integer apply(String s) {
                return Integer.parseInt(s);
            }
        };
    }

    public static Function<Integer, String> codonIndexToTLA() {
        return new Function<Integer, String>() {
            @Override
            public String apply(Integer integer) {
                return GeneticCode.getInstance().getCodonTLA(integer);
            }
        };
    }
}
