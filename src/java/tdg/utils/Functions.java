package tdg.utils;

import com.google.common.base.Function;


public class Functions {
   public static Function<String, Double> stringToDouble() {
       return new Function<String,Double>() {
           public Double apply(String s) { return Double.parseDouble(s); }
       };
   }

    public static Function<String, Integer> stringToInt() {
        return new Function<String, Integer>() {
            public Integer apply(String s) { return Integer.parseInt(s); }
        };
    }
}
