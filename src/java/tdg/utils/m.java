package tdg.utils;

import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

/**
 * Little static math utilities
 */
public class m {
    public static double[] repd(double x, int n) {
        double[] o = new double[n];
        for (int i = 0; i < n; i++) {
            o[i] = x;
        }
        return o;
    }

    public static int[] irep(int x, int n) {
        int[] o = new int[n];
        for (int i = 0; i < n; i++) {
            o[i] = x;
        }

        return o;
    }

    public static double[] seqd(double start, double end, double step) {
        double[] o = new double[(int) ((end - start) / step) + 1];
        for (int i = 0; i < o.length; i++) {
            o[i] = (step * (i + 1));
        }
        return o;
    }

    public static int[] seqi(int start, int end) {
        int[] o = new int[end - start + 1];
        for (int i = 0; i < o.length; i++) {
            o[i] = start + i;
        }
        return o;
    }

    public static void main(String[] args) {
        System.out.printf("%s\n", Ints.join(",", seqi(0, 19)));
    }
}
