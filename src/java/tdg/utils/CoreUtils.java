package tdg.utils;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class CoreUtils {

    public static void print2D(double[][] array, boolean flat) {
        if (flat) {
            for (int i = 0; i < array.length; i++) {
                for (int j = 0; j < array[i].length; j++) {
                    System.out.printf("%s,%s = %s\n", i, j, array[i][j]);
                }
            }
        } else {
            for (int i = 0; i < array.length; i++) {
                for (int j = 0; j < array[i].length; j++) {
                    System.out.printf("%s\t", array[i][j]);
                }
                System.out.println();
            }
        }
    }
    
    public static double sum(double[] list) {
        double t = 0.0;
        for (double d : list) {
            t += d;
        }
        return t;
    }

    public static int[] range(int from, int length) {
        int[] r = new int[length];

        for (int i = 0; i < r.length; i++) {
            r[i] = i + from;
        }

        return r;
    }

}