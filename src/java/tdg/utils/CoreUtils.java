package tdg.utils;

import java.io.IOException;
import java.net.ServerSocket;
import java.text.DecimalFormat;

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

    public static double[] alr(double[] p) {
        // additive log-ratio transformation
        // y = log(p1 / pn), log(p2 / pn), log(p3 / pn) ... log(pn-1 / pn)

        double sum = 0;
        for (double d : p) {
            sum += d;
        }

        double[] y = new double[p.length];

        for (int i = 0; i < p.length; i++) {
            y[i] = Math.log(p[i] / (1 - sum));
        }

        return y;
    }

    public static double[] alr_inv(double[] y) {
        // inverse additive log-ratio transformation
        // p_i = exp(y_i) / (1 + sum(exp(y_i)))
        double[] p = new double[y.length];

        double sum = 0;
        for (double d : y) {
            sum += Math.exp(d);
        }

        for (int i = 0; i < y.length; i++) {
            p[i] = Math.exp(y[i]) / (1 + sum);
        }

        return p;
    }

    public static int findFreePort() {
        ServerSocket socket = null;
        try {
            socket = new ServerSocket(0);
            socket.setReuseAddress(true);
            int port = socket.getLocalPort();
            try {
                socket.close();
            } catch (IOException e) {
// Ignore IOException on close()
            }
            return port;
        } catch (IOException e) {
        } finally {
            if (socket != null) {
                try {
                    socket.close();
                } catch (IOException e) {
                }
            }
        }
        throw new IllegalStateException("Could not find a free TCP/IP port.");
    }
}