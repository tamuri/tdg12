package tdg.analysis;

import com.google.common.base.Charsets;
import com.google.common.collect.Maps;
import com.google.common.io.Files;
import com.google.gson.Gson;
import com.google.gson.reflect.TypeToken;

import java.io.File;
import java.io.FileNotFoundException;
import java.lang.reflect.Type;
import java.util.*;

/**
 * @author Asif Tamuri
 * @version $Id$
 */
public class Defaults {
    static String FITNESS_DEFAULT_FILENAME = "fitness.defaults";
    static Map<Integer, Map<Integer, Double>> FITNESS_DEFAULTS = Maps.newHashMap();
    static {
        Gson g = new Gson();
        Type collectionType = new TypeToken<Map<Integer, Map<Integer, Double>>>(){}.getType();
        try {
            FITNESS_DEFAULTS = g.fromJson(Files.newReader(new File(FITNESS_DEFAULT_FILENAME), Charsets.UTF_8), collectionType);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    public static double getDefaultFitness(int site, int residue) {
        if (FITNESS_DEFAULTS.containsKey(site)) {
            return FITNESS_DEFAULTS.get(site).get(residue);
        } else {
            return 0.;
        }
    }

    public static Map<Integer, Double> getDefaultFitnesses(int site) {
        if (FITNESS_DEFAULTS.containsKey(site)) {
            return FITNESS_DEFAULTS.get(site);
        } else {
            return null;
        }
    }

    public static void main(String[] args) throws Exception {
        /*for (Map.Entry<Integer, Map<Integer, Double>> e : FITNESS_DEFAULTS.entrySet()) {
            System.out.printf("Site = %s\n", e.getKey());
            for (Map.Entry<Integer, Double> f : e.getValue().entrySet()) {
                System.out.printf("\tResidue = %s, Fitness = %s\n", f.getKey(), f.getValue());
            }
        }*/

        // compare two fitness JSON files
        /*
        Gson g = new Gson();
        Type collectionType = new TypeToken<Map<Integer, Map<Integer, Double>>>(){}.getType();

        Map<Integer, Map<Integer, Double>> approximateF;
        Map<Integer, Map<Integer, Double>> fullF;

        approximateF = g.fromJson(Files.newReader(new File("nm.fitness.approx"), Charsets.UTF_8), collectionType);
        fullF = g.fromJson(Files.newReader(new File("nm.fitness.full"), Charsets.UTF_8), collectionType);


        for (Map.Entry<Integer, Map<Integer, Double>> e : approximateF.entrySet()) {
            if (fullF.containsKey(e.getKey())) {
                Map<Integer, Double> full = fullF.get(e.getKey());
                for (Map.Entry<Integer, Double> f : e.getValue().entrySet()) {
                    double diff = Math.abs(f.getValue() - full.get(f.getKey()));
                    System.out.printf("%s %s %s %s %s\n", e.getKey(), f.getKey(), f.getValue(), full.get(f.getKey()), diff);
                }
            }
        }

        */

        List<String> l = Files.readLines(new File("/Users/atamuri/Documents/work/flu/evomodel/tdg10/conserved.sites"), Charsets.UTF_8);
        for (String line : l) {
            String[] parts = line.split(" ");
            System.out.printf("%s ", parts[0]);
            int aa = Integer.parseInt(parts[1]);

            for (int i = 0; i < 20; i++) System.out.printf("%s ", aa == i ? 1.0 : 0.0);

            System.out.println();

        }

    }
}
