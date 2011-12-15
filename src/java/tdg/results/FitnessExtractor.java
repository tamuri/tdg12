package tdg.results;

import com.google.common.base.Charsets;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;
import tdg.utils.GeneticCode;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Map;

/**
 * Given the full output of the TdG12 swMutSel0 analysis, writes a file of fitnesses, ordered by site
 *
 * TODO: Handle heterogeneous fitness output (e.g. Fitness_C1, Fitness_C2 etc.) Currently, you have to do those by hand
 *
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class FitnessExtractor {
    public static void main(String[] args) throws Exception {
        FitnessExtractor fe = new FitnessExtractor();
        fe.extract(args[0]);
    }

    private void extract(String resultsFile) throws Exception {

        Map<Integer, double[]> allFitnesses = Maps.newHashMap();

        BufferedReader reader = Files.newReader(new File(resultsFile), Charsets.US_ASCII);
        String line;

        while ((line = reader.readLine()) != null) {
            if (line.contains("Fitness:")) { // NOTE: Homogeneous fitness only!
                String[] parts = line.split("\\s+");

                double[] fitness = new double[GeneticCode.AMINO_ACID_STATES];
                for (int i = 5; i < 25; i++) {
                    fitness[i - 5] = Double.parseDouble(parts[i].replace(",", ""));
                }

                allFitnesses.put(Integer.parseInt(parts[1]), fitness);
            }
        }

        reader.close();
        
        ArrayList<Integer> orderedKeys = Lists.newArrayList(allFitnesses.keySet());
        Collections.sort(orderedKeys);

        BufferedWriter writer = Files.newWriter(new File("F.txt"), Charsets.US_ASCII);
        
        for (int site : orderedKeys) {
            writer.write(Doubles.join(" ", allFitnesses.get(site)));
            writer.write("\n");
        }

        writer.close();
    }
}
