package tdg.models;

import com.google.common.base.Charsets;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;
import tdg.utils.Functions;
import tdg.utils.GeneticCode;

import java.io.File;
import java.nio.charset.Charset;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

public class SimulatorTest {
    public static void main(String[] args) throws Exception {
        GeneticCode.initialise(GeneticCode.VERTEBRATE_MITOCHONDRIAL_CODE);

        Map<String, StringBuffer> seqout = Maps.newHashMap();

        double[] d = new double[20];
        Arrays.fill(d, 0.0);
        List<Double> fit = Doubles.asList(d);
        char[] res = new char[]{'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};

        List<String> fits = Files.readLines(new File("./fitness.for.synth.txt"), Charsets.UTF_8);

        for (int i = 0; i < 1000; i++) {

           d = Doubles.toArray(Lists.transform(Arrays.asList(fits.get(i).split(" ")), Functions.stringToDouble()));

            Simulator s = new Simulator();

            s.run("./tree.256.tree", // tree
                    1.0e-02, // tau
                    2.0, // kappa
                    new double[]{0.25, 0.25, 0.25, 0.25}, // s.pi
                    d, // fitnesses
                    res, // residues
                    1, // sites
                    1.0); // mu

            for (Map.Entry<String, String> x : s.seqout.entrySet()) {
                if (!seqout.containsKey(x.getKey())) {
                    seqout.put(x.getKey(), new StringBuffer());
                }
                seqout.get(x.getKey()).append(x.getValue());
            }
        }

        System.out.println();

        for (Map.Entry<String, StringBuffer> x : seqout.entrySet()) {
            System.out.printf("%s      %s\n", x.getKey(), x.getValue().toString());
        }

    }
}
