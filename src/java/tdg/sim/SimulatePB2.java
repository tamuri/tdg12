package tdg.sim;

import com.google.common.base.Charsets;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;
import tdg.utils.Functions;
import tdg.utils.GeneticCode;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

public class SimulatePB2 {
    public static void main(String[] args) throws Exception {
        //String nonhomogpath = "/Users/atamuri/Documents/work/mitochondria/110329_TdG_PB2_FullMaxEnt/pb2.exact.maxent/nonhomog/";
        String nonhomogpath = "/Users/atamuri/Documents/work/mitochondria/paper/response/pb2.no.pen/results.nonhomog/";

        String avpath = nonhomogpath + "fitness.av.sorted.txt";
        String hupath = nonhomogpath + "fitness.hu.sorted.txt";

        List<String> hufitnesses = Files.readLines(new File(hupath), Charsets.UTF_8);
        List<String> avfitnesses = Files.readLines(new File(avpath), Charsets.UTF_8);

        GeneticCode.setCode(GeneticCode.STANDARD_CODE);

        Map<String, StringBuffer> seqout = Maps.newHashMap();

        for (int i = 0; i < hufitnesses.size(); i++) {
            List<Double> avfit = Lists.transform(Arrays.asList(avfitnesses.get(i).split(" ")), Functions.stringToDouble());
            List<Double> hufit = Lists.transform(Arrays.asList(hufitnesses.get(i).split(" ")), Functions.stringToDouble());

            SimulatorHetero s = new SimulatorHetero();
            s.fitness1 = Doubles.toArray(avfit);
            s.fitness2 = Doubles.toArray(hufit);
            s.residues1 = new char[]{'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};
            s.residues2 = new char[]{'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};
            s.gamma = 0;
            s.gc = GeneticCode.getInstance();
            s.kappa = 7.8640425;
            s.mu = 3.132626667;
            s.pi = new double[]{2.36114375000000e-01, 1.95774375000000e-01, 3.65944375000000e-01, 2.02166875000000e-01};

            s.sites = 1;
            s.tau = 1.25010000000000e-02;
            s.tree = "/Users/atamuri/Documents/work/tdg10/etc/PB2_FMutSel0.parsedTree.out";

            s.run(s.tree, s.tau, s.kappa, s.pi, s.fitness1, s.residues1, s.fitness2, s.residues2, 1, s.mu);

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
