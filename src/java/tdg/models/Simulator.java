package tdg.models;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Doubles;
import pal.misc.Identifier;
import pal.tree.Node;
import pal.tree.Tree;
import tdg.Options;
import tdg.cli.CharArrayConverter;
import tdg.cli.DoubleArrayConverter;
import tdg.cli.DoubleConverter;
import tdg.cli.GeneticCodeConverter;
import tdg.models.parameters.Fitness;
import tdg.utils.GeneticCode;
import tdg.utils.PhyloUtils;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * @author Asif Tamuri
 * @version $Id: Simulator.java 152 2010-11-08 11:10:01Z tamuri $
 */
public class Simulator {

    private TDGCodonModel codonModel;
    private Map<Identifier, int[]> seq;
    private double[] Pt = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];

    public static void main(String... args) {
        Simulator s = new Simulator();
        JCommander jc = new JCommander(s);
        jc.setProgramName("tdg.models.Simulator");

        if (args.length == 0) {
            jc.usage();
            System.out.println("Options preceded by an asterisk are required.");
            System.out.println("Example: -t PB2.tree -sites 10 -tau 1e-6 -kappa 8.0004 -pi 0.21051,0.19380,0.40010,0.19559 -mu 3.0 -fitness 0,0,0 -characters E,K,R");
            System.exit(0);
        } else {
            jc.parse(args);
        }

        s.run(s.tree, s.tau, s.kappa, s.pi, s.fitness, s.residues, s.sites, s.mu);
    }

    private void run(String treeFile, double tau, double kappa, double[] pi, double[] fitness, char[] residues, int sites, double mu) {
        Tree tree = PhyloUtils.readTree(treeFile);

        // Initialise sequence store
        seq = Maps.newHashMap();
        seq.put(tree.getRoot().getIdentifier(), new int[sites]);
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            int[] store = new int[sites];
            seq.put(tree.getExternalNode(i).getIdentifier(), store);
        }
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            int[] store = new int[sites];
            seq.put(tree.getInternalNode(i).getIdentifier(), store);
        }

        TDGGlobals globals = new TDGGlobals(tau, kappa, pi, mu, gamma);

        List<Integer> aa = Lists.newArrayList();
        for (char c : residues) {
            aa.add(GeneticCode.getInstance().getAminoAcidIndexByChar(c));
        }

        // Initialise codon model (Q)
        codonModel = new TDGCodonModel(globals, new Fitness(fitness, false), aa);
        codonModel.updateModel();


        // Simulate sequence at the root
        Node root = tree.getRoot();
        for (int i = 0; i < sites; i++) {
            seq.get(root.getIdentifier())[i] = selectRandomCharacter(codonModel.getCodonFrequencies());
        }

        // Traverse tree
        downTree(root, sites);

        // Print tips
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            Identifier n = tree.getExternalNode(i).getIdentifier();
            System.out.printf("%s     ", n.getName());
            for (int j = 0; j < sites; j++) {
                System.out.printf("%s", GeneticCode.getInstance().getCodonTLA(seq.get(n)[j]));
            }
            System.out.println("");
        }

        System.err.printf("Pi: %s\n\n", Doubles.join(",", codonModel.getAminoAcidFrequencies()));
    }

    private void downTree(Node parent, int sites) {
        for (int i = 0; i < parent.getChildCount(); i++) {
            Node child = parent.getChild(i);

            codonModel.getProbabilityMatrix(Pt, child.getBranchLength());

            for (int j = 0; j < sites; j++) {
                int row = seq.get(parent.getIdentifier())[j];
                seq.get(child.getIdentifier())[j] = selectRandomCharacter(Arrays.copyOfRange(Pt, row * GeneticCode.CODON_STATES, row * GeneticCode.CODON_STATES + GeneticCode.CODON_STATES));
            }

            downTree(child, sites);
        }
    }

    private int selectRandomCharacter(double[] array) {
        double rnd = Math.random();
        //Random random = new Random(123456789);
        //random.nextDouble();
        int character = -1;
        double sum = 0.0;
        for (int i = 0; i < array.length; i++) {
            sum += array[i];
            if (sum > rnd) {
                character = i;
                break;
            }
        }
        return character;
    }

    // COMMAND LINE ARGUMENTS
    @Parameter(names = "-t", description = "Tree file in Newick format", required = true)
    private String tree;

    @Parameter(names = "-sites", description = "Number of sites to simulate.", required = true)
    private int sites;

    @Parameter(names = "-tau", description = "Rate for multiple substitutions.", converter = DoubleConverter.class, required = true)
    private double tau;

    @Parameter(names = "-kappa", description = "Transition/transversion bias.", converter = DoubleConverter.class, required = true)
    private double kappa;

    @Parameter(names = "-pi", description = "Comma-separated base nucleotide frequencies (T,C,A,G).", converter = DoubleArrayConverter.class, required = true)
    private double[] pi;

    @Parameter(names = "-fitness", description = "Comma-separated fitness coefficients.", converter = DoubleArrayConverter.class, required = true)
    private double[] fitness;

    @Parameter(names = "-characters", description = "Comma-separated amino acids (matching fitness coefficents).", converter = CharArrayConverter.class, required = true)
    private char[] residues;

    @Parameter(names = "-mu", description = "Branch/rate scaling factor.", converter = DoubleConverter.class, required = true)
    private double mu;

    @Parameter(names = "-gamma", description = "Error rate,", converter = DoubleConverter.class)
    private double gamma = 0.;

    @Parameter(names = "-gc", description = "Genetic code. 'standard' or 'vertebrate_mit'.", converter = GeneticCodeConverter.class, required = true)
    private GeneticCode gc;
    
}