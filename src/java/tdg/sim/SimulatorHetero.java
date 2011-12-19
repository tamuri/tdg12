package tdg.sim;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import pal.misc.Identifier;
import pal.tree.Node;
import pal.tree.Tree;
import tdg.cli.CharArrayConverter;
import tdg.cli.DoubleArrayConverter;
import tdg.cli.DoubleConverter;
import tdg.cli.GeneticCodeConverter;
import tdg.model.Fitness;
import tdg.model.TDGCodonModel;
import tdg.model.TDGGlobals;
import tdg.utils.GeneticCode;
import tdg.utils.PhyloUtils;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class SimulatorHetero {

    private TDGCodonModel codonModel1, codonModel2;
    private Map<Identifier, int[]> seq;
    private double[] Pt = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];

    public static void main(String... args) {
        SimulatorHetero s = new SimulatorHetero();
        JCommander jc = new JCommander(s);
        jc.setProgramName("tdg.model.Simulator");

        if (args.length == 0) {
            jc.usage();
            System.out.println("AnalyseOptions preceded by an asterisk are required.");
            System.out.println("Example: -t PB2.tree -sites 10 -tau 1e-6 -kappa 8.0004 -pi 0.21051,0.19380,0.40010,0.19559 -mu 3.0 -fitness 0,0,0 -characters E,K,R");
            System.exit(0);
        } else {
            jc.parse(args);
        }

        s.run(s.tree, s.tau, s.kappa, s.pi, s.fitness1, s.residues1, s.fitness2, s.residues2, s.sites, s.mu);
    }

    public void run(String treeFile, double tau, double kappa, double[] pi, double[] fitness1, char[] residues1, double[] fitness2, char[] residues2, int sites, double mu) {
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

        List<Integer> aa1 = Lists.newArrayList();
        for (char c : residues1) {
            aa1.add(GeneticCode.getInstance().getAminoAcidIndexByChar(c));
        }

        // Initialise codon model (Q)
        codonModel1 = new TDGCodonModel(globals, new Fitness(fitness1, false), aa1);
        codonModel1.updateModel();

        List<Integer> aa2 = Lists.newArrayList();
        for (char c : residues2) {
            aa2.add(GeneticCode.getInstance().getAminoAcidIndexByChar(c));
        }

        // Initialise codon model (Q)
        codonModel2 = new TDGCodonModel(globals, new Fitness(fitness2, false), aa2);
        codonModel2.updateModel();

        // Simulate sequence at the root
        Node root = tree.getRoot();
        for (int i = 0; i < sites; i++) {
            seq.get(root.getIdentifier())[i] = selectRandomCharacter(codonModel1.getCodonFrequencies());
        }

        // Traverse tree
        downTree(root, sites);

        seqout = Maps.newHashMap();
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            Identifier n = tree.getExternalNode(i).getIdentifier();
            seqout.put(n.getName(), GeneticCode.getInstance().getCodonTLA(seq.get(n)[0]));
        }

        // Print tips
        /*for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            Identifier n = tree.getExternalNode(i).getIdentifier();
            System.out.printf("%s     ", n.getName());
            for (int j = 0; j < sites; j++) {
                System.out.printf("%s", GeneticCode.getInstance().getCodonTLA(seq.get(n)[j]));
            }
            System.out.println("");
        }*/

        // System.err.printf("Pi: %s\n\n", Doubles.join(",", codonModel.getAminoAcidFrequencies()));
    }

    private void downTree(Node parent, int sites) {
        for (int i = 0; i < parent.getChildCount(); i++) {
            Node child = parent.getChild(i);

            // decide on which codon model to use
            if ((parent.getIdentifier().getName().length() > 0 && parent.getIdentifier().getName().substring(0, 2).equals("Hu"))
                    || child.getIdentifier().getName().substring(0, 2).equals("Hu")) {
                codonModel2.getProbabilityMatrix(Pt, child.getBranchLength());
            } else {
                codonModel1.getProbabilityMatrix(Pt, child.getBranchLength());
            }


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
    public String tree;

    @Parameter(names = "-sites", description = "Number of sites to simulate.", required = true)
    public int sites;

    @Parameter(names = "-tau", description = "Rate for multiple substitutions.", converter = DoubleConverter.class, required = true)
    public double tau;

    @Parameter(names = "-kappa", description = "Transition/transversion bias.", converter = DoubleConverter.class, required = true)
    public double kappa;

    @Parameter(names = "-pi", description = "Comma-separated base nucleotide frequencies (T,C,A,G).", converter = DoubleArrayConverter.class, required = true)
    public double[] pi;

    @Parameter(names = "-fitness1", description = "Comma-separated fitness coefficients.", converter = DoubleArrayConverter.class, required = true)
    public double[] fitness1;

    @Parameter(names = "-fitness2", description = "Comma-separated fitness coefficients.", converter = DoubleArrayConverter.class, required = true)
    public double[] fitness2;

    @Parameter(names = "-characters1", description = "Comma-separated amino acids (matching fitness coefficents).", converter = CharArrayConverter.class, required = true)
    public char[] residues1;

    @Parameter(names = "-characters2", description = "Comma-separated amino acids (matching fitness coefficents).", converter = CharArrayConverter.class, required = true)
    public char[] residues2;


    @Parameter(names = "-mu", description = "Branch/rate scaling factor.", converter = DoubleConverter.class, required = true)
    public double mu;

    @Parameter(names = "-gamma", description = "Error rate,", converter = DoubleConverter.class)
    public double gamma = 0.;

    @Parameter(names = "-gc", description = "Genetic code. 'standard' or 'vertebrate_mit'.", converter = GeneticCodeConverter.class, required = true)
    public GeneticCode gc;

    public Map<String, String> seqout;

}