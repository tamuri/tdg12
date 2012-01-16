package tdg.sim;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import com.beust.jcommander.ParametersDelegate;
import com.google.common.base.Charsets;
import com.google.common.base.Function;
import com.google.common.base.Joiner;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Files;
import tdg.Constants;
import tdg.cli.CharArrayConverter;
import tdg.cli.GeneticCodeOption;
import tdg.cli.GlobalsOptions;
import tdg.utils.Functions;
import tdg.utils.GeneticCode;

import java.io.BufferedReader;
import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class SingleSiteSimulator {

    public static void main(String[] args) throws Exception {
        SingleSiteSimulator sss = new SingleSiteSimulator();
        JCommander jc = new JCommander(sss);
        jc.setProgramName("java -cp " + Constants.PROGRAM_JAR + "tdg.sim.Simulator ");

        try {
            //jc.parse("-tree", "/Users/atamuri/Documents/work/tdg12/etc/PB2_FMutSel0.tree.out", "-sites", "10", "-heteroclades", "Av,Hu", "-fitness", "1,2,3,4,5,6,7,8", "-fitness", "1,2,3,4,5,6,7,8", "-characters", "A,R,N,D,C,Q,E,H", "-tau", "1e-2", "-kappa", "7.5", "-pi", "0.25,0.25,0.25,0.25", "-mu", "2.0", "-gc", "standard");
            jc.parse("-tree", "/Users/atamuri/Documents/work/tdg12/etc/PB2_FMutSel0.tree.out", "-heteroclades", "Av,Hu",
                    "-fitnessfile", "/Users/atamuri/Documents/work/mitochondria/paper/response/pb2.no.pen/results.nonhomog/fitness.av.sorted.txt",
                    "-fitnessfile", "/Users/atamuri/Documents/work/mitochondria/paper/response/pb2.no.pen/results.nonhomog/fitness.hu.sorted.txt",
                    "-tau", "1e-2", "-kappa", "7.5", "-pi", "0.25,0.25,0.25,0.25", "-mu", "2.0", "-gc", "standard");
            //jc.parse("-tree", "/Users/atamuri/Documents/work/tdg12/etc/PB2_FMutSel0.tree", "-sites", "10", "-fitness", "1,2,3,4,5,6,7,8", "-characters", "A,R,N,D,C,Q,E,H", "-tau", "1e-2", "-kappa", "7.5", "-pi", "0.25,0.25,0.25,0.25", "-mu", "2.0", "-gc", "standard");
        } catch (ParameterException pe) {
            System.out.printf("Error: %s\n\n", pe.getMessage());
            jc.usage();
            System.exit(0);
        }

        // fitnesses should have same size as characters (homogenous model) or size of characters * heteroclades (heterogeneous model)
        sss.validate();
        sss.run();

    }

    private void run() throws Exception {
        Simulator s = new Simulator();


        // If were simulating a single set of fitnesses
        if (this.fitness.size() > 0) {
            s.initialise(tree, sites, globalOptions);
            s.setClades(heteroClades);
            s.setAminoAcids(residues);
            for (int i = 0; i < heteroClades.size(); i++) {
                s.setCladeModel(heteroClades.get(i), fitness.subList(i * residues.length, (i + 1) * residues.length));
            }

            s.simulate();

            s.writeSimulatedData();
        } else {
            s.initialise(tree, 1, globalOptions);
            s.setClades(heteroClades);
            s.setAminoAcids(residues);

            Map<String, BufferedReader> cladeFitnessReaders = Maps.newHashMap();

            for (int i = 0; i < heteroClades.size(); i++) {
                cladeFitnessReaders.put(heteroClades.get(i), Files.newReader(new File(this.fitnessFiles.get(i)), Charsets.US_ASCII));
            }

            // Collect each simulated site here
            // Map<String, int[]> allSimulatedSites = Maps.newHashMap();
            ArrayListMultimap<String, Integer> allSimulatedSites = ArrayListMultimap.create();
            Map<String, int[]> thisSite;

            String line;
            while ((line = cladeFitnessReaders.get(heteroClades.get(0)).readLine()) != null) {
                // line is the fitness for the first clade
                s.setCladeModel(heteroClades.get(0), Lists.transform(Arrays.asList(line.split(" ")), Functions.stringToDouble()));

                for (int i = 1; i < heteroClades.size(); i++) {
                    s.setCladeModel(heteroClades.get(i), Lists.transform(Arrays.asList(cladeFitnessReaders.get(heteroClades.get(i)).readLine().split(" ")), Functions.stringToDouble()));
                }

                s.simulate();
                thisSite = s.getSimulatedData();

                for (Map.Entry<String, int[]> e : thisSite.entrySet()) {
                    allSimulatedSites.put(e.getKey(), e.getValue()[0]);
                }

            }

            // close all buffers
            for (BufferedReader br : cladeFitnessReaders.values()) {
                br.close();
            }

            // print the full simulated data
            Map<String, Collection<Integer>> c = allSimulatedSites.asMap();
            for (Map.Entry<String, Collection<Integer>> e : c.entrySet()) {
                List<String> codons = Lists.transform(allSimulatedSites.get(e.getKey()), new Function<Integer, String>() {
                    @Override
                    public String apply(Integer o) {
                        return GeneticCode.getInstance().getCodonTLA(o);
                    }
                });
                System.out.printf("%s     %s\n", e.getKey(), Joiner.on("").join(codons));
            }
            System.out.println();
        }


    }

    public void validate() {
        if (this.fitness.size() == 0 && this.fitnessFiles.size() == 0) {
            throw new RuntimeException("You must specify either -fitness or -fitnessfile.");
        }

        // We can't have both -fitness and -fitnessfile specified
        if (this.fitness.size() > 0 && this.fitnessFiles.size() > 0) {
            throw new RuntimeException("You can only use -fitness or -fitnessfile, not both together.");
        }


        // If user has specified fitness using the -fitness parameter
        if (this.fitness.size() > 0) {
            // If we're running the homogeneous model
            if (this.heteroClades.size() == 1) {
                // We expect the same number of fitnesses to characters.
                if (this.fitness.size() != this.residues.length) {
                    throw new RuntimeException(
                            String.format("You have %s fitnesses but %s characters. They should be the same number for each.\n",
                                    this.fitness.size(), this.residues.length));
                }
                // Otherwise, we have a heterogeneous model
            } else {
                // We expect the number of fitnesses to be (number of characters * number of clades)
                if (this.fitness.size() != (this.residues.length * this.heteroClades.size())) {
                    throw new RuntimeException(
                            String.format("You have %s fitnesses for %s clades and %s characters. You should have %s fitnesses (%s for each clade).\n",
                                    this.fitness.size(), this.heteroClades.size(), this.residues.length, (this.residues.length * this.heteroClades.size()), this.residues.length));
                }
            }
        } else {
            // we're using fitness files - check that we have the right number of files
            if (this.heteroClades.size() != this.fitnessFiles.size()) {
                throw new RuntimeException("You have defined " + this.heteroClades.size() + " clades but have " + this.fitnessFiles.size() + " fitness file(s).");
            }
        }
    }

    // COMMAND-LINE ARGUMENTS (for JCommander)
    @Parameter(names = "-tree", description = "Tree file in Newick format", required = true)
    public String tree;

    @Parameter(names = "-sites", description = "Number of times to simulate each set of fitnesses.", required = false)
    public int sites = 1;

    // -fitness option can be specified multiple times for heterogeneous models
    @Parameter(names = "-fitness", description = "Comma-separated fitness coefficients.", required = false)
    public List<Double> fitness = Lists.newArrayList();

    // -fitnessfile option can be specified multiple times for heterogenous models
    @Parameter(names = "-fitnessfile", description = "A file containing space-separated fitness coefficients, one row per site.", required = false)
    public List<String> fitnessFiles = Lists.newArrayList();

    // The fitnesses in -fitness can be specified in any order. Default is canonical amino acid order
    @Parameter(names = "-characters", description = "Comma-separated amino acids (matching -fitness order).", converter = CharArrayConverter.class, required = false)
    public char[] residues = GeneticCode.AMINO_ACIDS;

    // You must specify this for heterogeneous models
    @Parameter(names = "-heteroclades", description = "If simulating heterogeneous model, specify the clade labels.", required = false)
    public List<String> heteroClades = Lists.newArrayList("ALL");

    @ParametersDelegate
    public GlobalsOptions globalOptions = new GlobalsOptions();

    @ParametersDelegate
    public GeneticCodeOption gc = new GeneticCodeOption();
}
