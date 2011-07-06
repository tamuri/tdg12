package tdg.analysis;

import com.google.common.base.Function;
import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Doubles;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.MaxEvaluationsExceededException;
import org.apache.commons.math.optimization.GoalType;
import org.apache.commons.math.optimization.RealConvergenceChecker;
import org.apache.commons.math.optimization.RealPointValuePair;
import org.apache.commons.math.optimization.SimpleScalarValueChecker;
import org.apache.commons.math.optimization.direct.DirectSearchOptimizer;
import org.apache.commons.math.optimization.direct.MultiDirectional;
import org.apache.commons.math.optimization.direct.NelderMead;
import org.apache.commons.math.random.MersenneTwister;
import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;
import pal.alignment.Alignment;
import pal.tree.Tree;
import tdg.Options;
import tdg.models.LikelihoodCalculator;
import tdg.models.MinimisationParameters;
import tdg.models.TDGCodonModel;
import tdg.models.TDGGlobals;
import tdg.models.parameters.Fitness;
import tdg.models.parameters.Parameter;
import tdg.optim.EquivalentValueConvergenceChecker;
import tdg.optim.LikelihoodMaximiser;
import tdg.optim.PointAndValueConvergenceChecker;
import tdg.utils.GeneticCode;
import tdg.utils.PhyloUtils;

import javax.annotation.Generated;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Asif Tamuri
 * @version $Id: SiteAnalyser.java 157 2010-12-02 16:09:50Z atamuri $
 */
public class SiteAnalyser {
    private static final int INITIAL_PARAM_RANGE = 3;
    private double homogeneousLikelihood;
    private double nonHomogeneousLikelihood;
    private final Tree tree;
    private final Alignment alignment;
    private final TDGGlobals globals;
    private final int site;
    private final Options options;

    public SiteAnalyser(Tree tree, Alignment alignment, TDGGlobals globals, int site, Options options) {
        this.tree = tree;
        this.alignment = alignment;
        this.globals = globals;
        this.site = site;
        this.options = options;
    }

    public void run() {
        long startTime = System.currentTimeMillis();

        Map<String, Integer> sitePattern = PhyloUtils.getCodonsAtSite(alignment, site);

        // Remove any stop codons
        for (Map.Entry<String, Integer> e : sitePattern.entrySet()) {
            if (GeneticCode.getInstance().isUnknownCodonState(GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(e.getValue()))) {
                System.out.printf("Site %s - Sequence %s has stop codon (%s) - removing.\n", site, e.getKey(), GeneticCode.getInstance().getCodonTLA(e.getValue()));
                sitePattern.put(e.getKey(), -1);
            }
        }

        List<Integer> aminoAcidsAtSite = PhyloUtils.getDistinctAminoAcids(sitePattern.values());
        int observedResidueCount = aminoAcidsAtSite.size();

        // If we're not using the approximate method (collapsing the matrix)
        if (!options.approx) {
            // Optimise all 19 Fitness parameters
            for (int i = 0; i < GeneticCode.AMINO_ACID_STATES; i++) {
                if (!aminoAcidsAtSite.contains(i)) aminoAcidsAtSite.add(i);
            }
        }

        /* TODO: Default fitness
        if (aminoAcidsAtSite.size() > 1) {

            for (int i = 0; i < aminoAcidsAtSite.size(); i++) {
                f[i] = Defaults.getDefaultFitness(site, aminoAcidsAtSite.get(i));
            }
            // Fitness homogeneousFitness = new Fitness(new double[aminoAcidsAtSite.size()], true);
            homogeneousFitness = new Fitness(f, true);
            System.out.printf("Site %s - Initial fitness: { %s }\n", site, Doubles.join(", ", f));
        }
        */

        // ********************* HOMOGENEOUS MODEL *****************************
        LikelihoodCalculator homogeneousModel = new LikelihoodCalculator(tree, sitePattern);

        double[] f = new double[aminoAcidsAtSite.size()];
        Fitness homogeneousFitness = new Fitness(f, true);

        List<String> aaChar = Lists.transform(aminoAcidsAtSite, new Function<Integer,String>(){
            int pos = 1;
            public String apply(Integer input) {
                return String.format("%s:(%s, %s)", pos++, input, GeneticCode.getInstance().getAminoAcidCharByIndex(input));
            }
        });

        System.out.printf("Site %s - Residues: [%s/%s] { %s }\n", site, observedResidueCount, aminoAcidsAtSite.size(), Joiner.on(", ").join(aaChar));

        RandomData randomData = new RandomDataImpl(new MersenneTwister());

        // This makes a big difference to total optimisation time. We should think about which
        // convergence checked to use for optimising global parameters vs. site fitness parameters
        // Same as the standard scalar value convergence checker, but regards function has converged if
        // 50 consecutive evaluations do not change by more than 1E-6.
        //RealConvergenceChecker convergenceChecker = new EquivalentValueConvergenceChecker(1E-6, 50);
        //RealConvergenceChecker convergenceChecker = new PointAndValueConvergenceChecker(1E-6, 20, 1E-6);
        RealConvergenceChecker convergenceChecker = new SimpleScalarValueChecker(-1, 1E-6);
        
        int runs = options.optimRuns;

        Map<String, RealPointValuePair> optimiseRuns = Maps.newHashMap();


        TDGCodonModel tcm1 = new TDGCodonModel(globals, homogeneousFitness, aminoAcidsAtSite);
        homogeneousModel.addCladeModel("ALL", tcm1);

        for (int i = 0; i < runs; i++) {
            // The purpose of multiple runs is to test the optimisation routine with different initial parameters
            // So if we are performing multiple runs, then pick a set of random initial parameters for this run
            if (options.optimRuns > 1) {
                /* // this is the original...
                double[] newF = new double[aminoAcidsAtSite.size()];
                for (int j = 1; j < aminoAcidsAtSite.size(); j++) newF[j] = randomData.nextUniform(-INITIAL_PARAM_RANGE, INITIAL_PARAM_RANGE);
                homogeneousFitness.set(newF); //
                */
                // we're testing approx vs full:
                double[] newF = new double[aminoAcidsAtSite.size()];
                for (int j = 1; j < aminoAcidsAtSite.size(); j++) {
                    if (i == 0) {
                        // first run - all 0
                        newF[j] = 0;

                    } else if (i == 1) {
                        // second run - observed 0, unobserved -20
                        if (j < observedResidueCount) {
                            newF[j] = 0;
                        } else {
                            newF[j] = -20;
                        }
                    } else if (i == 2) {
                        // third run - random observed, unobserved -20
                        if (j < observedResidueCount) {
                            newF[j] = randomData.nextUniform(-INITIAL_PARAM_RANGE, INITIAL_PARAM_RANGE);
                        }

                        else {
                            newF[j] = -20;

                        }
                        // third run - all random
                        //newF[j] = randomData.nextUniform(-INITIAL_PARAM_RANGE, INITIAL_PARAM_RANGE);
                    }

                }
                //System.out.printf("newF = %s\n", Doubles.join(", ", newF));
                homogeneousFitness.set(newF);
            }

            homogeneousModel.setParameters(homogeneousFitness);

            // TODO: better way to handle this!!!?
            homogeneousModel.applyErrorMatrix(globals.getErrorMatrix());

            // Single residue observed at this site
            if (aminoAcidsAtSite.size() == 1) {
                System.out.printf("Site %s - Site is conserved.\n", site);
                homogeneousLikelihood = -homogeneousModel.function(new double[]{});

                if (!options.homogonly) {
                    nonHomogeneousLikelihood = homogeneousLikelihood;
                }

                System.out.printf("Site %s - Homogeneous model lnL: %s\n", site, homogeneousLikelihood);
                System.out.printf("Site %s - Fitness: { %s }\n", site, Doubles.join(", ", homogeneousFitness.get()));
                System.out.printf("Site %s - Pi: { %s }\n", site, Doubles.join(", ", tcm1.getAminoAcidFrequencies()));

                //TODO: we exit out of method here...what about the rest of the output??
                return;
            }
    
            MinimisationParameters mp = homogeneousModel.getMinimisationParameters();
            DirectSearchOptimizer dso = new NelderMead(); //new MultiDirectional(); //
            // Default convergence criteria : relative: 1.972152e-29, absolute: < 1.780059e-305
            dso.setConvergenceChecker(convergenceChecker);
            dso.setMaxEvaluations(10000);
            LikelihoodMaximiser maximiser = new LikelihoodMaximiser();
            maximiser.setLc(homogeneousModel);

            // System.out.printf("Site %s - Starting optimisation...", site);
            RealPointValuePair r;
            try {
                r = dso.optimize(maximiser, GoalType.MINIMIZE, mp.getParameters());
                String key = String.format("%.3f", -r.getValue());
                if (!optimiseRuns.containsKey(key)) optimiseRuns.put(key, r);
                System.out.printf("Site %s - Optimisation run %s (%s evaluations). lnL = %s, Params = {%s -> %s}\n", site, i + 1, dso.getEvaluations(), -r.getValue(), Doubles.join(", ", mp.getParameters()), Doubles.join(", ", r.getPoint()));
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
        }


        RealPointValuePair r;
        if (optimiseRuns.size() == 1) {
            r = optimiseRuns.values().iterator().next();
            if (runs > 1) System.out.printf("Site %s - %s runs converged to the same optima. (%s)\n", site, runs, -r.getValue());
        } else {
            System.out.printf("Site %s - %s runs converged to %s different optima. (%s)\n", site, runs, optimiseRuns.size(), Joiner.on(", ").join(optimiseRuns.keySet()));
            // get the best
            double best = Double.NEGATIVE_INFINITY;
            String bestKey = "";
            for (String k : optimiseRuns.keySet()) {
                double val = Double.parseDouble(k);
                if (val > best) {
                    bestKey = k;
                    best = val;
                }
            }

            r = optimiseRuns.get(bestKey);
        }

        homogeneousModel.function(r.getPoint());
        System.out.printf("Site %s - Homogeneous model lnL: %s\n", site, -r.getValue());

        // if we're optimising all 19 fitness parameters
        double[] orderedFitness = new double[GeneticCode.AMINO_ACID_STATES];
        Arrays.fill(orderedFitness, Double.NEGATIVE_INFINITY);
        for (int i = 0; i < GeneticCode.AMINO_ACID_STATES; i++) {
            if (aminoAcidsAtSite.contains(i)) orderedFitness[i] = homogeneousFitness.get()[aminoAcidsAtSite.indexOf(i)];
        }

        //System.out.printf("Site %s - Fitness: { %s }\n", site, Doubles.join(", ", homogeneousFitness.get()));
        System.out.printf("Site %s - Fitness: { %s }\n", site, Doubles.join(", ", orderedFitness).replaceAll("Infinity", "Inf"));

        System.out.printf("Site %s - Pi: { %s }\n", site, Doubles.join(", ", tcm1.getAminoAcidFrequencies()));
        homogeneousLikelihood = -r.getValue();

        /* // To print JSON string for fitness
        StringBuffer sb = new StringBuffer();
        sb.append(String.format("Site %s - JSON: {\"%s\":{", site, site));
        for (int i = 0; i < aminoAcidsAtSite.size(); i++) {
            sb.append(String.format("\"%s\":%s", aminoAcidsAtSite.get(i), homogeneousFitness.get()[i]));
            if (i < aminoAcidsAtSite.size() - 1) sb.append(String.format(", "));
        }
        sb.append(String.format("}}\n"));
        System.out.printf("%s", sb.toString());
        */


        if (options.homogonly) {
            return;
        }

        // ********************* NON-HOMOGENEOUS MODEL *************************
        LikelihoodCalculator nonHomogeneousModel = new LikelihoodCalculator(tree, sitePattern);

        // TODO: we should be reading the list of clade labels from command-line Options!
        List<String> clades = Lists.newArrayList("Av", "Hu", "Sw");
        List<Fitness> fitnesses = Lists.newArrayListWithCapacity(clades.size());
        List<TDGCodonModel> tdgModels = Lists.newArrayListWithCapacity(clades.size());
        for (int i = 0; i < clades.size(); i++) {
            fitnesses.add(i, new Fitness(homogeneousFitness.get().clone(), true));
            tdgModels.add(i, new TDGCodonModel(globals, fitnesses.get(i), aminoAcidsAtSite));
            nonHomogeneousModel.addCladeModel(clades.get(i), tdgModels.get(i));
        }

        /*nonHomogeneousModel.addCladeModel("Av", tdgModels.get(0));
        nonHomogeneousModel.addCladeModel("Hu", tdgModels.get(0));
        nonHomogeneousModel.addCladeModel("Sw", tdgModels.get(1));
        */
        
        nonHomogeneousModel.setParameters(fitnesses.toArray(new Parameter[fitnesses.size()]));


        /*
        // TODO: For some sites, initial parameters matter! (Flat likelihood surface?) e.g. PB2 site 199. Try it:
        // Fitness nonHomogeneousFitnessAv = new Fitness(new double[aminoAcidsAtSite.size()], true);
        // Fitness nonHomogeneousFitnessHu = new Fitness(new double[aminoAcidsAtSite.size()], true);
        */

        MinimisationParameters mp2 = nonHomogeneousModel.getMinimisationParameters();

        DirectSearchOptimizer dso2 = new NelderMead();
        dso2.setConvergenceChecker(convergenceChecker);

        LikelihoodMaximiser maximiser2 = new LikelihoodMaximiser();
        maximiser2.setLc(nonHomogeneousModel);

        RealPointValuePair r2;
        try {
            r2 = dso2.optimize(maximiser2, GoalType.MINIMIZE, mp2.getParameters());
        } catch (Exception e) {
            throw new RuntimeException(e);
        }

        nonHomogeneousModel.function(r2.getPoint());
        System.out.printf("Site %s - Non-homogeneous model lnL: %s\n", site, -r2.getValue());
        for (int i = 0; i < clades.size(); i++) {

            orderedFitness = new double[GeneticCode.AMINO_ACID_STATES];
            Arrays.fill(orderedFitness, Double.NEGATIVE_INFINITY);
            for (int j = 0; j < GeneticCode.AMINO_ACID_STATES; j++) {
                if (aminoAcidsAtSite.contains(j)) orderedFitness[j] = fitnesses.get(i).get()[aminoAcidsAtSite.indexOf(j)];
            }

            //System.out.printf("Site %s - Fitness %s: { %s }\n", site, clades.get(i), Doubles.join(", ", fitnesses.get(i).get()));
            System.out.printf("Site %s - Fitness %s: { %s }\n", site, clades.get(i), Doubles.join(", ", orderedFitness).replaceAll("Infinity", "Inf"));
            System.out.printf("Site %s - Pi %s: { %s }\n", site, clades.get(i), Doubles.join(", ", tdgModels.get(i).getAminoAcidFrequencies()));
        }
        nonHomogeneousLikelihood = -r2.getValue();

        System.out.printf("Site %s - Time: %s ms\n", site, System.currentTimeMillis() - startTime);

    }

    public double getNonHomogeneousLikelihood() {
        return nonHomogeneousLikelihood;
    }

    public double getHomogeneousLikelihood() {
        return homogeneousLikelihood;
    }

}