package tdg;

import com.google.common.base.Function;
import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Doubles;
import org.apache.commons.math.optimization.GoalType;
import org.apache.commons.math.optimization.RealConvergenceChecker;
import org.apache.commons.math.optimization.RealPointValuePair;
import org.apache.commons.math.optimization.SimpleScalarValueChecker;
import org.apache.commons.math.optimization.direct.DirectSearchOptimizer;
import org.apache.commons.math.optimization.direct.NelderMead;
import org.apache.commons.math.random.MersenneTwister;
import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;
import pal.alignment.Alignment;
import pal.tree.Tree;
import tdg.cli.AnalyseOptions;
import tdg.model.*;
import tdg.optim.LikelihoodFunctionWrapper;
import tdg.utils.GeneticCode;
import tdg.utils.PhyloUtils;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * Performs the MLE for a single site, calling the likelihood function. (i.e. optimising the fitness parameters for
 * the swMutSel0 model.)
 * <p/>
 * TODO: Bit of an organic mess. Need to clean/refactor this, Fitness, LikelihoodCalculator, TDGCodonModel + LikelihoodFunctionWrapper...how??
 *
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 * @see LikelihoodCalculator
 * @see TDGCodonModel
 */
public class SiteAnalyser {
    private double homogeneousLikelihood;
    private double heterogeneousLikelihood;
    private final Tree tree;
    private final Alignment alignment;
    private final TDGGlobals globals;
    private final int site;
    private final AnalyseOptions options;
    private final RandomData randomData = new RandomDataImpl(new MersenneTwister());

    public SiteAnalyser(Tree tree, Alignment alignment, TDGGlobals globals, int site, AnalyseOptions options) {
        this.tree = tree;
        this.alignment = alignment;
        this.globals = globals;
        this.site = site;
        this.options = options;
    }

    public void run() {
        long startTime = System.currentTimeMillis();

        // Get the codons observed at this site
        Map<String, Integer> sitePattern = PhyloUtils.getCodonsAtSite(alignment, site);

        // Remove any stop codons and treat them as gaps
        for (Map.Entry<String, Integer> e : sitePattern.entrySet()) {
            if (GeneticCode.getInstance().isUnknownCodonState(GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(e.getValue()))) {
                System.out.printf("Site %s - Sequence %s has unknown/stop codon (%s). Treating as gap.\n",
                        site, e.getKey(), GeneticCode.getInstance().getCodonTLA(e.getValue()));
                sitePattern.put(e.getKey(), GeneticCode.UNKNOWN_STATE);
            }
        }

        // Get a list of all amino acids observed at this site (from the given codons)
        List<Integer> aminoAcidsAtSite = PhyloUtils.getDistinctAminoAcids(sitePattern.values());
        int observedResidueCount = aminoAcidsAtSite.size();

        // If we're not using the approximate method (collapsing the matrix)
        if (!options.approx.useapprox) {
            // Optimise all 19 Fitness parameters, rather than just the observed amino acids
            // This will add the remaining amino acids at the end of the list of observed residues in canonical order
            for (int i = 0; i < GeneticCode.AMINO_ACID_STATES; i++) {
                if (!aminoAcidsAtSite.contains(i)) aminoAcidsAtSite.add(i);
            }
        }

        // Display which amino acids we've observed at this site and for which we're going to estimate the fitness.
        displayResidues(aminoAcidsAtSite, observedResidueCount);

        // TODO: Use specified initial fitness for a site from e.g. a file. Would be good for global parameter optimisation.

        // ********************* HOMOGENEOUS MODEL *****************************
        LikelihoodCalculator homogeneousModel = new LikelihoodCalculator(tree, sitePattern, options);

        // An object for the fitnesses we want to estimate
        Fitness homogeneousFitness = new Fitness(new double[aminoAcidsAtSite.size()], true);

        // Create an instance of the swMutSel0 model, homogeneous model only has one model for the entire tree
        TDGCodonModel tcm1 = new TDGCodonModel(globals, homogeneousFitness, aminoAcidsAtSite);
        homogeneousModel.addCladeModel("ALL", tcm1);

        // The purpose of multiple runs is to ensure good convergence of optimisation by using different initial parameters.
        // We run the procedure as many times as specified in the optimRuns.
        int runs = options.optimRuns;
        Map<String, RealPointValuePair> optimiseRuns = Maps.newHashMap();

        for (int i = 0; i < runs; i++) {
            // Set the initial values for the fitness parameters (see method for details)
            homogeneousFitness.set(getInitialFitnessParameters(aminoAcidsAtSite, observedResidueCount, i));

            // The model will use the parameters in this Fitness object
            homogeneousModel.setParameters(homogeneousFitness);

            // TODO: Remove this (or make it an option)
            // homogeneousModel.applyErrorMatrix(globals.getErrorMatrix());

            // Single residue observed at this site
            if (aminoAcidsAtSite.size() == 1) {
                System.out.printf("Site %s - Site is conserved.\n", site);
                homogeneousLikelihood = homogeneousModel.function(new double[]{});

                if (options.heteroClades != null && options.heteroClades.length() > 0) {
                    heterogeneousLikelihood = homogeneousLikelihood;
                }

                System.out.printf("Site %s - Homogeneous model lnL: %s\n", site, homogeneousLikelihood);
                System.out.printf("Site %s - Fitness: { %s }\n", site, Doubles.join(", ", getOrderedFitness(aminoAcidsAtSite, homogeneousFitness.get())));
                System.out.printf("Site %s - Pi: { %s }\n", site, Doubles.join(", ", tcm1.getAminoAcidFrequencies()));

                //TODO: we exit out of method here...what about the rest of the output (e.g. heterogeneous model)?
                return;
            }

            RealPointValuePair r = optimise(homogeneousModel);

            // Store the log-likelihood, to 3 decimal places, and point for this run
            String key = String.format("%.3f", r.getValue());
            if (!optimiseRuns.containsKey(key)) optimiseRuns.put(key, r);

        }

        // We've completed the optimisation runs - see what we have, pick the best
        RealPointValuePair r;
        if (optimiseRuns.size() == 1) {
            r = optimiseRuns.values().iterator().next();
            if (runs > 1)
                System.out.printf("Site %s - %s runs converged to the same optima. (%s)\n", site, runs, r.getValue());
        } else {
            System.out.printf("Site %s - %s runs converged to %s different optima. (%s)\n", site, runs, optimiseRuns.size(), Joiner.on(", ").join(optimiseRuns.keySet()));
            // Get the best
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
        System.out.printf("Site %s - Homogeneous model lnL: %s\n", site, r.getValue());
        System.out.printf("Site %s - Fitness: { %s }\n", site, Doubles.join(", ", getOrderedFitness(aminoAcidsAtSite, homogeneousFitness.get())));
        System.out.printf("Site %s - Pi: { %s }\n", site, Doubles.join(", ", tcm1.getAminoAcidFrequencies()));
        homogeneousLikelihood = r.getValue();

        if (options.heteroClades == null) {
            System.out.printf("Site %s - Time: %s ms\n", site, System.currentTimeMillis() - startTime);
            return;
        }

        // ********************* HETEROGENEOUS MODEL *************************
        LikelihoodCalculator heterogeneousModel = new LikelihoodCalculator(tree, sitePattern, options);

        List<String> clades = Lists.newArrayList(options.heteroClades.split(","));
        List<Fitness> fitnesses = Lists.newArrayListWithCapacity(clades.size());
        List<TDGCodonModel> tdgModels = Lists.newArrayListWithCapacity(clades.size());
        for (int i = 0; i < clades.size(); i++) {
            fitnesses.add(i, new Fitness(homogeneousFitness.get().clone(), true));
            tdgModels.add(i, new TDGCodonModel(globals, fitnesses.get(i), aminoAcidsAtSite));
            heterogeneousModel.addCladeModel(clades.get(i), tdgModels.get(i));
        }
        heterogeneousModel.setParameters(fitnesses.toArray(new Parameter[fitnesses.size()]));

        RealPointValuePair r2 = optimise(heterogeneousModel);
        heterogeneousModel.function(r2.getPoint());

        System.out.printf("Site %s - Non-homogeneous model lnL: %s\n", site, r2.getValue());
        for (int i = 0; i < clades.size(); i++) {
            System.out.printf("Site %s - Fitness %s: { %s }\n", site, clades.get(i), Doubles.join(", ", getOrderedFitness(aminoAcidsAtSite, fitnesses.get(i).get())));
            System.out.printf("Site %s - Pi %s: { %s }\n", site, clades.get(i), Doubles.join(", ", tdgModels.get(i).getAminoAcidFrequencies()));
        }
        heterogeneousLikelihood = r2.getValue();

        System.out.printf("Site %s - Time: %s ms\n", site, System.currentTimeMillis() - startTime);
    }

    private RealPointValuePair optimise(LikelihoodCalculator model) {
        MinimisationParameters mp = model.getMinimisationParameters();
        DirectSearchOptimizer dso = new NelderMead();
        dso.setMaxEvaluations(Constants.MAX_EVALUATIONS);

        // Big difference for total optimisation time. We should think about which
        // convergence criteria to use for optimising global parameters vs. site fitness parameters
        // Same as the standard scalar value convergence checker, but regards function has converged if
        // 50 consecutive evaluations do not change by more than 1E-6.
        //RealConvergenceChecker convergenceChecker = new EquivalentValueConvergenceChecker(1E-6, 50);
        //RealConvergenceChecker convergenceChecker = new PointAndValueConvergenceChecker(1E-6, 20, 1E-6);
        RealConvergenceChecker convergenceChecker = new SimpleScalarValueChecker(-1, Constants.CONVERGENCE_TOL);

        dso.setConvergenceChecker(convergenceChecker);

        LikelihoodFunctionWrapper wrapper = new LikelihoodFunctionWrapper();
        wrapper.setLc(model);

        RealPointValuePair pair;
        try {
            pair = dso.optimize(wrapper, GoalType.MAXIMIZE, mp.getParameters());
            System.out.printf("Site %s - Optimisation run (%s evaluations). lnL = %s, params = { %s -> %s }\n",
                    site,
                    dso.getEvaluations(),
                    pair.getValue(),
                    Doubles.join(", ", mp.getParameters()), // initial parameters
                    Doubles.join(", ", pair.getPoint()));
        } catch (Exception e) {
            throw new RuntimeException(e);
        }

        return pair;
    }

    private double[] getOrderedFitness(List<Integer> aminoAcidsAtSite, double[] unorderedFitnesses) {
        double[] orderedFitness = new double[GeneticCode.AMINO_ACID_STATES];
        Arrays.fill(orderedFitness, Double.NEGATIVE_INFINITY);

        for (int j = 0; j < GeneticCode.AMINO_ACID_STATES; j++)
            if (aminoAcidsAtSite.contains(j))
                orderedFitness[j] = unorderedFitnesses[aminoAcidsAtSite.indexOf(j)];

        return orderedFitness;
    }

    private double[] getInitialFitnessParameters(List<Integer> aminoAcidsAtSite, int observedResidueCount, int run) {
        double[] initialFitness = new double[aminoAcidsAtSite.size()];
        for (int i = 1; i < aminoAcidsAtSite.size(); i++) {
            // First run - all residues have equal fitness = 0
            if (run == 0) {
                initialFitness[i] = 0;
                // Second run - observed residues have equal fitness (= 0), unobserved have equal fitness (= 20)
            } else if (run == 1 && !options.approx.useapprox) {
                if (i < observedResidueCount) initialFitness[i] = 0;
                else initialFitness[i] = -Constants.FITNESS_BOUND;
                // All other runs - all residues have random fitness picked from uniform distribution in range
            } else {
                initialFitness[i] = randomData.nextUniform(-Constants.RANDOM_INITIAL_FITNESS_RANGE, Constants.RANDOM_INITIAL_FITNESS_RANGE);
                // could also be random observed, unobserved -FITNESS_BOUND e.g.:
                // if (i < observedResidueCount) { initialFitness[i] = randomData.nextUniform(-RANDOM_INITIAL_FITNESS_RANGE, RANDOM_INITIAL_FITNESS_RANGE); else initialFitness[i] = -FITNESS_BOUND;
            }
        }
        return initialFitness;
    }

    private void displayResidues(List<Integer> aminoAcidsAtSite, int observedResidueCount) {
        List<String> residueStrings = Lists.transform(aminoAcidsAtSite, new Function<Integer, String>() {
            int pos = 1;

            public String apply(Integer input) {
                return String.format("%s:(%s, %s)", pos++, input, GeneticCode.getInstance().getAminoAcidCharByIndex(input));
            }
        });
        System.out.printf("Site %s - Residues: [%s/%s] { %s }\n", site, observedResidueCount, aminoAcidsAtSite.size(), Joiner.on(", ").join(residueStrings));
    }

    public double getHeterogeneousLikelihood() {
        return heterogeneousLikelihood;
    }

    public double getHomogeneousLikelihood() {
        return homogeneousLikelihood;
    }

}