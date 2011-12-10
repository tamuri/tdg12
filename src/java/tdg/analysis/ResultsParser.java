package tdg.analysis;


import com.google.common.base.Charsets;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import com.google.common.primitives.Doubles;
import org.apache.commons.lang.text.StrTokenizer;
import pal.alignment.Alignment;
import tdg.models.TDGCodonModel;
import tdg.models.TDGGlobals;
import tdg.models.parameters.Fitness;
import tdg.utils.GeneticCode;
import tdg.utils.PhyloUtils;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ResultsParser {
    public static void main(String[] args) {
        ResultsParser rp = new ResultsParser();
        rp.run();
    }

    private void run() {
        // PB2
/*
        String alignmentPath = "/Users/atamuri/Documents/work/tdg10/etc/PB2.co";
        String rawResultsPath = "/Users/atamuri/Documents/work/mitochondria/110304_TdG_PB2_RawResults.txt";
        double tau = 1.25010000000000e-02;
        double kappa = 7.86404250000000e+00;
        double[] pi = new double[]{2.36114375000000e-01, 1.95774375000000e-01, 3.65944375000000e-01, 2.02166875000000e-01};
        double mu = 3.13262666700000e+00;
        double gamma = 0;
        GeneticCode.initialise(GeneticCode.STANDARD_CODE);
*/


        // Mit

        String alignmentPath = "/Users/atamuri/Documents/work/tdg10/etc/all.but.ND6.AUT.phys";
//        String rawResultsPath = "/Users/atamuri/Documents/work/mitochondria/110307_TdG_Mit_RawResults.txt";
        String rawResultsPath = "/Users/atamuri/Documents/work/mitochondria/110318_TdG_Mit_ApproxVsFull/full.results.out";
        double tau = 1.25000000000000e-02;
        double kappa = 7.06084250000000e+00;
        double[] pi = new double[]{1.85704375000000e-01, 2.73384375000000e-01, 4.78586875000000e-01, 6.23243750000000e-02};
        double mu = 2.41435441700000e+00;
        double gamma = 0;
        GeneticCode.initialise(GeneticCode.VERTEBRATE_MITOCHONDRIAL_CODE);


        // START

        final TDGGlobals tdgGlobals = new TDGGlobals(tau, kappa, pi, mu, gamma);
        final Alignment alignment = PhyloUtils.readAlignment(alignmentPath);

        // Neutral codon pi
        double[] neutralCodonPi = new double[GeneticCode.CODON_STATES];
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) neutralCodonPi[i] = tdgGlobals.getCodonPiProduct(i);
        System.out.printf("Neutral codon pi:\n%s\n", Doubles.join(",", neutralCodonPi));

        // Neutral amino acid pi
        double[] neutralAminoAcidPi = new double[20];
        double z = 0.;
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            if (GeneticCode.getInstance().isUnknownAminoAcidState(GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(i))) continue;
            neutralAminoAcidPi[GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(i)] += neutralCodonPi[i];
            z += neutralCodonPi[i];
        }
        for (int i = 0; i < GeneticCode.AMINO_ACID_STATES; i++) neutralAminoAcidPi[i] /= z;
        System.out.printf("Neutral amino acid pi (STOP codons ignored):\n%s\n", Doubles.join(",", neutralAminoAcidPi));

        // Neutral Q
        double[] neutralMutationQ = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                if (i == j) neutralMutationQ[i * GeneticCode.CODON_STATES + j] = -1;
                neutralMutationQ[i * GeneticCode.CODON_STATES + j] = tdgGlobals.getNeutralMutationRate(i, j) * tdgGlobals.getNu();
            }
        }
        System.out.printf("v (scaling for q):\n%s\n", tdgGlobals.getNu());
        System.out.printf("Neutral q:\n%s\n", Doubles.join(",", neutralMutationQ));

        // Neutral synonymous and non-synonymous rates
        double neutralNonSynonymousRate = 0;
        double neutralSynonymousRate = 0;
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            for (int j = 0; j < GeneticCode.CODON_STATES; j++) {

                if (i == j) continue;

                int aa_i = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(i);
                int aa_j = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(j);

                double rate = neutralCodonPi[i] * neutralMutationQ[i * GeneticCode.CODON_STATES + j];

                if (aa_i != aa_j) {
                    neutralNonSynonymousRate += rate;
                } else {
                    neutralSynonymousRate += rate;
                }
            }
        }

        System.out.printf("Neutral non-synonymous rate:\n%s\n", neutralNonSynonymousRate);
        System.out.printf("Neutral synonymous rate:\n%s\n", neutralSynonymousRate);

        LineProcessor<double[]> lp = new ObjectLineProcessor<double[]>(alignment, tdgGlobals, neutralSynonymousRate, neutralNonSynonymousRate);

        try {
            Files.readLines(new File(rawResultsPath), Charsets.UTF_8, lp);
        } catch (Exception e) {
            throw new RuntimeException(e);
        }

        /*double[] out = lp.getResult();
        System.out.printf("w = %s / wN = %s / wS = %s\n", (1.0 / (alignment.getSiteCount() / 3) * out[0]),
                (1.0 / (alignment.getSiteCount() / 3) * out[1]),
                (1.0 / (alignment.getSiteCount() / 3) * out[2]));
        */

    }

    private static class ObjectLineProcessor<K> implements LineProcessor<double[]> {
        private int currentSite;
        private double sum_site_omega;
        private double sum_site_omegaN;
        private double sum_site_omegaS;
        private Map<String, Integer> sitePattern;
        private List<Integer> aminoAcidsAtSite;
        private final Alignment alignment;
        private final TDGGlobals tdgGlobals;
        private final double neutralSynonymousRate;
        private final double neutralNonSynonymousRate;

        public ObjectLineProcessor(Alignment alignment, TDGGlobals tdgGlobals, double neutralSynonymousRate, double neutralNonSynonymousRate) {
            this.alignment = alignment;
            this.tdgGlobals = tdgGlobals;
            this.neutralSynonymousRate = neutralSynonymousRate;
            this.neutralNonSynonymousRate = neutralNonSynonymousRate;
            currentSite = -1;
            sum_site_omega = 0;
            sum_site_omegaN = 0;
            sum_site_omegaS = 0;
        }

        @Override
        public boolean processLine(String line) throws IOException {
            String[] parts = line.split(" ");

            int site = Integer.parseInt(parts[1]);
            if (site != currentSite) {
                // we're looking at a new site
                System.out.println();
                currentSite = site;
                sitePattern = PhyloUtils.getCodonsAtSite(alignment, site);
                // Remove any stop codons
                for (Map.Entry<String, Integer> e : sitePattern.entrySet()) {
                    if (GeneticCode.getInstance().isUnknownAminoAcidState(GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(e.getValue())) ) {
                        //System.out.printf("Site %s - Sequence %s has stop/unknown codon (%s = %s) - removing.\n", site, e.getKey(), e.getValue(), GeneticCode.getInstance().getCodonTLA(e.getValue()));
                        sitePattern.put(e.getKey(), -1);
                    }
                }
                aminoAcidsAtSite = PhyloUtils.getDistinctAminoAcids(sitePattern.values());
                // Optimise all 19 Fitness parameters
                for (int i = 0; i < GeneticCode.AMINO_ACID_STATES; i++) {
                    if (!aminoAcidsAtSite.contains(i)) aminoAcidsAtSite.add(i);
                }
                //System.out.printf("%s,Residues,%s,", site, aminoAcidsAtSite.size());
            }

            String label = parts[3];

            if (label.startsWith("Fitness")) {

                // make sure we have the expected number of fitnesses, based on observed residues
                double[] fitness = parseFitnessString(line);
                if (fitness.length != aminoAcidsAtSite.size()) throw new RuntimeException();

                // now create the full fitness vector, full Q matrix, full S matrix, based on these fitnesses
                printModelParameters(fitness);

            } else if (label.equals("Homogeneous")) {

                //System.out.printf("lnL,%s,", Double.parseDouble(parts[6]));

            } else if (label.equals("Non-homogeneous")) {

                //System.out.printf("lnL,%s,", Double.parseDouble(parts[6]));

            } else if (label.equals("Residues:")) {

                Matcher m = Pattern.compile("(\\d+)").matcher(parts[4]);
                m.find();
//                if (aminoAcidsAtSite.size() != Integer.parseInt(m.group(0))) throw new RuntimeException();

            }

            return true;
        }

        private void printModelParameters(double[] f) {
            Fitness fitness = new Fitness(f, true);
            TDGCodonModel tcm1 = new TDGCodonModel(tdgGlobals, fitness, aminoAcidsAtSite);
            tcm1.updateModel();

            // full F
            double[] fullF = new double[20];
            for (int i = 0; i < GeneticCode.AMINO_ACID_STATES; i++) {
                if (aminoAcidsAtSite.contains(i)) {
                    fullF[i] = f[aminoAcidsAtSite.indexOf(i)];
                } else {
                    fullF[i] = -21;
                }
            }

//            System.out.printf("F,%s,", Doubles.join(",", fullF));
            System.out.printf("%s", Doubles.join(",", fullF));

            // full Q
            double[] fullQ = tcm1.getFullQ();
            for (int i = 0; i < GeneticCode.CODON_STATES; i++) fullQ[i * GeneticCode.CODON_STATES + i] = -1;
            //System.out.printf("QS,%s,", Doubles.join(",", fullQ));
            //System.out.printf("PiS,%s,", Doubles.join(",", tcm1.getCodonFrequencies()));

            // full S
            double[] fullS = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];
            for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
                for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                    if (i == j) {
                        fullS[i * GeneticCode.CODON_STATES + j] = 0;
                    } else {
                        int aa_from = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(i);
                        int aa_to = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(j);
                        if (aa_from < 0 || aa_to < 0) {
                            // to or from stop codons
                            fullS[i * GeneticCode.CODON_STATES + j] = -21;
                        } else if (aminoAcidsAtSite.contains(aa_from) && aminoAcidsAtSite.contains(aa_to)) {
                            // both amino acids that occur at this site
                            fullS[i * GeneticCode.CODON_STATES + j] = fullF[aa_to] - fullF[aa_from];
                        } else if (aminoAcidsAtSite.contains(aa_from) && !aminoAcidsAtSite.contains(aa_to)) {
                            // from observed to unobserved
                            fullS[i * GeneticCode.CODON_STATES + j] = -21;
                        } else if (!aminoAcidsAtSite.contains(aa_from) && aminoAcidsAtSite.contains(aa_to)) {
                            // from unobserved to observed
                            //System.out.printf("%s -> %s\n", aa_from, aa_to);
                            fullS[i * GeneticCode.CODON_STATES + j] = 21;
                        } else {
                            // neither of these amino acids are observed at this site TODO: SHIT! Check this!!
                            fullS[i * GeneticCode.CODON_STATES + j] = -21;
                        }
                    }
                }
            }
            //System.out.printf("S,%s,", Doubles.join(",", fullS));

            double siteNonSynonymousRate = 0;
            double siteSynonymousRate = 0;
            // omega
            for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
                for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                    if (i == j) continue;
                    int aa_i = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(i);
                    int aa_j = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(j);

                    // double rate = neutralCodonPi[i] * neutralMutationQ[i * GeneticCode.CODON_STATES + j];
                    double rate = tcm1.getCodonFrequencies()[i] * fullQ[i * GeneticCode.CODON_STATES + j];
                    // System.out.printf("rate = %s\n", rate);

                    if (aa_i != aa_j) {
                        siteNonSynonymousRate += rate;
                    } else {
                        siteSynonymousRate += rate;
                    }
                }
            }
            double site_omega = (siteNonSynonymousRate / siteSynonymousRate) * (neutralSynonymousRate / neutralNonSynonymousRate);
            double site_omega_wN = siteNonSynonymousRate / neutralNonSynonymousRate;
            double site_omega_wS = siteNonSynonymousRate / neutralSynonymousRate;
            //System.out.printf("w,%s,wN,%s,wS,%s,", site_omega, site_omega_wN, site_omega_wS);



            // PB2 only: if this is a conserved site, output everything again for non-homogeneous model
            /*if (aminoAcidsAtSite.size() == 1) {
                // Avian
                System.out.printf("lnL,%s,F,%s,QS,%s,PiS,%s,S,%s,w,%s,wN,%s,wS,%s,", 0, Doubles.join(",", fullF), Doubles.join(",", fullQ), Doubles.join(",", tcm1.getCodonFrequencies()), Doubles.join(",", fullS), site_omega, site_omega_wN, site_omega_wS);
                // Human
                System.out.printf("F,%s,QS,%s,PiS,%s,S,%s,w,%s,wN,%s,wS,%s,", Doubles.join(",", fullF), Doubles.join(",", fullQ), Doubles.join(",", tcm1.getCodonFrequencies()), Doubles.join(",", fullS), site_omega, site_omega_wN, site_omega_wS);
            }*/
           /* sum_site_omega += site_omega;
            sum_site_omegaN += site_omega_wN;
            sum_site_omegaS += site_omega_wS;*/
        }

        @Override
        public double[] getResult() {
            return new double[]{sum_site_omega, sum_site_omegaN, sum_site_omegaS};
        }

        private double[] parseFitnessString(String line) {
            Matcher m = Pattern.compile("\\{\\s(.*)\\s\\}").matcher(line);
            m.find();
            StrTokenizer st = StrTokenizer.getCSVInstance(m.group(1));
            double[] fitness = new double[st.size()];
            for (int i = 0; i < st.size(); i++) fitness[i] = Double.parseDouble(st.nextToken());
            return fitness;
        }
    }
}
