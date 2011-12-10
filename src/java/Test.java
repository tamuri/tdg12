import com.google.common.primitives.Ints;
import tdg.models.TDGGlobals;
import tdg.utils.GeneticCode;

public class Test {
    public static void main(String[] args) {
        GeneticCode.initialise(GeneticCode.STANDARD_CODE);
        /*

-kappa
7.06084250000000e+00
-pi
1.85704375000000e-01,2.73384375000000e-01,4.78586875000000e-01,6.23243750000000e-02
-tau
1.25000000000000e-02
-mu
2.41435441700000e+00
-gc
vertebrate_mit

         */
        /*TDGGlobals g = new TDGGlobals(
                1.25000000000000e-02,
                7.06084250000000e+00,
                new double[]{1.85704375000000e-01, 2.73384375000000e-01, 4.78586875000000e-01, 6.23243750000000e-02},
                2.41435441700000e+00,
                0.0);*/

        /*

-tau
1.25010000000000e-02
-kappa
7.86404250000000e+00
-pi
2.36114375000000e-01,1.95774375000000e-01,3.65944375000000e-01,2.02166875000000e-01
-mu
3.13262666700000e+00
-gc
standard

        
         */
TDGGlobals g = new TDGGlobals(
                1.25000000000000e-02,
                7.86404250000000e+00,
                new double[]{2.36114375000000e-01,1.95774375000000e-01,3.65944375000000e-01,2.02166875000000e-01},
                3.13262666700000e+00,
                0.0);

        System.out.printf("%s\n", g.getNu());
        System.out.printf("%s / %s = %s\n", g.getTau() , g.getNu(), g.getTau() / g.getNu());
    }
    public static void main2(String[] args) {
        double[] fitnesses = new double[]{ 0.0, -2.973826358478222, -0.8495727798560446, -2.9473365716081483, -2.8783469421717935, -2.905381086132148, -3.5213452784524204, -2.4006269186728693, -3.76861482080034, -3.508493905830959, -1.338012571844407, -3.3016351422513743, -2.807265662187931, 0.9050451875512675 };
        int[] lookup = new int[]{6,11,15,16,10,5,14,0,12, 2,7, 9, 19, 4};
        double[] pi = new double[]{1.83840000000000e-01, 2.70890000000000e-01, 4.83530000000000e-01, 6.17400000000000e-02};
        double[] codonPi = new double[64];
        double z = 0;

        GeneticCode.initialise(GeneticCode.VERTEBRATE_MITOCHONDRIAL_CODE);

        for (int i = 0; i < 64; i++) {
            int aa = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(i);

            char[] nucs = GeneticCode.getInstance().getNucleotidesFromCodonIndex(i);

            if (GeneticCode.getInstance().isUnknownAminoAcidState(aa)) continue;

            if (Ints.contains(lookup, aa)) {
                double thisCodonFitness = fitnesses[Ints.indexOf(lookup, aa)];

                codonPi[i] = pi[GeneticCode.getInstance().getNucleotideIndexByChar(nucs[0])] *
                        pi[GeneticCode.getInstance().getNucleotideIndexByChar(nucs[1])] *
                        pi[GeneticCode.getInstance().getNucleotideIndexByChar(nucs[2])] *
                        Math.exp(thisCodonFitness);

                z += codonPi[i];
            }
        }

        // normalize


        
        for (int i = 0; i < 64 ; i++) {
            codonPi[i] /= z;
            System.out.printf("%s/%s -> %s\n", GeneticCode.getInstance().getCodonTLA(i),
                    GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(i) >= 0 ? GeneticCode.getInstance().getAminoAcidCharByIndex(GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(i)) : "STOP",
                    codonPi[i]);
        }

        // System.out.printf("%s\n",         Doubles.join("\n", codonPi));


        double[] freqs = new double[20];

        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            if (GeneticCode.getInstance().isUnknownAminoAcidState(GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(i))) continue;
            freqs[GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(i)] += codonPi[i];
        }

        freqs = new double[20];

        for (int i = 0; i < GeneticCode.AMINO_ACID_STATES; i++) {
            System.out.printf("%s/%s -> ", i, GeneticCode.getInstance().getAminoAcidCharByIndex(i));
            int[] codonIndexes = GeneticCode.getInstance().getCodonIndexFromAminoAcidIndex(i);


            for (int j = 0; j < codonIndexes.length; j++) {
                System.out.printf("%s/%s%s", codonIndexes[j], GeneticCode.getInstance().getCodonTLA(codonIndexes[j]), j < codonIndexes.length - 1 ? "," : "");
                freqs[i] += codonPi[codonIndexes[j]];
            }
            System.out.print(" -> ");
            for (int j = 0; j < codonIndexes.length; j++) {
                System.out.printf("%s%s", codonPi[codonIndexes[j]], j < codonIndexes.length - 1 ? " + " : " = ");
            }

            System.out.printf("%s\n", freqs[i]);


         }


        for (int i = 0; i < GeneticCode.AMINO_ACID_STATES; i++) {
            System.out.printf("%s/%s -> %s \n", i, GeneticCode.getInstance().getAminoAcidCharByIndex(i), freqs[i]);
        }
    }
}
