package tdg.utils;

import com.google.common.primitives.Chars;
import com.google.common.primitives.Ints;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * Copied some of the Nucleotide and Codon functions we're using from PAL because
 * PAL uses { A, C, G, T } for nucleotide states, but we want { T, C, A, G } !!
 *
 * From http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi:
 *
 * Standard Code:
 *     AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
 *   Starts = ---M---------------M---------------M----------------------------
 *   Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
 *   Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
 *   Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 *
 * Vertebrate Mitochondrial Code:
 *     AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG
 *   Starts = --------------------------------MMMM---------------M------------
 *   Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
 *   Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
 *   Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 *
 * @author Asif Tamuri
 * @version $Id: GeneticCode.java 151 2010-11-08 11:08:36Z tamuri $
 */
public class GeneticCode {
    public static final int UNKNOWN_STATE = -1;
    private static final char UNKNOWN_CHARACTER = '?';
    private static final String UNKNOWN_TLA = "???";
    private static final char[] NUCLEOTIDES = "TCAG".toCharArray();
    private static final int UT_STATE = 0;
    private static final int C_STATE  = 1;
    private static final int A_STATE  = 2;
    private static final int G_STATE  = 3;




    private static final char[] AMINO_ACIDS = "ARNDCQEGHILKMFPSTWYV".toCharArray();

    // From http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    public static final String STANDARD_CODE                 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    public static final String VERTEBRATE_MITOCHONDRIAL_CODE = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG";
    
    public static String CURRENT_CODE;

    public static final int CODON_STATES = 64;
    public static final int AMINO_ACID_STATES = AMINO_ACIDS.length;

    private static final char[][] CODONS_TLA_CHAR_ARRAY = new char[CODON_STATES][CODON_STATES];

    static {
        int i = 0;
        for (char n1 : NUCLEOTIDES) {
            for (char n2 : NUCLEOTIDES) {
                for (char n3 : NUCLEOTIDES) {
                    CODONS_TLA_CHAR_ARRAY[i++] = new char[]{n1, n2, n3};
                }
            }
        }
    }

    private static final String[] CODONS_TLA = new String[CODON_STATES];
    static {
        int i = 0;
        for (char[] c : CODONS_TLA_CHAR_ARRAY) {
            CODONS_TLA[i++] = new String(c);
        }
    }

    // TODO: Stop codons should map to a 'stop' amino acid, 20, rather than -1!
    // so they can be assigned a fitness for analyses (MdR)
    private static final int[] CODONS_TO_AMINO_ACIDS = new int[CODON_STATES];
    private static final int[][] AMINO_ACIDS_TO_CODONS = new int[AMINO_ACID_STATES][];

    public synchronized static void initialise(String code) {
        if (initialised)
            throw new RuntimeException("Genetic code already initialised!");

        GeneticCode gc = new GeneticCode();
        CURRENT_CODE = code;

        for (int i = 0; i < CODON_STATES; i++) {
            if (code.charAt(i) == '*') {
                CODONS_TO_AMINO_ACIDS[i] = UNKNOWN_STATE; // stop codon
            } else {
                CODONS_TO_AMINO_ACIDS[i] = gc.getAminoAcidIndexByChar(code.charAt(i));
            }
        }

        for (int i = 0; i < AMINO_ACID_STATES; i++) {
            List<Integer> l = new ArrayList<Integer>();
            int pos = 0;
            while (true) {
                pos = code.indexOf(AMINO_ACIDS[i], pos);
                if (pos == -1) {
                    break;
                } else {
                    l.add(pos);
                }
                pos++;
            }
            AMINO_ACIDS_TO_CODONS[i] = Ints.toArray(l);
        }

        System.out.printf("tdg.utils.GeneticCode - Initialised with %s\n", code.equals(VERTEBRATE_MITOCHONDRIAL_CODE) ? "VERTEBRATE_MITOCHONDRIAL_CODE" : "STANDARD_CODE");

        System.out.printf("tdg.utils.GeneticCode - Stop codons: ");
        for (int i = 0; i < CODON_STATES; i++) if (gc.isUnknownAminoAcidState(gc.getAminoAcidIndexFromCodonIndex(i))) System.out.printf("%s ", gc.getCodonTLA(i));
        System.out.println();
        System.out.printf("tdg.utils.GeneticCode - Methionine codon(s): ");
        for (int i : gc.getCodonIndexFromAminoAcidIndex(gc.getAminoAcidIndexByChar('M'))) System.out.printf("%s ", gc.getCodonTLA(i));
        System.out.println();

        initialised = true;
        INSTANCE = gc;
    }

    private GeneticCode(){}

    public int getNucleotideIndexByChar(char c) {
        return Chars.indexOf(NUCLEOTIDES, c);
    }

    public int getAminoAcidIndexByChar(char c) {
        return Chars.indexOf(AMINO_ACIDS, c);
    }

    public char getNucleotideCharByIndex(int i) {
        return NUCLEOTIDES[i];
    }

    public char getAminoAcidCharByIndex(int i) {
        if (isUnknownAminoAcidState(i)) return UNKNOWN_CHARACTER;
        return AMINO_ACIDS[i];
    }

    public boolean isUnknownCodonState(int i) {
        return (i < 0 || i >= CODON_STATES);
    }

    public String getCodonTLA(int i) {
        if (isUnknownCodonState(i)) return UNKNOWN_TLA;
        return CODONS_TLA[i];
    }

    public char[] getNucleotidesFromCodonIndex(int i) {
        if (isUnknownCodonState(i)) return UNKNOWN_TLA.toCharArray();
        return CODONS_TLA_CHAR_ARRAY[i];
    }

    public int getCodonIndexFromNucleotides(char[] nuc) {
        String n = new String(nuc);
        for (int i = 0; i < CODONS_TLA.length; i++) {
            if (CODONS_TLA[i].equals(n)) {
                return i;
            }
        }
        return UNKNOWN_STATE;
    }

    public int getAminoAcidIndexFromCodonIndex(int i) {
        if (isUnknownCodonState(i)) return UNKNOWN_STATE;
        return CODONS_TO_AMINO_ACIDS[i];
    }

    public int[] getCodonIndexFromAminoAcidIndex(int i) {
        if (isUnknownAminoAcidState(i)) return new int[]{UNKNOWN_STATE};
        return AMINO_ACIDS_TO_CODONS[i];
    }

    public boolean isUnknownAminoAcidState(int i) {
        return (i < 0 || i >= AMINO_ACID_STATES);
    }

    public boolean isTransitionByChar(char firstChar, char secondChar) {
        return isTransitionByState(getNucleotideIndexByChar(firstChar), getNucleotideIndexByChar(secondChar));
    }

    public boolean isTransitionByState(int firstState, int secondState) {
        switch(firstState) {
            case A_STATE: {
                if(secondState==G_STATE) {
                    return true;
                }
                return false;
            }
            case C_STATE : {
                if(secondState==UT_STATE) {
                    return true;
                }
                return false;
            }
            case G_STATE : {
                if(secondState==A_STATE) {
                    return true;
                }
                return false;
            }
            case UT_STATE : {
                if(secondState==C_STATE) {
                    return true;
                }
                return false;
            }
            default:
                return false;
        }
    }

    private static GeneticCode INSTANCE;
    private static boolean initialised = false;

    public static GeneticCode getInstance() {
        //  if (!initialised) throw new RuntimeException("Must initialise GeneticCode!");
        return INSTANCE;
    }
}
