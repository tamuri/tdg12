package tdg.utils;

import com.google.common.base.Charsets;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Files;
import com.google.common.primitives.Ints;
import pal.alignment.Alignment;
import pal.alignment.AlignmentReaders;
import pal.datatype.DataTypeTool;
import pal.tree.ReadTree;
import pal.tree.SimpleTree;
import pal.tree.Tree;

import java.io.File;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * Collection of static methods, including some that wrap common functions available in PAL.
 *
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class PhyloUtils {
    /**
     * Reads a PHYLIP formatted alignment file
     * @param path Full path to alignment file
     * @return a PAL alignment object
     */
    public static Alignment readAlignment(String path) {
        Alignment a;
        try {
            a = AlignmentReaders.readPhylipClustalAlignment(Files.newReader(new File(path), Charsets.UTF_8), DataTypeTool.getNucleotides());
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
        return a;
    }

    /**
     * Returns a collection of codons observed at a particular site
     * @param alignment the PAL alignment object
     * @param site the location of interest
     * @return a map, where key for entry is taxon name and value is the codon state
     */
    public static Map<String, Integer> getCodonsAtSite(Alignment alignment, int site) {
        // Alignment holds a nucleotide alignment
        int nucleotideSite = (site - 1) * 3;
        Map<String, Integer> states = Maps.newHashMap();

        for (int i = 0; i < alignment.getSequenceCount(); i++) {

            /*
            // check for bad codons
            if (GeneticCode.getInstance().getCodonIndexFromNucleotides(
                            new char[]{alignment.getData(i, nucleotideSite),
                                    alignment.getData(i, nucleotideSite + 1),
                                    alignment.getData(i, nucleotideSite + 2)}
                    ) == -1) {
                System.out.printf("%s has %s  !!!\n", alignment.getIdentifier(i).getName(), new String(new char[]{alignment.getData(i, nucleotideSite),
                        alignment.getData(i, nucleotideSite + 1),
                        alignment.getData(i, nucleotideSite + 2)}));
                }
            */

            states.put(alignment.getIdentifier(i).getName(),
                    GeneticCode.getInstance().getCodonIndexFromNucleotides(
                            new char[]{alignment.getData(i, nucleotideSite),
                                    alignment.getData(i, nucleotideSite + 1),
                                    alignment.getData(i, nucleotideSite + 2)}
                    ));
        }
        return states;
    }

    /**
     * Reads a Newick formatted tree
     * @param path path to the tree file
     * @return a PAL tree object
     */
    public static Tree readTree(String path) {
        SimpleTree tree;
        try {
            tree = new ReadTree(path);
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
        return tree;
    }

    /**
     * Gets all the codons that code for a particular set of amino acids
     * @param aminoAcids a collection of amino acid indexes
     * @return a list of codon indexes that code for the given amino acids
     */
    public static List<Integer> getCodonsFromAminoAcids(Collection<Integer> aminoAcids) {
        int[] collected = new int[0];
        for (int a : aminoAcids) {
            if (!GeneticCode.getInstance().isUnknownAminoAcidState(a)) {
                collected = Ints.concat(collected, GeneticCode.getInstance().getCodonIndexFromAminoAcidIndex(a));
            }
        }
        return Ints.asList(collected);
    }

    /**
     * Given a list of codons [observed at a given site], return a list of all observed residues, ordered by frequency.
     * i.e. the most common observed residue first, then the next most common etc. Amino acids that are observed an
     * equal number or times preserve the canonical amino acid order.
     *
     * @param codons a collection of codons (e.g. all codons observed at a particular location in alignment)
     * @return an ordered list of amino acid from most common -> least common, sub-ordered by canonical ordering.
     */
    public static List<Integer> getDistinctAminoAcids(Collection<Integer> codons) {
        // TODO: Christ almighty, there *must* be a better way to do this! Anyway, here we go...

        // A small class to store counts of amino acid residues from the list of codons
        class Residue implements Comparable<Residue>{
            Residue (int index) { this.index = index; }
            int index;
            int count = 0;
            // This sorts residues by the count field (descending) then the amino acid index
            public int compareTo(Residue residue) {
                if (this.count != residue.count)
                    return residue.count - this.count;
                return this.index - residue.index;
            }
        }

        // Create a list to store the 20 amino acids, initial count of 0
        List<Residue> allResidues = Lists.newArrayListWithCapacity(GeneticCode.AMINO_ACID_STATES);
        for (int i = 0; i < GeneticCode.AMINO_ACID_STATES; i++) allResidues.add(new Residue(i));
	
        // Loop through the list of codons
        for (int c : codons) {
            // If we find a valid amino acid, increment its count
            int a = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(c);
            if (!GeneticCode.getInstance().isUnknownAminoAcidState(a)) allResidues.get(a).count++;
        }

        Collections.sort(allResidues);

        // Actually, we only want to return a list of integers of observed residues, so, erm, here it is
        // TODO: Perhaps we should pass the list of Residues back!?
        List<Integer> distinctAAs = Lists.newArrayList();
        for (Residue r : allResidues) if (r.count > 0) distinctAAs.add(r.index);

        return distinctAAs;
    }
}
