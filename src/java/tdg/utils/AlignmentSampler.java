package tdg.utils;

import com.google.common.base.Charsets;
import com.google.common.base.Function;
import com.google.common.base.Joiner;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Collections2;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.io.Files;
import com.google.gson.Gson;
import org.apache.commons.math.random.MersenneTwister;
import org.apache.commons.math.random.RandomDataImpl;
import pal.alignment.Alignment;
import tdg.analysis.Defaults;

import java.io.File;
import java.io.Writer;
import java.util.Collection;
import java.util.Map;

/**
 * Given a sequence alignment, creates a given number of new alignments with
 * sites sampled, with replacement, from the original sequence alignment.
 *
 * Useful for bootstrapping
 *
 * @author Asif Tamuri
 * @version $Id$
 */
public class AlignmentSampler {

    private void run(String alignmentFile, int numberOfSamples) throws Exception {
        final Alignment alignment = PhyloUtils.readAlignment(alignmentFile);
        final int codons = alignment.getSiteCount() / 3;
        final RandomDataImpl rdi = new RandomDataImpl(new MersenneTwister());

        System.out.printf("Number of sequences: %s\nSites: %s\nCodons: %s\n", alignment.getSequenceCount(), alignment.getSiteCount(), codons);
        Defaults d = new Defaults();

        for (int i = 1; i <= numberOfSamples; i++) {
            Multimap<String, Integer> bootstrap = ArrayListMultimap.create();
            Map<Integer, Map<Integer, Double>> fitnessDefaults = Maps.newHashMap();

            for (int j = 1; j <= codons; j++) {
                int sampledSite = rdi.nextInt(1, codons);

                // add the codon TLA to our bootstrap
                Map<String, Integer> s = PhyloUtils.getCodonsAtSite(alignment, sampledSite);

                for (Map.Entry<String, Integer> e : s.entrySet()) {
                    bootstrap.put(e.getKey(), e.getValue());
                }

                // get good initial fitness parameters, to speed up convergence
                Map<Integer, Double> defaultF = Defaults.getDefaultFitnesses(sampledSite);
                if (defaultF != null) {
                    defaultF.put(-1, (double) sampledSite );
                    fitnessDefaults.put(j, defaultF);
                }
            }

            // check what we've got
            System.out.printf("Keys: %s (%s) %s\n", bootstrap.keySet().size(), bootstrap.keys().entrySet().iterator().next().getElement(),
                    bootstrap.get("Eubalaena_australis").size());

            // check number of codons for each key
            for (String k : bootstrap.keys()) {
                Collection<Integer> c = bootstrap.get(k);
                if (c.size() != codons) {
                    System.out.printf("%s != %s\n", c.size(), codons);
                }
            }



            // save the new alignment, but output in the same order as the original sequence

            Writer w = Files.newWriter(new File(alignmentFile + ".bootstrap." + i), Charsets.US_ASCII);
            w.append(String.format("%s %s\n", alignment.getSequenceCount(), alignment.getSiteCount()));

            for (int j = 0; j < alignment.getSequenceCount(); j++) {

                Collection<String> c = Collections2.transform(
                        bootstrap.get(alignment.getIdentifier(j).getName()),
                        new Function<Integer, String>() {
                            @Override
                            public String apply(Integer input) {
                                String s = GeneticCode.getInstance().getCodonTLA(input);
                                return s.equals("???") ? "---" : s;
                            }
                        });

                w.append(String.format("%s      %s\n", alignment.getIdentifier(j).getName(), Joiner.on("").join(c)));
            }

            w.close();
            System.out.printf("Created bootstrapped sample %s.\n", i);

            // write the new fitness file
            Gson g = new Gson();
            Files.write(g.toJson(fitnessDefaults), new File("fitness.defaults." + i), Charsets.US_ASCII);

        }
    }

    public static void main(String[] args) throws Exception {
        AlignmentSampler as = new AlignmentSampler();
        as.run(args[0], Integer.parseInt(args[1]));
    }
}
