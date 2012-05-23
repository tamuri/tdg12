package tdg;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import com.google.common.collect.Sets;
import pal.alignment.Alignment;
import pal.tree.Tree;
import tdg.cli.AnalyseOptions;
import tdg.model.TDGGlobals;
import tdg.utils.PhyloUtils;

import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * The main class for multithreaded (but not distributed) analysis. Most people should use this to run the TdG12 program.
 * Loads the tree, alignment and sets up the global parameters. Then calls SiteAnalyser for the site(s) you want to
 * analyse.
 *
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 * @see SiteAnalyser
 */
public class Analyse {
    AnalyseOptions options;

    public Analyse(AnalyseOptions options) {
        this.options = options;
    }

    public static void main(String... args) {
        AnalyseOptions options = new AnalyseOptions();
        JCommander jc = new JCommander(options);
        jc.setProgramName("java -cp " + Constants.PROGRAM_JAR + " tdg.Analyse ");

        try {
            jc.parse(args);
        } catch (ParameterException pe) {
            System.out.printf("Error: %s\n\n", pe.getMessage());
            jc.usage();
            System.out.println("Options preceded by an asterisk are required.");
            System.out.println("Example: -t HA.tree -s HA.co -site 123 -tau 1e-6 -kappa 8.0004 -pi 0.21051,0.19380,0.40010,0.19559 -mu 3.0");
            System.exit(0);
        }

        Analyse analyse = new Analyse(options);
        analyse.run();
    }

    private void run() {
        long startTime = System.currentTimeMillis();

        // Global parameters for the TdG model
        final TDGGlobals tdgGlobals = new TDGGlobals(options.globals.tau, options.globals.kappa, options.globals.pi, options.globals.mu, options.globals.gamma);
        final Tree tree = PhyloUtils.readTree(options.treeFile);
        final Alignment alignment = PhyloUtils.readAlignment(options.alignmentFile);
        final int sites = alignment.getSiteCount() / 3; // Alignment object is a nucleotide alignment
        System.out.printf("tdg.Analyse - %s alignment file has %s sequences, each with %s codon sites.\n", options.alignmentFile, alignment.getSequenceCount(), sites);

        // Something to collect results from the analysis of each site
        Set<Future<double[]>> results = Sets.newHashSet();

        // By default, we analyse the entire alignment
        int startSite = 1;
        int endSite = sites;

        // Analyse some specific site(s) in the alignment
        if (options.site != null) {
            if (options.site.first.equals(options.site.second)) {
                startSite = endSite = options.site.first;
            } else {
                startSite = options.site.first;
                endSite = options.site.second;
            }
        }

        System.out.printf("tdg.Analyse - Analysing location(s) %s to %s\n", startSite, endSite);

        final int threads = options.threads; // Runtime.getRuntime().availableProcessors();
        final ExecutorService threadPool = Executors.newFixedThreadPool(threads);
        System.out.printf("tdg.Analyse - Running with %s thread(s).\n", threads);

        // Add each site analysis to the thread pool
        for (int site = startSite; site <= endSite; site++) {
            Future<double[]> future = threadPool.submit(new SiteAnalyserThread(site, tree, alignment, tdgGlobals, options));
            results.add(future);
        }

        // Collect the results of the analysis
        double sumHomLnl = 0, sumNonHomLnl = 0;
        for (Future<double[]> future : results) {
            try {
                sumHomLnl += future.get()[0];
                sumNonHomLnl += future.get()[1];
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        long endTime = System.currentTimeMillis();

        System.out.printf("tdg.Analyse - Done!\n");
        System.out.printf("tdg.Analyse - Total homogeneous lnL: %s\n", sumHomLnl);
        System.out.printf("tdg.Analyse - Total non-homogeneous lnL: %s\n", sumNonHomLnl);
        System.out.printf("tdg.Analyse - Total time: %s ms (%.2f m).\n", endTime - startTime, (endTime - startTime) / 60000.0);
        threadPool.shutdown();
    }

    private class SiteAnalyserThread implements Callable<double[]> {
        private final int site;
        private final TDGGlobals globals;
        private final Alignment alignment;
        private final Tree tree;
        private final AnalyseOptions options;

        SiteAnalyserThread(int site, Tree tree, Alignment alignment, TDGGlobals globals, AnalyseOptions options) {
            this.site = site;
            this.tree = tree;
            this.alignment = alignment;
            this.globals = globals;
            this.options = options;
        }

        @Override
        public double[] call() {
            SiteAnalyser sa = new SiteAnalyser(tree, alignment, globals, site, options);
            sa.run();
            System.out.printf("Site %s - Done. Homogeneous lnL = %s. Non-homogeneous lnL = %s\n", this.site, sa.getHomogeneousLikelihood(), sa.getHeterogeneousLikelihood());
            return new double[]{sa.getHomogeneousLikelihood(), sa.getHeterogeneousLikelihood()};
        }
    }
}
