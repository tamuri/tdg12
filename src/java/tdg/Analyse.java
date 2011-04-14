package tdg;

import com.beust.jcommander.JCommander;
import com.google.common.collect.Sets;
import pal.alignment.Alignment;
import pal.tree.Tree;
import tdg.analysis.SiteAnalyser;
import tdg.models.TDGGlobals;
import tdg.utils.PhyloUtils;

import java.util.Set;
import java.util.concurrent.*;

/**
 * @author Asif Tamuri
 * @version $Id$
 */
public class Analyse {
    Options options;

    public Analyse(Options options) {
        this.options = options;
    }

    public static void main(String... args) {
        Options options = new Options();
        JCommander jc = new JCommander(options);

        if (args.length == 0) {
            jc.usage();
            System.out.println("Options preceded by an asterisk are required.");
            System.out.println("Example: -t HA.tree -s HA.co -site 123 -tau 1e-6 -kappa 8.0004 -pi 0.21051,0.19380,0.40010,0.19559 -mu 3.0 -homogonly");
            System.exit(0);
        } else {
            jc.parse(args);
        }

        Analyse analyse = new Analyse(options);
        analyse.run();
    }

    private void run() {
        long startTime = System.currentTimeMillis();

        // Global parameters for the TdG model
        final TDGGlobals tdgGlobals = new TDGGlobals(options.tau, options.kappa, options.pi, options.mu, options.gamma);
        final Tree tree = PhyloUtils.readTree(options.treeFile);
        final Alignment alignment = PhyloUtils.readAlignment(options.alignmentFile);
        final int sites = alignment.getSiteCount() / 3;
        System.out.printf("tdg.Analyse - %s alignment file has %s sequences, each with %s codon sites.\n", options.alignmentFile, alignment.getSequenceCount(), sites);

        // Something to collect results from the analysis of each site
        Set<Future<double[]>> results = Sets.newHashSet();

        int startSite = 1;
        int endSite = sites;

        // Only analyse the given site, not the entire alignment
        if (options.site != 0) {
            startSite = endSite = options.site;
            options.threads = 1;
        }

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

        System.out.printf("tdg.Analyse - Done!\n");
        System.out.printf("tdg.Analyse - Total homogeneous lnL: %s\n", sumHomLnl);
        System.out.printf("tdg.Analyse - Total non-homogeneous lnL: %s\n", sumNonHomLnl);
        System.out.printf("tdg.Analyse - Total time: %s ms (%.2f m).\n", System.currentTimeMillis() - startTime, (System.currentTimeMillis() - startTime) / 60000.0);
        threadPool.shutdown();
    }
 
    private class SiteAnalyserThread implements Callable<double[]> {
        private final int site;
        private final TDGGlobals globals;
        private final Alignment alignment;
        private final Tree tree;
        private final Options options;

        private SiteAnalyserThread(int site, Tree tree, Alignment alignment, TDGGlobals globals, Options options) {
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
            System.out.printf("Site %s - Done. Homogeneous lnL = %s. Non-homogeneous lnL = %s\n", this.site, sa.getHomogeneousLikelihood(), sa.getNonHomogeneousLikelihood());
            return new double[]{sa.getHomogeneousLikelihood(), sa.getNonHomogeneousLikelihood()};
        }
    }
}
