package tdg.cli;

import com.beust.jcommander.Parameter;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk
 */
public class ApproxOption {
    @Parameter(names = "-useapprox", description = "Use the approximate method to optimise the likelihood")
    public boolean useapprox = false;
}
