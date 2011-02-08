package tdg;

import com.beust.jcommander.Parameter;

/**
 * User: atamuri
 * Date: 07/02/11
 * Time: 17:47
 */
public class ClusterOptions extends Options {
    @Parameter(names = "-hostConnections", description = "Maximum number of connections from master to *each* slave. Default = 5", required = false)
    public int hostConnections = 5;

    @Parameter(names = "-hostTimeout", description = "Number of seconds to timeout a connection to a slave. Default = 6000000", required = false)
    public int hostTimeout = 6000000;

    @Parameter(names = "-hostThreadPool", description = "Number of threads (inc. idle) held by the master server. Default = 30", required = false)
    public int hostThreadPool = 30;
}
