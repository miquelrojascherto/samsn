package org.sams;

import org.openscience.cdk.tools.LoggingTool;

/**
 * Super class for <b>all</b> SAMS TestCase implementations that ensures that
 * the LoggingTool is configured. This is the JUnit4 version of SAMSTestCase.
 *
 * @author Miguel Rojas-Cherto
 */
public class SAMSTestCase {
	static {
        LoggingTool.configureLog4j();
    }
}
