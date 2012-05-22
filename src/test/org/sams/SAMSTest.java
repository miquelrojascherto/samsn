package org.sams;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;
import org.sams.io.CMLReaderTest;
import org.sams.io.CMLWriterTest;
import org.sams.io.MZXMLReaderTest;
import org.sams.io.PDFWriterTest;
import org.sams.manipulator.FingerMZDataManipulator;
import org.sams.manipulator.MZDataManipulatorTest;
import org.sams.spect.ScanMZXMLTest;


/**
 * TestSuite that runs all the JUnit tests for the formula module.
 *
 */
@RunWith(value=Suite.class)
@SuiteClasses(value={
	//
	FingerMZDataTest.class,
	MZDataTest.class,
	//
	CMLReaderTest.class,
    CMLWriterTest.class,
    MZXMLReaderTest.class,
    PDFWriterTest.class,
    //
    FingerMZDataManipulator.class,
    MZDataManipulatorTest.class,
    //
    ScanMZXMLTest.class
})

/**
 * TestCase for all classes.
 * 
 * @author Miguel Rojas-Cherto
 */
public class SAMSTest {}
