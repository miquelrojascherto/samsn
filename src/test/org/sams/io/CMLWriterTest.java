package org.sams.io;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.StringReader;
import java.io.StringWriter;

import junit.framework.Assert;

import org.junit.Test;
import org.openscience.cdk.ReactionScheme;
import org.sams.MZData;
import org.sams.SAMSTestCase;
import org.xmlcml.cml.element.CMLSpectrumList;

/**
 * TestCase for the CMLWriter class.
 * 
 * @author Miguel Rojas-Cherto
 */
public class CMLWriterTest extends SAMSTestCase {
	
	/**
	 *  Constructor for the CMLWriterTest object.
	 */
	public CMLWriterTest(){
		
	}
	
	/**
	 * @throws Exception
	 */
	@Test 
	public void testCMLWriter_null()throws Exception{
        MZData mzData = new MZData();
		
		StringWriter output = new StringWriter();
        CMLWriter cmlWriter = new CMLWriter(output);
		cmlWriter.write(mzData);
		cmlWriter.close();
		String cmlcode = output.toString();
		
		StringWriter output2 = new StringWriter();
		output2.write("<cml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"/>\n");
		
		BufferedReader reader = new BufferedReader(new StringReader(output2.toString()));
		BufferedReader reader2 = new BufferedReader(new StringReader(cmlcode));
		String str, str2;
		while ((str = reader.readLine()) != null) {
			str2 = reader2.readLine();
			Assert.assertEquals(str2,str);
		}
	}
	/**
	 * @throws Exception
	 */
	@Test 
	public void testCMLWriter()throws Exception{
        MZData mzData = getExampleMZData();
		
		StringWriter output = new StringWriter();
        CMLWriter cmlWriter = new CMLWriter(output);
		cmlWriter.write(mzData);
		cmlWriter.close();
		String cmlcode = output.toString();

		StringWriter output2 = new StringWriter();
		output2.write("<cml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n");
		output2.write("  <spectrumList xmlns=\"http://www.xml-cml.org/schema\"/>\n");
		output2.write("  <moleculeList convention=\"cdk:moleculeSet\" xmlns=\"http://www.xml-cml.org/schema\"/>\n");
		output2.write("  <reactionScheme xmlns=\"http://www.xml-cml.org/schema\"/>\n");
		output2.write("</cml>");
		
		BufferedReader reader = new BufferedReader(new StringReader(output2.toString()));
		BufferedReader reader2 = new BufferedReader(new StringReader(cmlcode));
		String str, str2;
		while ((str = reader.readLine()) != null) {
			str2 = reader2.readLine();
			Assert.assertEquals(str2,str);
		}
	}

	private MZData getExampleMZData() {
		MZData mzData = new MZData();
		CMLSpectrumList specList = new CMLSpectrumList();
		ReactionScheme rS0 = new ReactionScheme();
		mzData.setListSpectra(specList);
		mzData.setListReactions(rS0);
		return mzData;
	}

	
	/**
	 * CML containing spectra with peaks and atomcontainer information
	 * 
	 * @throws Exception
	 */
	@Test
	public void testCMLWriter2()throws Exception{

        String pathFile = "src/test/data/cml/WithInChIs.cml";
		
        // Original        
        FileInputStream fstream = new FileInputStream(pathFile);
        DataInputStream in = new DataInputStream(fstream);
        BufferedReader br = new BufferedReader(new InputStreamReader(in));
		
        // Process
        InputStream input = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader = new CMLReader(input);
		MZData mzData = new MZData();
		mzData = mzReader.read(mzData);
		
		StringWriter output2 = new StringWriter();
		CMLWriter cmlWriter2 = new CMLWriter(output2);
		cmlWriter2.write(mzData);
		cmlWriter2.close();
		
		BufferedReader br2 = new BufferedReader(new StringReader(output2.toString()));
		
		String line;
	    while ((line = br.readLine()) != null) {
	    	String line2 = br2.readLine();
	    	Assert.assertEquals(line, line2);
	    }

	    Assert.assertNull("Actual had more lines then the expected.", br2.readLine());
	    Assert.assertNull("Expected had more lines then the actual.", br.readLine());
	}
}