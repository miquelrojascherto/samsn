package org.sams.io;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;

import junit.framework.Assert;

import org.junit.Test;
import org.sams.MZData;
import org.sams.MZDataConstants;
import org.sams.SAMSTestCase;

/**
 * TestCase for the PDFWriter class.
 * 
 * @author Miguel Rojas-Cherto
 */
public class PDFWriterTest extends SAMSTestCase {
	
	/**
	 *  Constructor for the PDFWriterTest object.
	 */
	public PDFWriterTest(){
		
	}
	/**
	 * @throws Exception
	 */
	@Test 
	public void testPDFWriter1()throws Exception{
		String pathFile = "src/test/data/cml/F000787.Nom.gr.cml";
		InputStream input = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader = new CMLReader(input);
		MZData mzData = new MZData();
		mzData = mzReader.read(mzData);
		
		String output = "/tmp/output1.pdf";
		PDFWriter pdfWriter = new PDFWriter(output);
		pdfWriter.write(mzData);
		pdfWriter.close();
		
		Assert.assertTrue((new File(output).exists()));
	}

	/**
	 * @throws Exception
	 */
	@Test 
	public void testPDFWriter2()throws Exception{
		String pathFile = "src/test/data/cml/F000504-R.gr.cml";
		InputStream input = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader = new CMLReader(input);
		MZData mzData = new MZData();
		mzData = mzReader.read(mzData);
		
		String output = "/tmp/output2.pdf";
		PDFWriter pdfWriter = new PDFWriter(output);
		pdfWriter.write(mzData);
		pdfWriter.close();
		
		Assert.assertTrue((new File(output).exists()));
	}

	/**
	 * @throws Exception
	 */
	@Test 
	public void testPDFWriter3()throws Exception{
		String pathFile = "src/test/data/cml/F002169_C15H9O4.cml";
		InputStream input = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader = new CMLReader(input);
		MZData mzData = new MZData();
		mzData = mzReader.read(mzData);
		
		String output = "/tmp/output3.pdf";
		PDFWriter pdfWriter = new PDFWriter(output);
		pdfWriter.write(mzData);
		pdfWriter.close();
		
		Assert.assertTrue((new File(output).exists()));
	}
	/**
	 * Plot tree containing colors
	 * 
	 * @throws Exception
	 */
	@Test 
	public void testPDFColors()throws Exception{
		String pathFile = "src/test/data/cml/mzDataWithColors.cml";
		InputStream input = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader = new CMLReader(input);
		MZData mzData = new MZData();
		mzData = mzReader.read(mzData);
		
		String output = "/tmp/output4.pdf";
		PDFWriter pdfWriter = new PDFWriter(output);
		pdfWriter.write(mzData);
		pdfWriter.close();
		
		Assert.assertTrue((new File(output).exists()));
	}

	/**
	 * Plotting tree including inchi
	 * 
	 * @throws Exception
	 */
	@Test 
	public void testPDFInChI_set()throws Exception{
		String pathFile = "src/test/data/cml/mzDataWithColors.cml";
		InputStream input = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader = new CMLReader(input);
		MZData mzData = new MZData();
		mzData = mzReader.read(mzData);
		
		String output = "/tmp/output5.pdf";
		PDFWriter pdfWriter = new PDFWriter(output);
		String inchi = "InChI=1/C7H8N4O3/c1-10-3-4(8-6(10)13)11(2)7(14)9-5(3)12/h1-2H3,(H,8,13)(H,9,12,14)";
		pdfWriter.write(mzData,inchi);
		pdfWriter.close();
		
		Assert.assertTrue((new File(output).exists()));
	}


	/**
	 * Plotting tree including inchi
	 * 
	 * @throws Exception
	 */
	@Test 
	public void testPDFLosses()throws Exception{
		String pathFile = "src/test/data/cml/F002169_C15H9O4.cml";
		InputStream input = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader = new CMLReader(input);
		MZData mzData = new MZData();
		mzData = mzReader.read(mzData);
		
		String output = "/tmp/output6.pdf";
		PDFWriter pdfWriter = new PDFWriter(output);
		pdfWriter.setType(MZDataConstants.LOSSES);
		pdfWriter.write(mzData);
		pdfWriter.close();
		
		Assert.assertTrue((new File(output).exists()));
	}
	/**
	 * Plot tree containing colors
	 * 
	 * @throws Exception
	 */
	@Test 
	public void testPDFLossesColors()throws Exception{
		String pathFile = "src/test/data/cml/mzDataWithColors.cml";
		InputStream input = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader = new CMLReader(input);
		MZData mzData = new MZData();
		mzData = mzReader.read(mzData);
		
		String output = "/tmp/output7.pdf";
		PDFWriter pdfWriter = new PDFWriter(output);
		pdfWriter.setType(MZDataConstants.LOSSES);
		pdfWriter.write(mzData);
		pdfWriter.close();
		
		Assert.assertTrue((new File(output).exists()));
	}
}
