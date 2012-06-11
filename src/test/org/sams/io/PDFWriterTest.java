package org.sams.io;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.StringWriter;

import junit.framework.Assert;

import org.junit.Test;
import org.sams.MZData;
import org.sams.MZDataConstants;
import org.sams.SAMSTestCase;
import org.sams.manipulator.MZDataManipulator;

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

	/**
	 * Plot tree containing colors
	 * 
	 * @throws Exception
	 */
	@Test 
	public void testPDFOnlyLossesColors()throws Exception{
		String pathFile1 = "/tmp/170.04.cml";
		String pathFile2 = "/tmp/581570ac84464d7bb32f57d801776a89-1.cml";
		InputStream input1 = new BufferedInputStream(new FileInputStream(pathFile1));
		CMLReader mzReader1 = new CMLReader(input1);
		MZData mzData1 = new MZData();
		mzData1 = mzReader1.read(mzData1);
		InputStream input2 = new BufferedInputStream(new FileInputStream(pathFile2));
		CMLReader mzReader2 = new CMLReader(input2);
		MZData mzData2 = new MZData();
		mzData2 = mzReader2.read(mzData2);
		
		MZData mzDataCol21 = MZDataManipulator.getColored(mzData2,mzData1);
		MZData mzDataCol12 = MZDataManipulator.getColored(mzData1,mzData2);

		String output1 = "/tmp/output8.1.pdf";
		PDFWriter pdfWriter1L = new PDFWriter(output1);
		pdfWriter1L.setType(MZDataConstants.LOSSES);
		pdfWriter1L.write(mzDataCol21);
		pdfWriter1L.close();
		
		String output2 = "/tmp/output8.2.pdf";
		PDFWriter pdfWriter2L = new PDFWriter(output2);
		pdfWriter2L.setType(MZDataConstants.LOSSES);
		pdfWriter2L.write(mzDataCol12);
		pdfWriter2L.close();
		
//		StringWriter output = new StringWriter();
//        CMLWriter cmlWriter = new CMLWriter(output);
//		cmlWriter.write(mzDataCol21);
//		cmlWriter.close();
//		String cmlcode = output.toString();
//		System.out.println(cmlcode);
		
	}
}
