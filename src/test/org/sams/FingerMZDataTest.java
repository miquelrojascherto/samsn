package org.sams;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import junit.framework.Assert;
import junit.framework.TestCase;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.sams.io.CMLReader;
import org.sams.manipulator.MZDataManipulator;

public class FingerMZDataTest  extends TestCase {
	static DefaultChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
	String pathFile = "src/test/data/cml/spect_react.cml";
	
	public void testGetFingerprint0()throws Exception{
		List<String> value1 = Arrays.asList("C11H10NO2","C11H10NO2||C10H10N","C11H10NO2||C10H10NO","C11H10NO2||C10H8N","C11H10NO2||C11H8NO","C11H10NO2||C9H8NO","C11H10NO2||C9H8NO||C8H8N");
		
		
		InputStream input1 = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader1 = new CMLReader(input1);
		MZData mzData1 = new MZData();
		mzData1 = MZDataManipulator.group(mzReader1.read(mzData1),0.4);
	    
		FingerMZData fingerprinter = new FingerMZData();
		List<String> bs = fingerprinter.getFingerprintL0(mzData1);
		
		Assert.assertEquals(value1.size(), bs.size());
		for(String value : value1){
			Assert.assertTrue(bs.contains(value));
		}
	}
	public void testGetFingerprint0_loss()throws Exception{
		List<String> value1 = Arrays.asList("C2H2O","C2H2O||CO","CO2","CO","CH2O2","H2O");
		
		InputStream input1 = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader1 = new CMLReader(input1);
		MZData mzData1 = new MZData();
		mzData1 = MZDataManipulator.group(mzReader1.read(mzData1),0.4);
	    
		FingerMZData fingerprinter = new FingerMZData();
		List<String> bs = fingerprinter.getFingerprintL0_loss(mzData1);
		
		Assert.assertEquals(value1.size(), bs.size());
		for(String value : value1){
			Assert.assertTrue(bs.contains(value));
		}
	}
	public void testGetFingerprintL1()throws Exception{
		List<String> value1 = Arrays.asList("C11H10NO2","C10H10N","C10H10NO","C10H8N","C11H8NO","C9H8NO","C8H8N");
		
		InputStream input1 = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader1 = new CMLReader(input1);
		MZData mzData1 = new MZData();
		mzData1 = MZDataManipulator.group(mzReader1.read(mzData1),0.4);
	    
		FingerMZData fingerprinter = new FingerMZData();
		List<String> bs = fingerprinter.getFingerprintL1(mzData1);
		
		Assert.assertEquals(value1.size(), bs.size());
		for(String value : value1){
			Assert.assertTrue(bs.contains(value));
		}
	}
	public void testGetFingerprintL1_Loss()throws Exception{
		List<String> value1 = Arrays.asList("C2H2O","CO","CO2","CH2O2","H2O");
		
		InputStream input1 = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader1 = new CMLReader(input1);
		MZData mzData1 = new MZData();
		mzData1 = MZDataManipulator.group(mzReader1.read(mzData1),0.4);
	    
		FingerMZData fingerprinter = new FingerMZData();
		List<String> bs = fingerprinter.getFingerprintL1_Loss(mzData1);
		
		Assert.assertEquals(value1.size(), bs.size());
		for(String value : value1){
			Assert.assertTrue(bs.contains(value));
		}
	}
	public void testGetFingerprintL2()throws Exception{
		List<String> value1 = Arrays.asList("C11H10NO2@C10H10N","C11H10NO2@C10H10NO","C11H10NO2@C10H8N","C11H10NO2@C11H8NO","C11H10NO2@C9H8NO","C9H8NO@C8H8N");
		
		InputStream input1 = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader1 = new CMLReader(input1);
		MZData mzData1 = new MZData();
		mzData1 = MZDataManipulator.group(mzReader1.read(mzData1),0.4);
	    
		FingerMZData fingerprinter = new FingerMZData();
		List<String> bs = fingerprinter.getFingerprintL2(mzData1);
		
		Assert.assertEquals(value1.size(), bs.size());
		for(String value : value1){
			Assert.assertTrue(bs.contains(value));
		}
	}

	public void testGetFingerprintL2_Loss()throws Exception{
		List<String> value1 = Arrays.asList("C2H2O@CO");
		
		InputStream input1 = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader1 = new CMLReader(input1);
		MZData mzData1 = new MZData();
		mzData1 = MZDataManipulator.group(mzReader1.read(mzData1),0.4);
	    
		FingerMZData fingerprinter = new FingerMZData();
		List<String> bs = fingerprinter.getFingerprintL2_Loss(mzData1);
		
		Assert.assertEquals(value1.size(), bs.size());
		for(String value : value1){
			Assert.assertTrue(bs.contains(value));
		}
	}
	public void testGetFingerprintL3()throws Exception{
		List<String> value1 = Arrays.asList("C11H10NO2@C9H8NO@C8H8N");
		
		InputStream input1 = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader1 = new CMLReader(input1);
		MZData mzData1 = new MZData();
		mzData1 = MZDataManipulator.group(mzReader1.read(mzData1),0.4);
	    
		FingerMZData fingerprinter = new FingerMZData();
		List<String> bs = fingerprinter.getFingerprintL3(mzData1);
		
		Assert.assertEquals(value1.size(), bs.size());
		for(String value : value1){
			Assert.assertTrue(bs.contains(value));
		}
	}
	public void testGetFingerprintV2()throws Exception{

		List<String> value1 = Arrays.asList("C11H10NO2@[C10H10N||C10H10NO]","C11H10NO2@[C10H10N||C10H8N]","C11H10NO2@[C10H10N||C11H8NO]",
				"C11H10NO2@[C10H10N||C9H8NO]","C11H10NO2@[C10H10NO||C10H8N]","C11H10NO2@[C10H10NO||C11H8NO]","C11H10NO2@[C10H10NO||C9H8NO]",
				"C11H10NO2@[C10H8N||C11H8NO]","C11H10NO2@[C10H8N||C9H8NO]","C11H10NO2@[C11H8NO||C9H8NO]");
		
		InputStream input1 = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader1 = new CMLReader(input1);
		MZData mzData1 = new MZData();
		mzData1 = MZDataManipulator.group(mzReader1.read(mzData1),0.4);
	    
		FingerMZData fingerprinter = new FingerMZData();
		List<String> bs = fingerprinter.getFingerprintV2(mzData1);
		
		Assert.assertEquals(value1.size(), bs.size());
		for(String value : value1){
			Assert.assertTrue(bs.contains(value));
		}
		
//		List<String> bs2a = fingerprinter.getFingerprintV0(mzData1);
//		List<String> bs2 = new ArrayList<String>(); 
//		for(String bsa : bs2a){
//			List<String> bs_2 = FingerMZData.getPack(bsa,2);
//		}
//		Assert.assertEquals(value1.size(), bs2.size());
//		for(String value : value1){
//			Assert.assertTrue(bs2.contains(value));
//		}
	}

	public void testGetFingerprintV2_Loss()throws Exception{

		List<String> value1 = new ArrayList<String>();
		
		InputStream input1 = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader1 = new CMLReader(input1);
		MZData mzData1 = new MZData();
		mzData1 = MZDataManipulator.group(mzReader1.read(mzData1),0.4);
	    
		FingerMZData fingerprinter = new FingerMZData();
		List<String> bs = fingerprinter.getFingerprintV2_Loss(mzData1);
		
		Assert.assertEquals(value1.size(), bs.size());
	}
	
	public void testGetFingerprintV3()throws Exception{
		List<String> value1 = Arrays.asList("C11H10NO2@[C10H10N||C10H10NO||C10H8N]","C11H10NO2@[C10H10N||C10H10NO||C11H8NO]",
				"C11H10NO2@[C10H10N||C10H10NO||C9H8NO]","C11H10NO2@[C10H10NO||C10H8N||C11H8NO]","C11H10NO2@[C10H10NO||C10H8N||C9H8NO]",
				"C11H10NO2@[C10H8N||C11H8NO||C9H8NO]");
		
		InputStream input1 = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader1 = new CMLReader(input1);
		MZData mzData1 = new MZData();
		mzData1 = MZDataManipulator.group(mzReader1.read(mzData1),0.4);
	    
		FingerMZData fingerprinter = new FingerMZData();
		List<String> bs = fingerprinter.getFingerprintV3(mzData1);
				
		Assert.assertEquals(value1.size(), bs.size());
		for(String value : value1){
			Assert.assertTrue(bs.contains(value));
		}
	}
	public void testGetFingerprintV0()throws Exception{
		List<String> value1 = Arrays.asList("C11H10NO2@[C10H10N||C10H10NO||C10H8N||C11H8NO||C9H8NO]");
		
		InputStream input1 = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader1 = new CMLReader(input1);
		MZData mzData1 = new MZData();
		mzData1 = MZDataManipulator.group(mzReader1.read(mzData1),0.4);
	    
		FingerMZData fingerprinter = new FingerMZData();
		List<String> bs = fingerprinter.getFingerprintV0(mzData1);
		
		Assert.assertEquals(value1.size(), bs.size());
		for(String value : value1){
			Assert.assertTrue(bs.contains(value));
		}
	}

	public void testGetPackList3()throws Exception{
		String pack = "C11H10NO2@[C10H10N||C10H10NO||C10H8N||C11H8NO||C9H8NO]";
		List<String> value1 = Arrays.asList("C11H10NO2@[C10H10N||C10H10NO||C10H8N]","C11H10NO2@[C10H10N||C10H10NO||C11H8NO]",
				"C11H10NO2@[C10H10N||C10H10NO||C9H8NO]","C11H10NO2@[C10H10NO||C10H8N||C11H8NO]","C11H10NO2@[C10H10NO||C10H8N||C9H8NO]",
				"C11H10NO2@[C10H8N||C11H8NO||C9H8NO]");
		
		FingerMZData fingerprinter = new FingerMZData();
		List<String> bs = fingerprinter.getPackList(pack, 3);
		
		Assert.assertEquals(value1.size(), bs.size());
		for(String value : value1){
			Assert.assertTrue(bs.contains(value));
		}
	}
}
