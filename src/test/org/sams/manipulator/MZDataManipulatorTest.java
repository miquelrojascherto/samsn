package org.sams.manipulator;

import java.io.BufferedInputStream;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.UUID;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.interfaces.IMoleculeSet;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.openscience.cdk.tools.manipulator.ReactionSchemeManipulator;
import org.sams.MZData;
import org.sams.MZDataConstants;
import org.sams.SAMSTestCase;
import org.sams.io.CMLReader;
import org.sams.io.CMLWriter;
import org.xmlcml.cml.element.CMLPeak;

/**
 * TestCase for the MZDataManipulator class.
 * 
 * @author Miguel Rojas-Cherto
 */
public class MZDataManipulatorTest extends SAMSTestCase {
	static DefaultChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
	
	/**
	 *  Constructor for the MZDataManipulatorTest object.
	 */
	public MZDataManipulatorTest(){
		
	}

	/**
	 * Using the table with minimal information
	 * 
	 * @throws Exception
	 * @sams-bug SAMS-27
	 */
	@Test 
	public void testGetMZDataMatrixFromFileShort()throws Exception{
		
		double[][] matrix = MZDataManipulator.getMZDataMatrixFromFile("src/test/data/txt/matrixPeaksShort.txt");
		MZData mzData = new MZData();
		
		mzData = MZDataManipulator.getMZData(matrix, mzData);
		
		Assert.assertEquals(5, mzData.getListSpectra().getSpectrumElements().size());
		
	}
	
	/**
	 * CML containing spectra with peaks
	 * 
	 * @throws Exception
	 */
	@Test 
	public void testReadProperties()throws Exception{
		double[][] matrix = MZDataManipulator.getMZDataMatrixFromFile("src/test/data/txt/dilu_a.txt");
		MZData mzData = new MZData();
		
		mzData = MZDataManipulator.getMZData2(matrix,mzData,"src/test/data/mzXML/dilu_a.mzXML","","",1);

		Map<Object, Object> properties = mzData.getProperties();
		Assert.assertEquals(5, properties.get(MZDataConstants.NUM_GROUPS));
		Assert.assertEquals(143, properties.get(MZDataConstants.SCANS));
		
	}
	
	@Test 
	public void testGroup_MZData_MZData_int()throws Exception{
		String pathFileA = "src/test/data/cml/groupC15H9O4_A.cml";
		String pathFileB = "src/test/data/cml/groupC15H9O4_B.cml";
		InputStream inputA = new BufferedInputStream(new FileInputStream(pathFileA));
		CMLReader mzReaderA = new CMLReader(inputA);
		MZData mzData1 = new MZData();
		mzData1 = mzReaderA.read(mzData1);
		mzData1.setID(UUID.randomUUID().toString());
		
		InputStream inputB = new BufferedInputStream(new FileInputStream(pathFileB));
		CMLReader mzReaderB = new CMLReader(inputB);
		MZData mzData2 = new MZData();
		mzData2 = mzReaderB.read(mzData2);
		mzData2.setID(UUID.randomUUID().toString());
		
		List<MZData> mzDataList = new ArrayList<MZData>();
		mzDataList.add(mzData1);
		mzDataList.add(mzData2);

		MZData mzData3 = MZDataManipulator.group(mzDataList,0.5);
		
		CMLPeak peak1 = mzData1.getListSpectra().getSpectrumElements().get(0).getPeakListElements().get(0).getPeakElements().get(0);
		CMLPeak peak2 = mzData2.getListSpectra().getSpectrumElements().get(0).getPeakListElements().get(0).getPeakElements().get(0);
		CMLPeak peak3 = mzData3.getListSpectra().getSpectrumElements().get(0).getPeakListElements().get(0).getPeakElements().get(0);
		
		Assert.assertEquals(223.07635498046875,peak1.getXValue(),0.001);
		Assert.assertEquals(4.512167E7,peak1.getYValue(),0.1);

		Assert.assertEquals(223.08046875,peak2.getXValue(),0.001);
		Assert.assertEquals(4.2167E7,peak2.getYValue(),0.1);

		Assert.assertEquals(223.0784118652344,peak3.getXValue(),0.001);
		Assert.assertEquals(4.3644335E7,peak3.getYValue(),0.1);
		
	}

	@Test 
	public void testGetTheoretical_MZData()throws Exception{
		String pathFileA = "src/test/data/cml/groupC15H9O4_A.cml";
		InputStream inputA = new BufferedInputStream(new FileInputStream(pathFileA));
		CMLReader mzReaderA = new CMLReader(inputA);
		MZData mzData1 = new MZData();
		mzData1 = mzReaderA.read(mzData1);
		
		CMLPeak peak1 = mzData1.getListSpectra().getSpectrumElements().get(0).getPeakListElements().get(0).getPeakElements().get(0);
		double mass1 = peak1.getXValue(); //223.07635498046875
		
		MZData mzData2 = MZDataManipulator.getTheoretical(mzData1);
		
		// Peak 1
		String form = "C15H11O2"; //223.07535601209074
		IMolecularFormula formula = MolecularFormulaManipulator.getMajorIsotopeMolecularFormula(form,builder);
	    formula.setCharge(-1);
	    double mass = MolecularFormulaManipulator.getTotalExactMass(formula);
	    
	    
		CMLPeak peak2 = mzData2.getListSpectra().getSpectrumElements().get(0).getPeakListElements().get(0).getPeakElements().get(0);
		double mass2 = peak2.getXValue();
	    Assert.assertEquals(mass,mass2,0.0000001);
	    Assert.assertEquals(223.07635498046875,mass1,0.000001);
		Assert.assertNotSame(peak1.getXValue(),mass2);


		// Peak 2
		CMLPeak peak3 = mzData2.getListSpectra().getSpectrumElements().get(1).getPeakListElements().get(0).getPeakElements().get(0);
		String form3 = "C14H11O1"; //195.08153855190926
		IMolecularFormula formula3 = MolecularFormulaManipulator.getMajorIsotopeMolecularFormula(form3,builder);
	    formula3.setCharge(-1);
	    double mass3 = MolecularFormulaManipulator.getTotalExactMass(formula3);
	    Assert.assertEquals(peak3.getXValue(),mass3,0.000001);

	    
	}
	@Test 
	public void testGetMZData_colored()throws Exception{
		String pathFileA = "src/test/data/cml/groupC15H9O4_A.cml";
		InputStream inputA = new BufferedInputStream(new FileInputStream(pathFileA));
		CMLReader mzReaderA = new CMLReader(inputA);
		MZData mzDataA = new MZData();
		mzDataA = mzReaderA.read(mzDataA);
		
		String pathFileB = "src/test/data/cml/groupC15H9O4_A.cml";
		InputStream inputB = new BufferedInputStream(new FileInputStream(pathFileB));
		CMLReader mzReaderB = new CMLReader(inputB);
		MZData mzDataB = new MZData();
		mzDataB = mzReaderB.read(mzDataB);
		
		MZData mzDataColA = MZDataManipulator.getColored(mzDataA,mzDataB);
		
		HashMap<String, String> ha = new HashMap<String, String>();
		ha.put("21","2");
		ha.put("24","3");
		ha.put("5","2");
		ha.put("6","3");
		ha.put("1","1");
		ha.put("11","2");
		ha.put("4","2");
		
		IMoleculeSet moleculesSet = ReactionSchemeManipulator.getAllMolecules(mzDataColA.getListReactions());
		for(IAtomContainer mol:moleculesSet.molecules()){
			if(ha.containsKey(mol.getID())){
				Assert.assertSame(ha.get(mol.getID()), mol.getProperty(MZDataConstants.COLOR));
			}
		}
		
	}
	@Test 
	public void testGetMZData_colored2()throws Exception{
		String pathFileA = "src/test/data/cml/HMDB01982_pos02.gr.cml";
		InputStream inputA = new BufferedInputStream(new FileInputStream(pathFileA));
		CMLReader mzReaderA = new CMLReader(inputA);
		MZData mzDataA = new MZData();
		mzDataA = mzReaderA.read(mzDataA);
		
		String pathFileB = "src/test/data/cml/HMDB01982_pos02-R5.gr.cml";
		InputStream inputB = new BufferedInputStream(new FileInputStream(pathFileB));
		CMLReader mzReaderB = new CMLReader(inputB);
		MZData mzDataB = new MZData();
		mzDataB = mzReaderB.read(mzDataB);
		
		MZData mzDataColA = MZDataManipulator.getColored(mzDataA,mzDataB);
			
		HashMap<String, String> ha = new HashMap<String, String>();
		ha.put("4","2");
		ha.put("17","3");
		ha.put("15","3");
		ha.put("5","0");
		ha.put("6","0");
		
		IMoleculeSet moleculesSet = ReactionSchemeManipulator.getAllMolecules(mzDataColA.getListReactions());
		for(IAtomContainer mol:moleculesSet.molecules()){
			if(ha.containsKey(mol.getID())){
				Assert.assertSame(ha.get(mol.getID()), mol.getProperty(MZDataConstants.COLOR));
			}
		}
	}

	@Test 
	public void testSetINCHIs()throws Exception{
		String pathFileA = "src/test/data/cml/HMDB01982_pos02.gr.cml";
		InputStream inputA = new BufferedInputStream(new FileInputStream(pathFileA));
		CMLReader mzReaderA = new CMLReader(inputA);
		MZData mzDataA = new MZData();
		mzDataA = mzReaderA.read(mzDataA);

		HashMap<String, String> ha = new HashMap<String, String>();
		ha.put("4","InChI=1");
		ha.put("17","InChI=2");
		ha.put("15","InChI=3");
		ha.put("5","InChI=4");
		ha.put("6","InChI=5");
		
		MZData mzDataInch = MZDataManipulator.SetInChIs(mzDataA,ha);
		
		IMoleculeSet moleculesSet = ReactionSchemeManipulator.getAllMolecules(mzDataInch.getListReactions());
		for(IAtomContainer mol:moleculesSet.molecules()){
			if(ha.containsKey(mol.getID())){
				Assert.assertSame(ha.get(mol.getID()), mol.getProperty(CDKConstants.INCHI));
			}else{
				Assert.assertNull(mol.getProperty(CDKConstants.INCHI));
			}
		}
	}

	@Test
	public void testGetLevels()throws Exception{
		List<Integer> lev = Arrays.asList(0,1,1,2);
		double[][] matrix = new double[4][5];
		//0:ID, 1:pred, 2:m/z, 3:int, 4:charge
		matrix[0][0] = 1.0;matrix[0][1] = 0.0;matrix[0][2] = 134.98;matrix[0][3] = 123442334.98;matrix[0][4] = 1.0;
		matrix[1][0] = 2.0;matrix[1][1] = 1.0;matrix[1][2] = 97.2;matrix[1][3] = 5534534.98;matrix[1][4] = 1.0;
		matrix[2][0] = 3.0;matrix[2][1] = 1.0;matrix[2][2] = 45.45;matrix[2][3] = 989898.555;matrix[2][4] = 1.0;
		matrix[3][0] = 4.0;matrix[3][1] = 3.0;matrix[3][2] = 25.45;matrix[3][3] = 9898.555;matrix[3][4] = 1.0;
		
		List<Integer> levelsL = MZDataManipulator.extractLevels(matrix);
		for(int i = 0 ; i < levelsL.size(); i++){
			Assert.assertEquals(lev.get(i),levelsL.get(i));
		}
	}
	
	@Test
	public void testSplit()throws Exception{
		
		String pathFileA = "src/test/data/cml/F002169_C15H9O4.cml";
		InputStream inputA = new BufferedInputStream(new FileInputStream(pathFileA));
		CMLReader mzReaderA = new CMLReader(inputA);
		MZData mzData = new MZData();
		mzData = mzReaderA.read(mzData);

		List<MZData> mzDataList = MZDataManipulator.splitMZData(mzData);
		
		int numGroup = 13;
		Assert.assertEquals(numGroup, mzDataList.size());

		MZData mzData1 = mzDataList.get(0);
		Map<Object, Object> properties = mzData1.getProperties();
		int numGroupB = Integer.parseInt((String)properties.get(MZDataConstants.NUM_GROUPS));
		Assert.assertEquals(1, numGroupB);
	
		Assert.assertEquals(6, mzData1.getListSpectra().getSpectrumElements().size());

		MZData mzData2 = mzDataList.get(1);
		Map<Object, Object> properties2 = mzData2.getProperties();
		int numGroupC = Integer.parseInt((String)properties2.get(MZDataConstants.NUM_GROUPS));
		Assert.assertEquals(1, numGroupC);
	
		Assert.assertEquals(7, mzData2.getListSpectra().getSpectrumElements().size());

	}
	
	@Test
	public void testGetRt()throws Exception{
		
		String pathFileA = "src/test/data/cml/F002169_C15H9O4.cml";
		InputStream inputA = new BufferedInputStream(new FileInputStream(pathFileA));
		CMLReader mzReaderA = new CMLReader(inputA);
		MZData mzData = new MZData();
		mzData = mzReaderA.read(mzData);

		List<MZData> mzDataList = MZDataManipulator.splitMZData(mzData);
		
		MZData mzData1 = mzDataList.get(0);
		Assert.assertEquals(58.45, MZDataManipulator.getRt(mzData1),0.01);
		
		MZData mzData2 = mzDataList.get(1);
		Assert.assertEquals(58.38, MZDataManipulator.getRt(mzData2),0.01);

	}

}
