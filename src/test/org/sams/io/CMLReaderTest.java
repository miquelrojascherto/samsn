package org.sams.io;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.StringWriter;
import java.util.List;
import java.util.Map;

import junit.framework.Assert;

import org.junit.Test;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.interfaces.IMolecularFormulaSet;
import org.openscience.cdk.interfaces.IMoleculeSet;
import org.openscience.cdk.interfaces.IReactionScheme;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.openscience.cdk.tools.manipulator.ReactionSchemeManipulator;
import org.sams.MZData;
import org.sams.MZDataConstants;
import org.sams.SAMSTestCase;
import org.sams.manipulator.MZDataManipulator;
import org.xmlcml.cml.base.CMLElements;
import org.xmlcml.cml.element.CMLConditionList;
import org.xmlcml.cml.element.CMLMetadata;
import org.xmlcml.cml.element.CMLMetadataList;
import org.xmlcml.cml.element.CMLPeak;
import org.xmlcml.cml.element.CMLScalar;
import org.xmlcml.cml.element.CMLSpectrum;

/**
 * TestCase for the CMLreader class.
 * 
 * @author Miguel Rojas-Cherto
 */
public class CMLReaderTest extends SAMSTestCase {
	
	/**
	 *  Constructor for the CMLreaderTest object.
	 */
	public CMLReaderTest(){
		
	}
	/**
	 * CML containing only spectra
	 * 
	 * @throws Exception
	 */
	@Test 
	public void testCMLReaderSpectra()throws Exception{
        String pathFile = "src/test/data/cml/multSpectra.cml";
		InputStream input = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader = new CMLReader(input);
		MZData mzData = new MZData();
		mzData = mzReader.read(mzData);

		Assert.assertNotNull(mzData.getListSpectra());
		Assert.assertNotNull(mzData.getListReactions());
		
		Assert.assertEquals(0,ReactionSchemeManipulator.getAllReactions(mzData.getListReactions()).getReactionCount());
		Assert.assertEquals(18, mzData.getListSpectra().getSpectrumElements().size());
		
	}

	/**
	 * CML containing only reactions
	 * 
	 * @throws Exception
	 */
	@Test 
	public void testCMLReaderReactions()throws Exception{
        String pathFile = "src/test/data/cml/multReact.cml";
		InputStream input = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader = new CMLReader(input);
		MZData mzData = new MZData();
		mzData = mzReader.read(mzData);

		Assert.assertNotNull(mzData.getListSpectra());
		Assert.assertNotNull(mzData.getListReactions());
		
		Assert.assertEquals(32,ReactionSchemeManipulator.getAllReactions(mzData.getListReactions()).getReactionCount());
		Assert.assertEquals(0, mzData.getListSpectra().getSpectrumElements().size());
		
	}
	

	/**
	 * CML containing reactions and spectra
	 * 
	 * @throws Exception
	 */
	@Test 
	public void testCMLReaderSpectra_React()throws Exception{
        String pathFile = "src/test/data/cml/spect_react.cml";
		InputStream input = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader = new CMLReader(input);
		MZData mzData = new MZData();
		mzData = mzReader.read(mzData);
		
		Assert.assertEquals(32,ReactionSchemeManipulator.getAllReactions(mzData.getListReactions()).getReactionCount());
		Assert.assertEquals(18, mzData.getListSpectra().getSpectrumElements().size());
		
	}
	
	/**
	 * CML containing spectra with metadata
	 * 
	 * @throws Exception
	 */
	@Test 
	public void testCMLReader_Metadata()throws Exception{
        String pathFile = "src/test/data/cml/spect_react.cml";
        String[] dictRefList = {MZDataConstants.SCANS,MZDataConstants.NUM_GROUPS,MZDataConstants.PROTOCOL,MZDataConstants.INSTRUMENT,MZDataConstants.PARENT_FILE,MZDataConstants.PARENT_COMPOUND,MZDataConstants.MZGAP,MZDataConstants.SNTHRESH,MZDataConstants.RINT,MZDataConstants.ACCURACY,MZDataConstants.RULES,MZDataConstants.ELEMENTS,MZDataConstants.VERSION_SAMS};
		
        InputStream input = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader = new CMLReader(input);
		MZData mzData = new MZData();
		mzData = mzReader.read(mzData);
		
		List<CMLMetadata> metadataList = mzData.getListSpectra().getMetadataListElements().get(0).getMetadataDescendants();
		
		Assert.assertEquals(dictRefList.length,metadataList.size());
		
		Assert.assertEquals(dictRefList[0],metadataList.get(0).getDictRef());
		Assert.assertEquals("18",metadataList.get(0).getContent());
		Assert.assertEquals(mzData.getListSpectra().getSpectrumElements().size(),Integer.parseInt(metadataList.get(0).getContent()));
		
		Assert.assertEquals(dictRefList[1],metadataList.get(1).getDictRef());
		Assert.assertEquals("6",metadataList.get(1).getContent());
		
		Assert.assertEquals(dictRefList[2],metadataList.get(2).getDictRef());
		Assert.assertEquals("leid.pos.v0.1",metadataList.get(2).getContent());

		Assert.assertEquals(dictRefList[3],metadataList.get(3).getDictRef());
		Assert.assertEquals("leiden.Orbi",metadataList.get(3).getContent());

		Assert.assertEquals(dictRefList[4],metadataList.get(4).getDictRef());
		Assert.assertEquals("negHMDB00039.mzXML",metadataList.get(4).getContent());
//		Assert.assertEquals("mzXMLData",metadataList.get(4).getConventionAttribute());

		Assert.assertEquals(dictRefList[5],metadataList.get(5).getDictRef());
		Assert.assertEquals("InChI=1S/C18H32O3/c1-2-3-4-5-7-10-13-16-17(21-16)14-11-8-6-9-12-15-18(19)20/h7,10,16-17H,2-6,8-9,11-15H2,1H3,(H,19,20)/b10-7-",metadataList.get(5).getContent());
		
		Assert.assertEquals(dictRefList[6],metadataList.get(6).getDictRef());
		Assert.assertEquals("0.4",metadataList.get(6).getContent());

		Assert.assertEquals(dictRefList[7],metadataList.get(7).getDictRef());
		Assert.assertEquals("0.2",metadataList.get(7).getContent());

		Assert.assertEquals(dictRefList[8],metadataList.get(8).getDictRef());
		Assert.assertEquals("0.0",metadataList.get(8).getContent());
	}

	/**
	 * CML containing spectra with metadata
	 * 
	 * @throws Exception
	 */
	@Test 
	public void testCMLReader_Spectrum_Conditions()throws Exception{
        String pathFile = "src/test/data/cml/spect_react.cml";
        String[] dictRefList = {MZDataConstants.PRECURSOR_MZ,MZDataConstants.WINDOW,MZDataConstants.ACTIVATION_METHOD,MZDataConstants.POLARITY};
		
        InputStream input = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader = new CMLReader(input);
		MZData mzData = new MZData();
		mzData = mzReader.read(mzData);
		
		CMLSpectrum spect = mzData.getListSpectra().getSpectrumElements().get(0);

		Assert.assertEquals(1,spect.getConditionListElements().size());
		
		CMLConditionList conditionListEl = spect.getConditionListElements().get(0);
		
		Assert.assertEquals(dictRefList.length,conditionListEl.getScalarElements().size());
		
		CMLElements<CMLScalar> conditionScEl = conditionListEl.getScalarElements();
		
		Assert.assertEquals(dictRefList[0],conditionScEl.get(0).getDictRef());
		Assert.assertEquals("0",conditionScEl.get(0).getValue());
		
		Assert.assertEquals(dictRefList[1],conditionScEl.get(1).getDictRef());
		Assert.assertEquals("1",conditionScEl.get(1).getValue());
		
		Assert.assertEquals(dictRefList[2],conditionScEl.get(2).getDictRef());
		Assert.assertEquals("CID",conditionScEl.get(2).getValue());

		Assert.assertEquals(dictRefList[3],conditionScEl.get(3).getDictRef());
		Assert.assertEquals("positive",conditionScEl.get(3).getValue());
	}

	/**
	 * CML containing spectra with metadata
	 * 
	 * @throws Exception
	 */
	@Test 
	public void testCMLReader_Spectrum_Metadata()throws Exception{
        String pathFile = "src/test/data/cml/spect_react.cml";
        String[] dictRefList = {MZDataConstants.SCAN_NUM,MZDataConstants.MS_LEVEL,MZDataConstants.RT,MZDataConstants.GROUP_PEAK_MSN};
		
        InputStream input = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader = new CMLReader(input);
		MZData mzData = new MZData();
		mzData = mzReader.read(mzData);
		
		CMLSpectrum spect = mzData.getListSpectra().getSpectrumElements().get(0);

		Assert.assertEquals(1,spect.getConditionListElements().size());
		
		CMLMetadataList metadataListEl = spect.getMetadataListElements().get(0);
		
		Assert.assertEquals(dictRefList.length,metadataListEl.getMetadataElements().size());
		
		CMLElements<CMLMetadata> metadataScEl = metadataListEl.getMetadataElements();
		
		Assert.assertEquals(dictRefList[1],metadataScEl.get(1).getDictRef());
		Assert.assertEquals("1",metadataScEl.get(1).getContent());
		
		Assert.assertEquals(dictRefList[2],metadataScEl.get(2).getDictRef());
		Assert.assertEquals("607.215",metadataScEl.get(2).getContent());

		Assert.assertEquals(dictRefList[3],metadataScEl.get(3).getDictRef());
		Assert.assertEquals("1",metadataScEl.get(3).getContent());
		
		Assert.assertEquals(dictRefList[0],metadataScEl.get(0).getDictRef());
		Assert.assertEquals("0",metadataScEl.get(0).getContent());
	}
	/**
	 * CML containing spectra with peaks
	 * 
	 * @throws Exception
	 */
	@Test 
	public void testCMLReader_Spectrum_Peaks()throws Exception{
        String pathFile = "src/test/data/cml/spect_react.cml";
		
        InputStream input = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader = new CMLReader(input);
		MZData mzData = new MZData();
		mzData = mzReader.read(mzData);
		
		CMLSpectrum spect = mzData.getListSpectra().getSpectrumElements().get(0);
		
		Assert.assertEquals(1,spect.getPeakListElements().size());
		
		CMLPeak peakE = spect.getPeakListElements().get(0).getPeakElements().get(0);
		
		Assert.assertEquals("1.1",peakE.getId());
		Assert.assertEquals(188.07041931152344,peakE.getXValue(),0.001);
		Assert.assertEquals(42010.92578125,peakE.getYValue(),0.001);
		Assert.assertEquals("1",peakE.getMoleculeRefs()[0]);
	}
	

	/**
	 * CML containing spectra with peaks
	 * 
	 * @throws Exception
	 */
	@Test 
	public void testReadProperties()throws Exception{
        String pathFile = "src/test/data/cml/spect_react.cml";
		
        InputStream input = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader = new CMLReader(input);
		MZData mzData = new MZData();
		mzData = mzReader.read(mzData);
		
		Map<Object, Object> properties = mzData.getProperties();
		Assert.assertEquals("0.4", properties.get(MZDataConstants.MZGAP));
		Assert.assertEquals("0.2", properties.get(MZDataConstants.SNTHRESH));
		Assert.assertEquals("0.0", properties.get(MZDataConstants.RINT));
		Assert.assertEquals("6", properties.get(MZDataConstants.NUM_GROUPS));
		Assert.assertEquals("18", properties.get(MZDataConstants.SCANS));
		Assert.assertEquals("0.14", properties.get(MZDataConstants.VERSION_SAMS));
		Assert.assertEquals("15,15,15,15,15", properties.get(MZDataConstants.ACCURACY));
		Assert.assertEquals("nitrogenR,RDBER", properties.get(MZDataConstants.RULES));
		Assert.assertEquals("C1..50,H1..100,N0..30,O1..30", properties.get(MZDataConstants.ELEMENTS));
		Assert.assertEquals("InChI=1S/C18H32O3/c1-2-3-4-5-7-10-13-16-17(21-16)14-11-8-6-9-12-15-18(19)20/h7,10,16-17H,2-6,8-9,11-15H2,1H3,(H,19,20)/b10-7-", properties.get(MZDataConstants.PARENT_COMPOUND));
		
	}
	
	/**
	 * CML containing spectra with peaks
	 * 
	 * @throws Exception
	 */
	@Test 
	public void testReadPropertiesGroup()throws Exception{
        String pathFile = "src/test/data/cml/spect_react.cml";
		
        InputStream input = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader = new CMLReader(input);
		MZData mzData = new MZData();
		mzData = mzReader.read(mzData);
		mzData = MZDataManipulator.group(mzData, 0.4);
		
		Map<Object, Object> properties = mzData.getProperties();
		Assert.assertEquals("0.4", properties.get(MZDataConstants.MZGAP));
		Assert.assertEquals("0.2", properties.get(MZDataConstants.SNTHRESH));
		Assert.assertEquals("0.0", properties.get(MZDataConstants.RINT));
		Assert.assertEquals("1", properties.get(MZDataConstants.NUM_GROUPS));
		Assert.assertEquals("18", properties.get(MZDataConstants.SCANS));
		Assert.assertEquals("0.14", properties.get(MZDataConstants.VERSION_SAMS));
		Assert.assertEquals("15,15,15,15,15", properties.get(MZDataConstants.ACCURACY));
		Assert.assertEquals("nitrogenR,RDBER", properties.get(MZDataConstants.RULES));
		Assert.assertEquals("C1..50,H1..100,N0..30,O1..30", properties.get(MZDataConstants.ELEMENTS));
		Assert.assertEquals("InChI=1S/C18H32O3/c1-2-3-4-5-7-10-13-16-17(21-16)14-11-8-6-9-12-15-18(19)20/h7,10,16-17H,2-6,8-9,11-15H2,1H3,(H,19,20)/b10-7-", properties.get(MZDataConstants.PARENT_COMPOUND));
		
	}
	
	/**
	 * Read formula properties of the an atomContainer from CML containing 
	 * spectra with peaks and atomcontainer information
	 * 
	 * bug: SAMS-30
	 * 
	 * @throws Exception
	 */
	@Test 
	public void testReadPropertiesMolecules()throws Exception{
        String pathFile = "src/test/data/cml/WithInChIs.cml";
        String id107 = "107";
        String formula107 = "C13H9O";
		
        InputStream input = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader = new CMLReader(input);
		MZData mzData = new MZData();
		mzData = mzReader.read(mzData);
		
		IReactionScheme rea = mzData.getListReactions();
		IMoleculeSet allM = ReactionSchemeManipulator.getAllMolecules(rea);
		for(IAtomContainer mol : allM.molecules()){
			if(id107.equals(mol.getID())){
				IMolecularFormula formMol = (IMolecularFormula)mol.getProperty(CDKConstants.FORMULA);
				Assert.assertEquals(formula107,MolecularFormulaManipulator.getString(formMol));
			}
		}
	}

	
	/**
	 * Read formula properties of the an atomContainer from CML containing 
	 * spectra with peaks and atomcontainer information. Applied the grouping.
	 * 
	 * bug: SAMS-30
	 * 
	 * @throws Exception
	 */
	@Test 
	public void testReadPropertiesMoleculesGroup()throws Exception{
        String pathFile = "src/test/data/cml/WithFormulas.cml";
        String id107 = "107";
        String formula107 = "C13H9O";
		
        InputStream input = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader = new CMLReader(input);
		MZData mzData = new MZData();
		mzData = mzReader.read(mzData);
		
		MZData mzDataGr = MZDataManipulator.group(mzData, 0.4);
		IReactionScheme rea = mzDataGr.getListReactions();
		IMoleculeSet allM = ReactionSchemeManipulator.getAllMolecules(rea);
		for(IAtomContainer mol : allM.molecules()){
			if(id107.equals(mol.getID())){
				IMolecularFormula formMol = (IMolecularFormula)mol.getProperty(CDKConstants.FORMULA);
				Assert.assertEquals(formula107,MolecularFormulaManipulator.getString(formMol));
			}
		}
	}
	
	/**
	 * Read formula properties of the an atomContainer from CML containing 
	 * spectra with peaks and atomcontainer information. One of the fragments
	 * contains multiple formula added
	 * 
	 * bug: SAMS-39
	 * 
	 * @throws Exception
	 */
	@Test 
	public void testReadPropertiesMoleculesMultiFormul()throws Exception{
        String pathFile = "src/test/data/cml/WithFormulasMulti.cml";
        String id107 = "107";
        String formula107A = "C13H9O";
        String formula107B = "C13H9O2";
		
        InputStream input = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader = new CMLReader(input);
		MZData mzData = new MZData();
		mzData = mzReader.read(mzData);
		
		IReactionScheme rea = mzData.getListReactions();
		IMoleculeSet allM = ReactionSchemeManipulator.getAllMolecules(rea);
		for(IAtomContainer mol : allM.molecules()){
			if(id107.equals(mol.getID())){
				IMolecularFormulaSet formMolSet = (IMolecularFormulaSet)mol.getProperty(CDKConstants.FORMULA);
				IMolecularFormula formA = formMolSet.getMolecularFormula(0);
				IMolecularFormula formB = formMolSet.getMolecularFormula(1);
				Assert.assertEquals(formula107A,MolecularFormulaManipulator.getString(formA));
				Assert.assertEquals(formula107B,MolecularFormulaManipulator.getString(formB));
			}
		}
	}
	
	@Test 
	public void testReadPropertiesMoleculesSpectra()throws Exception{
		String pathFile = "src/test/data/cml/WithInChIs.cml";
        String id107 = "107";
        String parentID = "101";
        String mass = "181.06541442871094";
        String intensity = "156335.0";
        String inchi = "InChI=1S2/C15H10O4/c16-10-7-5-9(6-8-10)15-14(18)13(17)11-3-1-2-4-12(11)19-15/h1-8,16,18H";
        String group = "4";
        String comment = "That is a comment";
        
        InputStream input = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader = new CMLReader(input);
		MZData mzData = new MZData();
		mzData = mzReader.read(mzData);
		
		
		IReactionScheme rea = mzData.getListReactions();
		IMoleculeSet allM = ReactionSchemeManipulator.getAllMolecules(rea);
		for(IAtomContainer mol : allM.molecules()){
			if(id107.equals(mol.getID())){
				Assert.assertEquals(id107,mol.getID());
				Assert.assertEquals(group,mol.getProperty(MZDataConstants.GROUP_PEAK_MSN));
				Assert.assertEquals(mass,mol.getProperty(MZDataConstants.MASS));
				Assert.assertEquals(intensity,mol.getProperty(MZDataConstants.INTENSITIY));
				Assert.assertEquals(inchi,mol.getProperty(CDKConstants.INCHI));
				Assert.assertEquals(parentID,mol.getProperty(MZDataConstants.PRECURSOR_ID));
				Assert.assertEquals(comment,mol.getProperty(CDKConstants.ANNOTATIONS));
			}
		}		
	}

	
	@Test 
	public void testReadPropertiesMoleculesSpectraGroup()throws Exception{
		String pathFile = "src/test/data/cml/WithInChIs.cml";
        String id107 = "107";
        String parentID = "101";
        String mass = "181.06541442871094";
        String intensity = "156335.0";
        String inchi = "InChI=1S2/C15H10O4/c16-10-7-5-9(6-8-10)15-14(18)13(17)11-3-1-2-4-12(11)19-15/h1-8,16,18H";
        String group = "1";
        String comment = "That is a comment";
        
        InputStream input = new BufferedInputStream(new FileInputStream(pathFile));
		CMLReader mzReader = new CMLReader(input);
		MZData mzData = new MZData();
		mzData = mzReader.read(mzData);

		MZData mzDataGr = MZDataManipulator.group(mzData, 0.4);
		IReactionScheme rea = mzDataGr.getListReactions();
		IMoleculeSet allM = ReactionSchemeManipulator.getAllMolecules(rea);
		for(IAtomContainer mol : allM.molecules()){
			if(id107.equals(mol.getID())){
				Assert.assertEquals(id107,mol.getID());
				Assert.assertEquals(group,mol.getProperty(MZDataConstants.GROUP_PEAK_MSN));
				Assert.assertEquals(mass,mol.getProperty(MZDataConstants.MASS));
				Assert.assertEquals(intensity,mol.getProperty(MZDataConstants.INTENSITIY));
				Assert.assertEquals(inchi,mol.getProperty(CDKConstants.INCHI));
				Assert.assertEquals(parentID,mol.getProperty(MZDataConstants.PRECURSOR_ID));
				Assert.assertEquals(comment,mol.getProperty(CDKConstants.ANNOTATIONS));
			}
		}
		
	}
}

