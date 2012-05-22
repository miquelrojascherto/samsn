package org.sams.io;

import java.io.File;
import java.io.FileInputStream;
import java.util.List;
import java.util.Map;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.tools.manipulator.ReactionSchemeManipulator;
import org.sams.MZData;
import org.sams.MZDataConstants;
import org.sams.SAMSTestCase;
import org.xmlcml.cml.element.CMLMetadata;

/**
 * TestCase for the MZXMLReader class.
 * 
 * @author Miguel Rojas-Cherto
 */
public class MZXMLReaderTest extends SAMSTestCase {

	/**
	 *  Constructor for the CMLWriterTest object.
	 */
	public MZXMLReaderTest(){
		
	}
	/**
	 * @throws Exception
	 */
	@Test 
	public void testMZXMLReader()throws Exception{
		String pathFile = "src/test/data/mzXML/Ade_AM_1.mzXML";
		MZXMLReader mzReader = new MZXMLReader(new FileInputStream(new File(pathFile)));
		MZData mzData = new MZData();
		
		mzReader.setMZGap(0.2);
		mzReader.setSNThresh(1.0);
		mzReader.setRInt(0.4);
		mzData = mzReader.read(mzData);
		
		Assert.assertNotNull(mzData.getListSpectra());
		Assert.assertNotNull(mzData.getListReactions());
		
		Assert.assertEquals(0,ReactionSchemeManipulator.getAllReactions(mzData.getListReactions()).getReactionCount());
		Assert.assertEquals(44, mzData.getListSpectra().getSpectrumElements().size());
		
	}

	/**
	 * @throws Exception
	 */
	@Test 
	public void testReadProperties()throws Exception{
		String pathFile = "src/test/data/mzXML/Ade_AM_1.mzXML";
		MZXMLReader mzReader = new MZXMLReader(new FileInputStream(new File(pathFile)));
		MZData mzData = new MZData();
		
		double mzgap = 0.2;
		mzReader.setMZGap(mzgap);
		double snthresh = 1.0;
		mzReader.setSNThresh(snthresh);
		double rint = 0.4;
		mzReader.setRInt(rint);
		mzData = mzReader.read(mzData);
		
		Map<Object, Object> properties = mzData.getProperties();
		Assert.assertEquals(mzgap, properties.get(MZDataConstants.MZGAP));
		Assert.assertEquals(snthresh, properties.get(MZDataConstants.SNTHRESH));
		Assert.assertEquals(rint, properties.get(MZDataConstants.RINT));
		Assert.assertEquals(15, properties.get(MZDataConstants.NUM_GROUPS));
		Assert.assertEquals(44, properties.get(MZDataConstants.SCANS));
		
	}


	/**
	 * @throws Exception
	 */
	@Test 
	public void testReadMetadata()throws Exception{
		String[] dictRefList = {MZDataConstants.SCANS,MZDataConstants.NUM_GROUPS,MZDataConstants.PROTOCOL,MZDataConstants.INSTRUMENT,MZDataConstants.PARENT_FILE,MZDataConstants.MZGAP,MZDataConstants.SNTHRESH,MZDataConstants.RINT};
		String pathFile = "src/test/data/mzXML/Ade_AM_1.mzXML";
		MZXMLReader mzReader = new MZXMLReader(new FileInputStream(new File(pathFile)));
		MZData mzData = new MZData();
		
		double mzgap = 0.2;
		mzReader.setMZGap(mzgap);
		double snthresh = 1.0;
		mzReader.setSNThresh(snthresh);
		double rint = 0.4;
		mzReader.setRInt(rint);
		mzData = mzReader.read(mzData);
		
		List<CMLMetadata> metadataList = mzData.getListSpectra().getMetadataListElements().get(0).getMetadataDescendants();
		
		Assert.assertEquals(6,metadataList.size());
		
		Assert.assertEquals(dictRefList[0],metadataList.get(0).getDictRef());
		Assert.assertEquals("44",metadataList.get(0).getContent());
		
		Assert.assertEquals(dictRefList[1],metadataList.get(1).getDictRef());
		Assert.assertEquals("15",metadataList.get(1).getContent());
		
		
	}
	/**
	 * @throws Exception
	 */
	@Test 
	public void testReadEmpty()throws Exception{
		String pathFile = "src/test/data/mzXML/HMDB01939_neg02.mzXML";
		MZXMLReader mzReader = new MZXMLReader(new FileInputStream(new File(pathFile)));
		MZData mzData = new MZData();
		
		double mzgap = 0.2;
		mzReader.setMZGap(mzgap);
		double snthresh = 1.0;
		mzReader.setSNThresh(snthresh);
		double rint = 0.4;
		mzReader.setRInt(rint);
		mzData = mzReader.read(mzData);
		
		List<CMLMetadata> metadataList = mzData.getListSpectra().getMetadataListElements().get(0).getMetadataDescendants();
		
		Assert.assertEquals(6,metadataList.size());
		
		Assert.assertNotNull(mzData.getListSpectra());
		Assert.assertNotNull(mzData.getListReactions());
		
		Assert.assertEquals(0,ReactionSchemeManipulator.getAllReactions(mzData.getListReactions()).getReactionCount());
		Assert.assertEquals(0, mzData.getListSpectra().getSpectrumElements().size());
		
		
	}
}
