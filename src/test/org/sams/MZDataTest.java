package org.sams;

import java.util.Map;

import junit.framework.Assert;

import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;

public class MZDataTest extends SAMSTestCase {
	static DefaultChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
	
	/**
	 *  Constructor for the MZDataTest object.
	 */
	public MZDataTest(){
		
	}
	
	@Test 
	public void testClone()throws Exception{
		
	}
	
	@Test 
	public void testGetPropertiesNull()throws Exception{
		MZData mzData = new MZData();
		
		Map<Object, Object> properties = mzData.getProperties();
		Assert.assertEquals(0, properties.size());
		
	}
	
	@Test 
	public void testSetProperties()throws Exception{
		MZData mzData = new MZData();
		
		mzData.setProperty(MZDataConstants.NUM_GROUPS, "1");
		Assert.assertEquals(1, mzData.getProperties().size());
		
		mzData.setProperty(MZDataConstants.PARENT_COMPOUND, "InChI=bla");
		Assert.assertEquals(2, mzData.getProperties().size());
		
	}

	@Test 
	public void testUniqueSetProperties()throws Exception{
		MZData mzData = new MZData();
		
		mzData.setProperty(MZDataConstants.PARENT_COMPOUND, "InChI=bla");
		Assert.assertEquals(1, mzData.getProperties().size());
		Assert.assertEquals("InChI=bla", mzData.getProperty(MZDataConstants.PARENT_COMPOUND));
		
		mzData.setProperty(MZDataConstants.PARENT_COMPOUND, "InChI=bla2");
		Assert.assertEquals(1, mzData.getProperties().size());
		Assert.assertEquals("InChI=bla2", mzData.getProperty(MZDataConstants.PARENT_COMPOUND));
		
	}
}
