package org.sams.spect;
import org.junit.Assert;
import org.junit.Test;
import org.sams.SAMSTestCase;

/**
 * TestCase for the ScanMZXML class.
 * 
 * @author Miguel Rojas-Cherto
 */
public class ScanMZXMLTest extends SAMSTestCase {
	String num = "1";
	String polarity = "+1";
	String collisionEnergy = "34";

	/**
	 *  Constructor for the ScanMZXMLTest object.
	 */
	public ScanMZXMLTest(){
		
	}
	/**
	 * 
	 * @throws Exception
	 */
	@Test 
	public void testNull() throws Exception{
		ScanMZXML e = new ScanMZXML(null,null,null);
		Assert.assertNotNull(e);
		
	}

	@Test 
	public void testGetNum() throws Exception{
		ScanMZXML e = new ScanMZXML(num,null,null);
		Assert.assertEquals(num, e.getNum());
		
	}
	@Test 
	public void testGetpolarity() throws Exception{
		ScanMZXML e = new ScanMZXML(null,polarity,null);
		Assert.assertEquals(polarity, e.getpolarity());
		
	}
	@Test 
	public void testGetCollisionEnergy() throws Exception{
		ScanMZXML e = new ScanMZXML(null,null,collisionEnergy);
		Assert.assertEquals(collisionEnergy, e.getcollisionEnergy());
		
	}
}