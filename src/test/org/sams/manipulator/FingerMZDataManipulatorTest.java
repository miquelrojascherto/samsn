package org.sams.manipulator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import junit.framework.Assert;
import junit.framework.TestCase;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.sams.MZData;

public class FingerMZDataManipulatorTest  extends TestCase {
	static DefaultChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();

	/**
	 *  Constructor for the FingerMZDataManipulatorTest object.
	 */
	public FingerMZDataManipulatorTest(){
		
	}
	
	/**
	 * @throws Exception
	 */
	public void testCompare()throws Exception{
		List<String> value1 = new ArrayList<String>();
		value1.add("C6H15N2O3@C6H10NO2");
		value1.add("C6H15N2O3@C6H10NO1");
		value1.add("C6H15N2O3@C6H10NO2@C5H8N");
		value1.add("C6H15N2O3@C6H10NO2@C5H9N");
		value1.add("C6H15N2O3@C6H11N2O@C5H11N2");
		value1.add("C6H15N2O3@C6H11N2O@C5H12N2");
		
		List<String> value2 = new ArrayList<String>();
		value2.add("C6H15N2O3@C6H10NO2");
		value2.add("C6H15N2O3@C6H10NO1");
		
		Assert.assertEquals(1.00, FingerMZDataManipulator.compare(value1, value2),0.01);
		Assert.assertEquals(0.33, FingerMZDataManipulator.compare(value2, value1),0.01);

	}
	/**
	 * @throws Exception
	 */
	public void testCompare2()throws Exception{
		List<String> value1 = new ArrayList<String>();
		value1.add("66");
		value1.add("67");
		value1.add("92");
		value1.add("105");
		value1.add("111");
		value1.add("145");
		value1.add("151");
		value1.add("160");
		value1.add("161");
		value1.add("162");
		value1.add("163");
		value1.add("164");
		value1.add("165");
		
		List<String> value2 = new ArrayList<String>();
		value2.add("66");
		value2.add("87");
		value2.add("92");
		value2.add("105");
		value2.add("132");
		value2.add("144");
		value2.add("151");
		value2.add("160");
		value2.add("161");
		value2.add("162");

		Assert.assertEquals(0.7, FingerMZDataManipulator.compare(value1, value2),0.01);
		Assert.assertEquals(0.53, FingerMZDataManipulator.compare(value2, value1),0.01);

	}
	/**
	 * @throws Exception
	 */
	public void testTanimmoto()throws Exception{
		List<String> value1 = new ArrayList<String>();
		value1.add("66");
		value1.add("67");
		value1.add("92");
		value1.add("105");
		value1.add("111");
		value1.add("145");
		value1.add("151");
		value1.add("160");
		value1.add("161");
		value1.add("162");
		value1.add("163");
		value1.add("164");
		value1.add("165");
		
		List<String> value2 = new ArrayList<String>();
		value2.add("66");
		value2.add("87");
		value2.add("92");
		value2.add("105");
		value2.add("132");
		value2.add("144");
		value2.add("151");
		value2.add("160");
		value2.add("161");
		value2.add("162");
		
		Assert.assertEquals(0.4375, FingerMZDataManipulator.tanimotto(value1, value2),0.0001);
		Assert.assertEquals(0.4375, FingerMZDataManipulator.tanimotto(value2, value1),0.0001);
	}
	/**
	 * Bug FP-3
	 * 
	 * @throws Exception
	 */
	public void testTanimmoto2()throws Exception{
		List<String> value1 = new ArrayList<String>();
		value1.add("C6H15N2O3");
		value1.add("C6H15N2O3@C6H10NO2");
		value1.add("C6H15N2O3@C6H10NO2@C5H8N");
		value1.add("C6H15N2O3@C6H10NO2@C5H9N");
		value1.add("C6H15N2O3@C6H11N2O");
		value1.add("C6H15N2O3@C6H11N2O@C5H11N2");
		value1.add("C6H15N2O3@C6H11N2O@C5H12N2");
		
		List<String> value2 = new ArrayList<String>();
		value2.add("C6H15N2O3");
		value2.add("C6H15N2O3@C6H10NO2");
		value2.add("C6H15N2O3@C6H10NO2@C5H8N");
		value2.add("C6H15N2O3@C6H10NO2@C5H9N");
		
		Assert.assertEquals(0.66, FingerMZDataManipulator.tanimotto2(value1, value2),0.01);
		Assert.assertEquals(0.66, FingerMZDataManipulator.tanimotto2(value2, value1),0.01);

	}

	/**
	 * comparing with its self
	 * 
	 * @throws Exception
	 */
	public void testTanimmoto2a()throws Exception{
		List<String> value1 = new ArrayList<String>();
		value1.add("C6H15N2O3");
		value1.add("C6H15N2O3@C6H10NO2");
		value1.add("C6H15N2O3@C6H10NO2@C5H8N");
		value1.add("C6H15N2O3@C6H10NO2@C5H9N");
		value1.add("C6H15N2O3@C6H11N2O");
		value1.add("C6H15N2O3@C6H11N2O@C5H11N2");
		value1.add("C6H15N2O3@C6H11N2O@C5H12N2");
		
		
		Assert.assertEquals(1.0, FingerMZDataManipulator.tanimotto2(value1, value1),0.01);

	}
	/**
	 * comparing with its self
	 * 
	 * @throws Exception
	 */
	public void testTanimmoto3A()throws Exception{
		List<String> value1 = new ArrayList<String>();
		value1.add("A@B");
		value1.add("B@C");
		value1.add("A@D");
		value1.add("D@E");
		value1.add("A@B@C");
		value1.add("A@D@E");
		
		List<String> value2 = new ArrayList<String>();
		value2.add("A@B");
		value2.add("B@C");
		value2.add("A@H");
		value2.add("H@I");
		value2.add("A@B@C");
		value2.add("A@H@I");

		//Finall features: AD, DE, ABC, AH, HI
		Assert.assertEquals(0.2, FingerMZDataManipulator.tanimotto3(value1, value2),0.01);
		Assert.assertEquals(0.2, FingerMZDataManipulator.tanimotto3(value2, value1),0.01);

	}

	/**
	 * comparing with different order
	 * 
	 * @throws Exception
	 * @sams-bug SAMS-16
	 */
	public void testTanimmoto3B()throws Exception{
		String[] valu1A = {"171@143","171@153","171@[143||153]","181@153","18@28","18@44","18@72","18@[72||44]","199@171","199@171@143","199@171@153","199@181","199@181@153","199@[171||181]","209@153","209@181","209@181@153","209@[153||181]","237@165","237@193","237@[165||193]","255@133","255@199","255@199@171","255@199@171@143","255@199@171@153","255@199@181","255@199@181@153","255@209","255@209@153","255@209@181","255@209@181@153","255@237","255@237@165","255@237@193","255@[199||133]","255@[199||209]","255@[199||237]","255@[209||133]","255@[209||237]","255@[237||133]","28@18","28@28","28@[28||18]","46@28","46@28@28","46@56","46@[56||28]","56@18","56@18@28","56@28","56@28@18","56@28@28","56@[28||18]"};
		List<String> value1a = FingerMZDataManipulator.removeDuplicates(Arrays.asList(valu1A));
		String[] value2A = {"171@153","181@153","199@171","199@171@153","199@181","199@181","199@181@153","199@181@153","199@[171||181]","199@[171||181]","209@153","209@153","209@181","209@181","209@[153||181]","209@[153||181]","237@165","237@165","237@193","237@193","237@[165||193]","237@[165||193]","255@133","255@133","255@199","255@199","255@199@171","255@199@171","255@199@171@153","255@199@171@153","255@199@181","255@199@181","255@199@181@153","255@199@181@153","255@209","255@209","255@209@153","255@209@153","255@209@181","255@209@181","255@237","255@237","255@237@165","255@237@165","255@237@193","255@237@193","255@[133||199]","255@[133||199]","255@[133||209]","255@[133||209]","255@[133||237]","255@[133||237]","255@[199||209]","255@[199||209]","255@[199||237]","255@[199||237]","255@[209||237]","255@[209||237]"};
		List<String> value2a = Arrays.asList(value2A);

		String[] valu1B = {"171@153","171@[143||153]","181@153","18@28","18@44","18@72","18@[72||44]","199@171","199@171@143","199@171@153","199@181","199@181@153","199@[171||181]","209@153","209@181","209@181@153","209@[153||181]","237@165","237@193","237@[165||193]","255@133","255@199","255@199@171","255@199@171@143","255@199@171@153","255@199@181","255@199@181@153","255@209","255@209@153","255@209@181","255@209@181@153","255@237","255@237@165","255@237@193","255@[199||133]","255@[199||209]","255@[199||237]","255@[209||133]","255@[209||237]","255@[237||133]","28@18","28@28","28@[28||18]","46@28","46@28@28","46@56","46@[56||28]","56@18","56@18@28","56@28","56@28@18","56@28@28","56@[28||18]","171@143"};
		List<String> value1b = FingerMZDataManipulator.removeDuplicates(Arrays.asList(valu1B));
		
		double result1 = FingerMZDataManipulator.tanimotto3(value1a, value2a);
		double result2 = FingerMZDataManipulator.tanimotto3(value1b, value2a);
		Assert.assertEquals(result1, result2,0.0001);

	}

	/**
	 * comparing with different order
	 * 
	 * @throws Exception
	 * @sams-bug SAMS-16
	 */
	public void testTanimmoto3C()throws Exception{
		String[] valu1A = {"171@143","171@153","171@[143||153]","181@153","18@28","18@44","18@72","18@[72||44]","199@171","199@171@143","199@171@153","199@181","199@181@153","199@[171||181]","209@153","209@181","209@181@153","209@[153||181]","237@165","237@193","237@[165||193]","255@133","255@199","255@199@171","255@199@171@143","255@199@171@153","255@199@181","255@199@181@153","255@209","255@209@153","255@209@181","255@209@181@153","255@237","255@237@165","255@237@193","255@[199||133]","255@[199||209]","255@[199||237]","255@[209||133]","255@[209||237]","255@[237||133]","28@18","28@28","28@[28||18]","46@28","46@28@28","46@56","46@[56||28]","56@18","56@18@28","56@28","56@28@18","56@28@28","56@[28||18]"};
		List<String> value1a = FingerMZDataManipulator.removeDuplicates(Arrays.asList(valu1A));
		String[] value2A = {"171@153","181@153","199@171","199@171@153","199@181","199@181","199@181@153","199@181@153","199@[171||181]","199@[171||181]","209@153","209@153","209@181","209@181","209@[153||181]","209@[153||181]","237@165","237@165","237@193","237@193","237@[165||193]","237@[165||193]","255@133","255@133","255@199","255@199","255@199@171","255@199@171","255@199@171@153","255@199@171@153","255@199@181","255@199@181","255@199@181@153","255@199@181@153","255@209","255@209","255@209@153","255@209@153","255@209@181","255@209@181","255@237","255@237","255@237@165","255@237@165","255@237@193","255@237@193","255@[133||199]","255@[133||199]","255@[133||209]","255@[133||209]","255@[133||237]","255@[133||237]","255@[199||209]","255@[199||209]","255@[199||237]","255@[199||237]","255@[209||237]","255@[209||237]"};
		List<String> value2a = FingerMZDataManipulator.removeDuplicates(Arrays.asList(value2A));

		double result1 = FingerMZDataManipulator.tanimotto3(value1a, value2a);
		double result2 = FingerMZDataManipulator.tanimotto3(value2a, value1a);
		Assert.assertEquals(result1, result2,0.0001);

	}
	
	public void testGetDuplicates() throws Exception{

		String[] arr1 = {"A","B","C","C","D"};
		List<String> list1 = Arrays.asList(arr1);
		Set<String> set1 = new HashSet<String>(list1);

		String[] arr2 = {"A","B","C","D"};
		List<String> list2 = Arrays.asList(arr2);
		Set<String> set2 = new HashSet<String>(list2);
		
		Assert.assertTrue(set1.size() < list1.size());
		Assert.assertFalse(set2.size() < list2.size());
	}
	
	/**
	 * Testing features are unique
	 * 
	 * @throws Exception
	 * @sams-bug SAMS-17
	 */
	public void testGetbs24_24_1()throws Exception{
		
		MZData mzdata = new MZData();
		List<String> list = FingerMZDataManipulator.getbs24_24_1(mzdata);
		Set<String> set = new HashSet<String>(list);
		Assert.assertFalse(set.size() < list.size());
	}
	/**
	 * Testing features are unique
	 * 
	 * @throws Exception
	 * @sams-bug SAMS-17
	 */
	public void testGetbs24_24_1Nom()throws Exception{

		MZData mzdata = new MZData();
		List<String> list = FingerMZDataManipulator.getbs24_24_1Nom(mzdata);
		Set<String> set = new HashSet<String>(list);
		Assert.assertFalse(set.size() < list.size());
	}

	/**
	 * Testing features are unique
	 * 
	 * @throws Exception
	 * @sams-bug SAMS-17
	 */
	public void testGetNominalPath()throws Exception{

		MZData mzdata = new MZData();
		List<String> list0 = new ArrayList<String>();
		List<String> list = FingerMZDataManipulator.getNominalPath(list0,"#");
		Set<String> set = new HashSet<String>(list);
		Assert.assertFalse(set.size() < list.size());
	}
}
