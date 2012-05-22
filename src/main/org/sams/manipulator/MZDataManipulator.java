package org.sams.manipulator;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.ReactionScheme;
import org.openscience.cdk.formula.IsotopeContainer;
import org.openscience.cdk.formula.IsotopePattern;
import org.openscience.cdk.formula.IsotopePatternGenerator;
import org.openscience.cdk.formula.IsotopePatternSimilarity;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.interfaces.IMolecularFormulaSet;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IMoleculeSet;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.IReactionScheme;
import org.openscience.cdk.interfaces.IReactionSet;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.openscience.cdk.tools.manipulator.ReactionManipulator;
import org.openscience.cdk.tools.manipulator.ReactionSchemeManipulator;
import org.sams.FingerMZData;
import org.sams.MEFgenerator.Polarity;
import org.sams.MZData;
import org.sams.MZDataConstants;
import org.sams.main.PackageVersion;
import org.sams.spect.ParentIon;
import org.sams.spect.ScanMZXML;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;
import org.xmlcml.cml.base.CMLElements;
import org.xmlcml.cml.element.CMLConditionList;
import org.xmlcml.cml.element.CMLMetadata;
import org.xmlcml.cml.element.CMLMetadataList;
import org.xmlcml.cml.element.CMLPeak;
import org.xmlcml.cml.element.CMLPeakList;
import org.xmlcml.cml.element.CMLScalar;
import org.xmlcml.cml.element.CMLSpectrum;
import org.xmlcml.cml.element.CMLSpectrumList;
/**
 * Class to manipulate MZData objects.
 * 
 * @author Miguel Rojas-Cherto
 */
public class MZDataManipulator {

	static DefaultChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
	private static double minAbISO = 0.01;
	private static double toleranceISO = 4;
	/**
	 * Grouping the different spectral trees based on the elemental formula.
	 * 
	 * @param mzData   The mzData object
	 * @param minOccur The minimun occurrence to accept
	 * @return         The mzData object grouped by the occurrence
	 */
	public static MZData group(MZData mzData, double minOccur) {
		List<MZData> mzDataList = new ArrayList<MZData>();
		mzDataList.add(mzData);
		
		Map<String, List<IAtomContainer>> listPath = MZDataManipulator.groupA(mzDataList,minOccur);
		MZData newMZD = MZDataManipulator.getMZData(listPath,mzDataList);
		newMZD.setProperties(mzData.getProperties());
		newMZD.setProperty(MZDataConstants.NUM_GROUPS, "1");
		return newMZD;
	}
	/**
	 * Grouping the different spectral trees based on the nominal mass.
	 * 
	 * @param mzData   The mzData object
	 * @param minOccur The minimun occurrence to accept
	 * @return         The mzData object grouped by the occurrence
	 */
	public static MZData groupNominal(MZData mzData, double minOccur) {
		Map<String, List<CMLPeak>> listPath = MZDataManipulator.groupN(mzData,minOccur);
		MZData newMZD = MZDataManipulator.getMZDataNominal(listPath,mzData);
		newMZD.setProperties(mzData.getProperties());
		newMZD.setProperty(MZDataConstants.NUM_GROUPS, "1");
		return newMZD;
	}
	/**
	 * Grouping the different spectral trees based on the elemental formula.
	 * 
	 * @param mzDataList   A List with mzData objects
	 * @param minOccur     The minimun occurrence to accept
	 * @return             The mzData object grouped by the occurrence
	 */
	public static MZData group(List<MZData> mzDataList, double minOccur) {
		Map<String, List<IAtomContainer>> listPath = MZDataManipulator.groupA(mzDataList,minOccur);
		MZData newMZD = MZDataManipulator.getMZData(listPath,mzDataList);
		newMZD.setProperty(MZDataConstants.NUM_GROUPS, "1");
		return newMZD;
	}
	/**
	 * At the moment only return a List of Strings
	 * 
	 * @param mzData
	 * @param minOccur
	 * @return
	 */
	private static Map<String, List<IAtomContainer>> groupA(List<MZData> mzDataList, double minOccur) {
		HashMap<String, List<IAtomContainer>> hashMapAC = new HashMap<String, List<IAtomContainer>>();
		
		List<String> listPath = new ArrayList<String>();
		for(MZData mzData:mzDataList)
			listPath.addAll(getListPath(mzData));
		
		int countL1 = 0;
		for(String path:listPath){
			int levelI = path.split("\\|\\|").length;
//			System.out.println(levelI+" "+path);
			if(levelI == 1)
				countL1 ++;
			
		}
		
		Map<String, List<IAtomContainer>> hashMap = getHashMapPath(mzDataList);
		Set<Entry<String, List<IAtomContainer>>> set = hashMap.entrySet();
		for( Iterator<Entry<String, List<IAtomContainer>>> it = set.iterator(); it.hasNext();){
			Entry<String, List<IAtomContainer>> key = it.next();
			List<IAtomContainer> ll = key.getValue();
			String keyS = key.getKey();
			double value = ((double)ll.size())/((double)countL1);
//			System.out.println(keyS+" "+ll.size()+"/"+countL1+" > "+value);
			if(value < minOccur){
//				System.out.println("   R");
//				hashMap.remove(keyS);
			}else{
				hashMapAC.put(keyS, ll);
			}
			
		}
		// sort
		Map<String, List<IAtomContainer>> hashMapS = new TreeMap<String, List<IAtomContainer>>(hashMapAC);
		return hashMapS;
	}

	/**
	 * At the moment only return a List of Strings
	 * 
	 * @param mzData
	 * @param minOccur
	 * @return
	 */
	private static Map<String, List<CMLPeak>> groupN(MZData mzData, double minOccur) {
		HashMap<String, List<CMLPeak>> hashMapAC = new HashMap<String, List<CMLPeak>>();
		List<String> listPath = getListPathNominal(mzData);
		int countL1 = 0;
		for(String path:listPath){
			int levelI = path.split("\\|\\|").length;
			if(levelI == 1)
				countL1 ++;
			
		}
		
		Map<String, List<CMLPeak>> hashMap = getHashMapPathN(mzData);
		Set<Entry<String, List<CMLPeak>>> set = hashMap.entrySet();
		for( Iterator<Entry<String, List<CMLPeak>>> it = set.iterator(); it.hasNext();){
			Entry<String, List<CMLPeak>> key = it.next();
			List<CMLPeak> ll = key.getValue();
			String keyS = key.getKey();
			double value = ((double)ll.size())/((double)countL1);
//			System.out.println(keyS+" "+ll.size()+"/"+countL1+" > "+value);
			if(value < minOccur){
//				System.out.println("   R");
//				hashMap.remove(keyS);
			}else{
				hashMapAC.put(keyS, ll);
			}
			
		}
		// sort
		Map<String, List<CMLPeak>> hashMapS = new TreeMap<String, List<CMLPeak>>(hashMapAC);
		return hashMapS;
	}
	
	public static List<String> getListPath(MZData mzData){
		IReactionScheme reactionScheme = mzData.getListReactions();
		// looking for molecule 1_ . The top
	    List<String> listID = new ArrayList<String>();
	    if(reactionScheme == null)
	    	return listID;
		IMoleculeSet molSetI = ReactionSchemeManipulator.getAllMolecules(reactionScheme);
		IMoleculeSet molSet = builder.newMoleculeSet();
		for (IAtomContainer atomContainer : molSetI.molecules()) {
			molSet.addAtomContainer(atomContainer);
		}
		for (IAtomContainer atomContainer : molSetI.molecules()) {
	    	if(atomContainer.getID().contains("Loss"))
	    		molSet.removeAtomContainer(atomContainer);
	    }
		List<IAtomContainer> topMolecules = new ArrayList<IAtomContainer>();
		IReactionSet reactionS = ReactionSchemeManipulator.extractTopReactions(reactionScheme);
	    for(IReaction reaction:reactionS.reactions()){
	    	IAtomContainer molecule = reaction.getReactants().getAtomContainer(0);
	    	if(!listID.contains(molecule.getID())){
	    		listID.add(molecule.getID());
	    		topMolecules.add(molecule);
	    	}
	    }
	    for(IReaction reaction:reactionScheme.reactions()){
	    	IAtomContainer molecule = reaction.getReactants().getAtomContainer(0);
	    	if(!listID.contains(molecule.getID())){
	    		listID.add(molecule.getID());
	    		topMolecules.add(molecule);
	    	}
	    }
	    ArrayList<String> listPath = new ArrayList<String>();
	    for(IAtomContainer molecule : topMolecules){
		    ArrayList<String> listPathInt = new ArrayList<String>();
	    	if(((IMolecule)molecule).getProperty(CDKConstants.FORMULA) instanceof ArrayList){
	    		ArrayList<String> mfL = ((ArrayList<String>)((IMolecule)molecule).getProperty(CDKConstants.FORMULA));
	    		String path = "";
				for(int i = 0 ; i < mfL.size(); i++){
					path += mfL.get(i).replace(" ", "");
					if(i!=mfL.size()-1)
						path += "@";
	    		}
//				path += "||";
				listPathInt.add(path);
			}else if(((IMolecule)molecule).getProperty(CDKConstants.FORMULA) instanceof IMolecularFormula){
				listPathInt.add((MolecularFormulaManipulator.getString((IMolecularFormula)((IMolecule)molecule).getProperty(CDKConstants.FORMULA))).replace(" ", ""));
			}else if(((IMolecule)molecule).getProperty(CDKConstants.FORMULA) instanceof IMolecularFormulaSet){
				IMolecularFormulaSet formulaSet = (IMolecularFormulaSet)((IMolecule)molecule).getProperty(CDKConstants.FORMULA);
				String path = "";
				for(int i = 0; i < formulaSet.size(); i++){
					path += (MolecularFormulaManipulator.getString(formulaSet.getMolecularFormula(i))).replace(" ", "");
					if(i != formulaSet.size()-1)
						path += "@";
				}
//				path += "||";
				listPathInt.add(path);
			}
	    	
	    	for (IAtomContainer atomContainer : molSet.molecules()) {
		    	if(atomContainer.getID().contains("Loss"))
		    		continue;
		    	if(topMolecules.contains(atomContainer))
		    		continue;

		    	String path = getPath((IMolecule)molecule, (IMolecule)atomContainer,reactionScheme);
		    	if(path!=null && !path.endsWith("||"))
					listPathInt.add(path);

		    }
	    	Collections.sort(listPathInt);
	    	listPath.addAll(listPathInt);
	    	
	    }
	    return listPath;
	}

	public static List<String> getListPathNominal(MZData mzData){
		return getListPathNominal(mzData,"#0");
	}
	
	public static List<String> getListPathNominal(MZData mzData, String nom){
	    ArrayList<String> listPath = new ArrayList<String>();

		String str = "";
		CMLSpectrumList cmlSpectList = mzData.getListSpectra();
		List<CMLMetadata> metadataList = cmlSpectList.getMetadataListElements().get(0).getMetadataDescendants();
		int numGroups = 0;
		for(CMLMetadata metadata:metadataList){
			if(metadata.getDictRef().equals(MZDataConstants.NUM_GROUPS))
				numGroups = Integer.parseInt(metadata.getContent());
		}
		CMLElements<CMLSpectrum> specElem = cmlSpectList.getSpectrumElements();
		for(int i = 0 ; i < numGroups; i++){
			List<CMLSpectrum> newCMLSpeList = new ArrayList<CMLSpectrum>();
			String scanIDP = "";
			for(CMLSpectrum spectrum:specElem){
				List<CMLMetadata> ml = spectrum.getMetadataListElements().get(0).getMetadataDescendants();
				boolean validS = false;
				for(CMLMetadata metadata:ml){
					if(metadata.getDictRef().equals(MZDataConstants.GROUP_PEAK_MSN) && metadata.getContent().equals((new Integer(i+1)).toString())){
						newCMLSpeList.add(spectrum);
						validS = true;
					}
				}
				if(validS){
					for(CMLMetadata metadata:ml){
						if(metadata.getDictRef().equals(MZDataConstants.MS_LEVEL) && metadata.getContent().equals("1")){
							scanIDP = spectrum.getId();
						}
					}
				}
			}
			/////////////////////////////////////////////////////////
		    ParentIon ionFinal = extractParentIonN(newCMLSpeList,scanIDP,mzData);
		    List<String> listPathInt = new ArrayList<String>();
	    	
		    DecimalFormat Currency = new DecimalFormat(nom);
		    if(ionFinal == null)
		    	return null;
	        String formated_1 = Currency.format(ionFinal.getMass());
			String nmPath1 = formated_1;
			listPathInt.add(nmPath1);
		    
		    for(int f2 = 0 ; f2 < ionFinal.getFragments().size(); f2++){
	        	ParentIon ion2 = ionFinal.getFragments().get(f2);
		        String formated_2 = Currency.format(ion2.getMass());
	        	String nmPath2 = formated_1+"||"+formated_2;
        		listPathInt.add(nmPath2);
//        		System.out.println("nmPath2: "+nmPath2);
	        	
        		for(int f3 = 0 ; f3 < ion2.getFragments().size(); f3++){
	            	ParentIon ion3 = ion2.getFragments().get(f3);
    		        String formated_3 = Currency.format(ion3.getMass());
	            	String nmPath3 = formated_1+"||"+formated_2+"||"+formated_3;
	            	listPathInt.add(nmPath3);
//	        		System.out.println("nmPath3: "+nmPath3);
		        	for(int f4 = 0 ; f4 < ion3.getFragments().size(); f4++){
	                	ParentIon ion4 = ion3.getFragments().get(f4);
	    		        String formated_4 = Currency.format(ion4.getMass());
	                	String nmPath4 = formated_1+"||"+formated_2+"||"+formated_3+"||"+formated_4;
	                	listPathInt.add(nmPath4);
	    	        	
	                	for(int f5 = 0 ; f5 < ion4.getFragments().size(); f5++){
	                    	ParentIon ion5 = ion4.getFragments().get(f5);
    	    		        String formated_5 = Currency.format(ion5.getMass());
	                    	String nmPath5 = formated_1+"||"+formated_2+"||"+formated_3+"||"+formated_4+"||"+formated_5;
	                    	listPathInt.add(nmPath5);
	        	        	
	                    	for(int f6 = 0 ; f6 < ion5.getFragments().size(); f6++){
	                        	ParentIon ion6 = ion5.getFragments().get(f6);
	                        	String formated_6 = Currency.format(ion6.getMass());
        	    				String nmPath6 = formated_1+"||"+formated_2+"||"+formated_3+"||"+formated_4+"||"+formated_5+"||"+formated_6;
        	    				listPathInt.add(nmPath6);
        	    	        	for(int f7 = 0 ; f7 < ion6.getFragments().size(); f7++){
	                            	ParentIon ion7 = ion6.getFragments().get(f7);
	                            	String formated_7 = Currency.format(ion7.getMass());
	            	    		    String nmPath7 = formated_1+"||"+formated_2+"||"+formated_3+"||"+formated_4+"||"+formated_5+"||"+formated_6+"||"+formated_7;
	            	    		    listPathInt.add(nmPath7);
	            		        	
	                        	}
	                    	}
	                	}
	            	}
	        	}
		    }
	    	Collections.sort(listPathInt);
	    	listPath.addAll(listPathInt);
		}
        return listPath;
	}
	public static List<String> getListPathNominal_loss(MZData mzData){
		return getListPathNominal(mzData,"#0");
	}
	
	public static List<String> getListPathNominal_loss(MZData mzData, String nom){
	    ArrayList<String> listPath = new ArrayList<String>();

		String str = "";
		CMLSpectrumList cmlSpectList = mzData.getListSpectra();
		List<CMLMetadata> metadataList = cmlSpectList.getMetadataListElements().get(0).getMetadataDescendants();
		int numGroups = 0;
		for(CMLMetadata metadata:metadataList){
			if(metadata.getDictRef().equals(MZDataConstants.NUM_GROUPS))
				numGroups = Integer.parseInt(metadata.getContent());
		}
		CMLElements<CMLSpectrum> specElem = cmlSpectList.getSpectrumElements();
		for(int i = 0 ; i < numGroups; i++){
			List<CMLSpectrum> newCMLSpeList = new ArrayList<CMLSpectrum>();
			String scanIDP = "";
			for(CMLSpectrum spectrum:specElem){
				List<CMLMetadata> ml = spectrum.getMetadataListElements().get(0).getMetadataDescendants();
				boolean validS = false;
				for(CMLMetadata metadata:ml){
					if(metadata.getDictRef().equals(MZDataConstants.GROUP_PEAK_MSN) && metadata.getContent().equals((new Integer(i+1)).toString())){
						newCMLSpeList.add(spectrum);
						validS = true;
					}
				}
				if(validS){
					for(CMLMetadata metadata:ml){
						if(metadata.getDictRef().equals(MZDataConstants.MS_LEVEL) && metadata.getContent().equals("1")){
							scanIDP = spectrum.getId();
						}
					}
				}
			}
			/////////////////////////////////////////////////////////
		    ParentIon ionFinal = extractParentIonN(newCMLSpeList,scanIDP,mzData);
		    List<String> listPathInt = new ArrayList<String>();
	    	
		    DecimalFormat Currency = new DecimalFormat(nom);
		    if(ionFinal == null)
		    	return null;
//	        String formated_1 = Currency.format(ionFinal.getMass());
//			String nmPath1 = formated_1;
//			listPathInt.add(nmPath1);
		    
		    for(int f2 = 0 ; f2 < ionFinal.getFragments().size(); f2++){
	        	ParentIon ion2 = ionFinal.getFragments().get(f2);
		        String formated_2 = Currency.format(ionFinal.getMass()-ion2.getMass());
//	        	String nmPath2 = formated_1+"||"+formated_2;
	        	String nmPath2 = formated_2;
        		listPathInt.add(nmPath2);
//        		System.out.println("nmPath2: "+nmPath2);
	        	
        		for(int f3 = 0 ; f3 < ion2.getFragments().size(); f3++){
	            	ParentIon ion3 = ion2.getFragments().get(f3);
    		        String formated_3 = Currency.format(ion2.getMass()-ion3.getMass());
//	            	String nmPath3 = formated_1+"||"+formated_2+"||"+formated_3;
	            	String nmPath3 = formated_2+"||"+formated_3;
	            	listPathInt.add(nmPath3);
//	        		System.out.println("nmPath3: "+nmPath3);
		        	for(int f4 = 0 ; f4 < ion3.getFragments().size(); f4++){
	                	ParentIon ion4 = ion3.getFragments().get(f4);
	    		        String formated_4 = Currency.format(ion3.getMass()-ion4.getMass());
//	                	String nmPath4 = formated_1+"||"+formated_2+"||"+formated_3+"||"+formated_4;
	                	String nmPath4 = formated_2+"||"+formated_3+"||"+formated_4;
	                	listPathInt.add(nmPath4);
	    	        	
	                	for(int f5 = 0 ; f5 < ion4.getFragments().size(); f5++){
	                    	ParentIon ion5 = ion4.getFragments().get(f5);
    	    		        String formated_5 = Currency.format(ion4.getMass()-ion5.getMass());
//	                    	String nmPath5 = formated_1+"||"+formated_2+"||"+formated_3+"||"+formated_4+"||"+formated_5;
	                    	String nmPath5 = formated_2+"||"+formated_3+"||"+formated_4+"||"+formated_5;
	                    	listPathInt.add(nmPath5);
	        	        	
	                    	for(int f6 = 0 ; f6 < ion5.getFragments().size(); f6++){
	                        	ParentIon ion6 = ion5.getFragments().get(f6);
	                        	String formated_6 = Currency.format(ion5.getMass()-ion6.getMass());
//        	    				String nmPath6 = formated_1+"||"+formated_2+"||"+formated_3+"||"+formated_4+"||"+formated_5+"||"+formated_6;
        	    				String nmPath6 = formated_2+"||"+formated_3+"||"+formated_4+"||"+formated_5+"||"+formated_6;
        	    				listPathInt.add(nmPath6);
        	    	        	for(int f7 = 0 ; f7 < ion6.getFragments().size(); f7++){
	                            	ParentIon ion7 = ion6.getFragments().get(f7);
	                            	String formated_7 = Currency.format(ion6.getMass()-ion7.getMass());
//	            	    		    String nmPath7 = formated_1+"||"+formated_2+"||"+formated_3+"||"+formated_4+"||"+formated_5+"||"+formated_6+"||"+formated_7;
	            	    		    String nmPath7 = formated_2+"||"+formated_3+"||"+formated_4+"||"+formated_5+"||"+formated_6+"||"+formated_7;
	            	    		    listPathInt.add(nmPath7);
	            		        	
	                        	}
	                    	}
	                	}
	            	}
	        	}
		    }
	    	Collections.sort(listPathInt);
	    	listPath.addAll(listPathInt);
		}
        return listPath;
	}
	
	public static List<String> getListPathLoss(MZData mzData){
		IReactionScheme reactionScheme = mzData.getListReactions();
		IMoleculeSet molSetI = ReactionSchemeManipulator.getAllMolecules(reactionScheme);
		IMoleculeSet molSet = builder.newMoleculeSet();
		for (IAtomContainer atomContainer : molSetI.molecules()) {
			molSet.addAtomContainer(atomContainer);
		}
		for (IAtomContainer atomContainer : molSetI.molecules()) {
	    	String charge = (String) atomContainer.getProperty("cdk:partialCharge");
	    	if(charge.equals("0.0")){
	    		molSet.removeAtomContainer(atomContainer);
	    	}
	    }
		// looking for molecule 1_ . The top
	    List<String> listID = new ArrayList<String>();
		List<IAtomContainer> topMolecules = new ArrayList<IAtomContainer>();
		IReactionSet reactionS = ReactionSchemeManipulator.extractTopReactions(reactionScheme);
	    for(IReaction reaction:reactionS.reactions()){
	    	IAtomContainer molecule = reaction.getReactants().getAtomContainer(0);
	    	if(!listID.contains(molecule.getID())){
	    		listID.add(molecule.getID());
	    		topMolecules.add(molecule);
	    	}
	    }
	    for(IReaction reaction:reactionScheme.reactions()){
	    	IAtomContainer molecule = reaction.getReactants().getAtomContainer(0);
	    	if(!listID.contains(molecule.getID())){
	    		listID.add(molecule.getID());
	    		topMolecules.add(molecule);
	    	}
	    }
	    List<String> listPath = new ArrayList<String>();
	    for(IAtomContainer molecule : topMolecules){

	    	for (IAtomContainer atomContainer : molSet.molecules()) {
	    		
				String path = getPathLoss((IMolecule)molecule, (IMolecule)atomContainer,reactionScheme);
				if(path!=null)
					listPath.add(path);
		    }
	    	
	    }
	    return listPath;
	}

	public static List<String> getListPathVLoss(MZData mzData){
		List<String> listPath = MZDataManipulator.getListPathLoss(mzData);
		List<String> fingerList = new ArrayList<String>();
		for(String pathParent:listPath){
			List<String> listNei = new ArrayList<String>();
			String[] fm = pathParent.split("\\|\\|");
			
			for(String pathKid:listPath){

				String parent = "";
				if(!pathKid.contains("||") || pathKid.equals("||"))
					continue;
				else 
					parent = pathKid.substring(0, pathKid.lastIndexOf("||"));
				
				if(parent.equals(pathParent)){
					String[] fmK = pathKid.split("\\|\\|");
//					System.out.print("fmK: "+pathKid+" -> ");
					listNei.add(fmK[fmK.length-1]);
//					System.out.println(fmK[fmK.length-1]);
				}
					
			}
			if(listNei.size() < 2)
				continue;
			Collections.sort(listNei);
//			System.out.println(listNei.size());
			String sG = "["+listNei.get(0);
			for(int j = 1 ; j < listNei.size(); j++){
				if(j == listNei.size()-1){
					sG = sG+"||"+listNei.get(j)+"]";
				}else{
					sG += "||"+listNei.get(j);
				}
        	}
			fingerList.add(fm[fm.length-1]+"@"+sG);
		}
		return fingerList;
	}
	public static List<String> getListPathV(MZData mzData){
		List<String> listPath = MZDataManipulator.getListPath(mzData);
		List<String> fingerList = new ArrayList<String>();
		for(String pathParent:listPath){
			List<String> listNei = new ArrayList<String>();
			String[] fm = pathParent.split("\\|\\|");
			
			for(String pathKid:listPath){

				String parent = "";
				if(!pathKid.contains("||"))
					continue;
				else 
					parent = pathKid.substring(0, pathKid.lastIndexOf("||"));
				
				if(parent.equals(pathParent)){
					String[] fmK = pathKid.split("\\|\\|");
					listNei.add(fmK[fmK.length-1]);
//					System.out.println("fmK: "+pathParent+" -> "+fmK[fmK.length-1]);
				}
					
			}
			if(listNei.size() < 2)
				continue;
			Collections.sort(listNei);
//			System.out.println(listNei.size());
			String sG = "["+listNei.get(0);
			for(int j = 1 ; j < listNei.size(); j++){
				if(j == listNei.size()-1){
					sG = sG+"||"+listNei.get(j)+"]";
				}else{
					sG += "||"+listNei.get(j);
				}
        	}
			fingerList.add(fm[fm.length-1]+"@"+sG);
		}
		return fingerList;
	}
	private static String getPath(IMolecule origenMol, IMolecule finalMol, IReactionScheme reactionScheme) {
		ArrayList<IMoleculeSet> molSet = ReactionSchemeManipulator.getMoleculeSet(origenMol, finalMol, reactionScheme);
		if(molSet.size() == 0)
			return null;
		String path = "";
		int numM = molSet.get(0).getAtomContainerCount();
		int countM = 1;
		for(IAtomContainer molecule : molSet.get(0).molecules()){
			if(((IMolecule)molecule).getProperty(CDKConstants.FORMULA) instanceof ArrayList){
				ArrayList<String> mfL = ((ArrayList<String>)((IMolecule)molecule).getProperty(CDKConstants.FORMULA));
				for(int i = 0 ; i < mfL.size(); i++){
					path += mfL.get(i).replace(" ", "");
					if(i!=mfL.size()-1)
						path += "@";
	    		}
				if(numM != countM){
					path += "||";
					countM++;
				}else
					countM++;
			}else if(((IMolecule)molecule).getProperty(CDKConstants.FORMULA) instanceof IMolecularFormula){
				path += (MolecularFormulaManipulator.getString((IMolecularFormula)((IMolecule)molecule).getProperty(CDKConstants.FORMULA))).replace(" ", "");
				if(numM != countM){
					path += "||";
					countM++;
				}else
					countM++;
			}else if(((IMolecule)molecule).getProperty(CDKConstants.FORMULA) instanceof IMolecularFormulaSet){
				IMolecularFormulaSet formulaSet = (IMolecularFormulaSet)((IMolecule)molecule).getProperty(CDKConstants.FORMULA);
				
				for(int i = 0; i < formulaSet.size(); i++){
					path += (MolecularFormulaManipulator.getString(formulaSet.getMolecularFormula(i))).replace(" ", "");
					if(i != formulaSet.size()-1)
						path += "@";
				}
				if(numM != countM){
					path += "||";
					countM++;
				}else
					countM++;
			}
		}
		
		return path;
	}

	private static String getPathLoss(IMolecule origenMol, IMolecule finalMol, IReactionScheme reactionScheme) {
		ArrayList<IMoleculeSet> molSet = ReactionSchemeManipulator.getMoleculeSet(origenMol, finalMol, reactionScheme);
		IReactionSet reactions = ReactionSchemeManipulator.getAllReactions(reactionScheme);
		if(molSet.size() == 0)
			return null;
		String path = "";
		int numM = molSet.get(0).getMoleculeCount();
		int countMM = 1;
		int countM = 0;
		for(IAtomContainer molecule : molSet.get(0).molecules()){
			if(countM == 0){ // the first is always the top ion
				countM ++;
				continue;
			}
			
			// find the lost.
			boolean flag = false;
			for(IReaction reaction:reactions.reactions()){
				IMoleculeSet molecules = ReactionManipulator.getAllProducts(reaction);
				for(IAtomContainer mo:molecules.molecules()){
					if(mo.equals(molecule)){
						flag = true;
						break;
					}
				}
				if(flag){
					for(IAtomContainer mo:molecules.molecules()){
						if(!mo.equals(molecule)){
							molecule = mo;
							break;
						}
					}
					break;
				}
			}
			if(((IMolecule)molecule).getProperty(CDKConstants.FORMULA) instanceof ArrayList){
				ArrayList<String> mfL = ((ArrayList<String>)((IMolecule)molecule).getProperty(CDKConstants.FORMULA));
				for(int i = 0 ; i < mfL.size(); i++){
					path += mfL.get(i).replace(" ", "");
					if(i!=mfL.size()-1)
						path += "@";
	    		}
				if(numM != countMM){
					path += "||";
					countMM++;
				}else
					countMM++;

			}else if(((IMolecule)molecule).getProperty(CDKConstants.FORMULA) instanceof IMolecularFormula){
				path += (MolecularFormulaManipulator.getString((IMolecularFormula)((IMolecule)molecule).getProperty(CDKConstants.FORMULA))).replace(" ", "");
				if(numM != countM){
					path += "||";
					countM++;
				}else
					countM++;
			}else if(((IMolecule)molecule).getProperty(CDKConstants.FORMULA) instanceof IMolecularFormulaSet){
				IMolecularFormulaSet formulaSet = (IMolecularFormulaSet)((IMolecule)molecule).getProperty(CDKConstants.FORMULA);
				
				for(int i = 0; i < formulaSet.size(); i++){
					path += (MolecularFormulaManipulator.getString(formulaSet.getMolecularFormula(i))).replace(" ", "");
					if(i != formulaSet.size()-1)
						path += "@";
				}
				if(numM != countMM){
					path += "||";
					countMM++;
				}else
					countMM++;
			}
		}
		return path.substring(0, path.length()-2);
	}
	public static Map<String, List<IAtomContainer>> getHashMapPath(List<MZData> mzDataList){
		HashMap<String,List<IAtomContainer>> hashPath = new HashMap<String,List<IAtomContainer>>();
		
		for(MZData mzData:mzDataList){
		
			IReactionScheme reactionScheme = mzData.getListReactions();
			
			if(reactionScheme == null)
				continue;
			
			IMoleculeSet molSet = ReactionSchemeManipulator.getAllMolecules(reactionScheme);
			for (IAtomContainer atomContainer : molSet.molecules()) {
		    	if(atomContainer.getID().contains("_Loss"))
		    		molSet.removeAtomContainer(atomContainer);
		    }
			for (IAtomContainer atomContainer : molSet.molecules()) {
				atomContainer.setProperty("mzData", mzData.getID());
		    }
			// looking for molecule 1_ . The top
		    List<String> listID = new ArrayList<String>();
			List<IAtomContainer> topMolecules = new ArrayList<IAtomContainer>();
			IReactionSet reactionS = ReactionSchemeManipulator.extractTopReactions(reactionScheme);
		    for(IReaction reaction:reactionS.reactions()){
		    	IAtomContainer molecule = reaction.getReactants().getAtomContainer(0);
		    	if(!listID.contains(molecule.getID())){
		    		listID.add(molecule.getID());
		    		topMolecules.add(molecule);
		    	}
		    }
		    for(IReaction reaction:reactionScheme.reactions()){
		    	IAtomContainer molecule = reaction.getReactants().getAtomContainer(0);
		    	if(!listID.contains(molecule.getID())){
		    		listID.add(molecule.getID());
		    		topMolecules.add(molecule);
		    	}
		    }
		    for(IAtomContainer molecule : topMolecules){
		    	
		    	String finalPath = efExtractor(molecule);
		    	
		    	if(!hashPath.containsKey(finalPath)){
		    		List<IAtomContainer> listIA = new ArrayList<IAtomContainer>();
		    		listIA.add(molecule);
		    		hashPath.put(finalPath, listIA);
		    	}else{
		    		List<IAtomContainer> listIA = hashPath.get(finalPath);
		    		listIA.add(molecule);
		    		hashPath.put(finalPath, listIA);
		    	}
		    	
		    	for (IAtomContainer atomContainer : molSet.molecules()) {
			    	if(atomContainer.getID().contains("Loss"))
			    		continue;
			    	if(atomContainer.equals(molecule))
			    		continue;
			    	
					String path = getPath((IMolecule)molecule, (IMolecule)atomContainer,reactionScheme);
					if(path!=null){
	//					System.out.println(atomContainer.getID()+" path: "+path);
						if(!hashPath.containsKey(path)){
				    		List<IAtomContainer> listIA = new ArrayList<IAtomContainer>();
				    		listIA.add(atomContainer);
				    		hashPath.put(path, listIA);
				    	}else{
				    		List<IAtomContainer> listIA = hashPath.get(path);
				    		listIA.add(atomContainer);
				    		hashPath.put(path, listIA);
				    	}
					}
			    }
		    }
//		    break;
		}
		Map<String, List<IAtomContainer>> hashMapS = new TreeMap<String, List<IAtomContainer>>(hashPath);
	    return hashMapS;
	}

	public static Map<String, List<CMLPeak>> getHashMapPathN(MZData mzData){
		HashMap<String,List<CMLPeak>> hashPath = new HashMap<String,List<CMLPeak>>();
		HashMap<CMLPeak,String> listPathInt = new HashMap<CMLPeak, String>();
		
		/////////////////////////////////////////////////////////////////////////////

		String str = "";
		CMLSpectrumList cmlSpectList = mzData.getListSpectra();
		List<CMLMetadata> metadataList = cmlSpectList.getMetadataListElements().get(0).getMetadataDescendants();
		int numGroups = 0;
		for(CMLMetadata metadata:metadataList){
			if(metadata.getDictRef().equals(MZDataConstants.NUM_GROUPS))
				numGroups = Integer.parseInt(metadata.getContent());
		}
		CMLElements<CMLSpectrum> specElem = cmlSpectList.getSpectrumElements();
		for(int i = 0 ; i < numGroups; i++){
			List<CMLSpectrum> newCMLSpeList = new ArrayList<CMLSpectrum>();
			String scanIDP = "";
			for(CMLSpectrum spectrum:specElem){
				List<CMLMetadata> ml = spectrum.getMetadataListElements().get(0).getMetadataDescendants();
				boolean validS = false;
				for(CMLMetadata metadata:ml){
					if(metadata.getDictRef().equals(MZDataConstants.GROUP_PEAK_MSN) && metadata.getContent().equals((new Integer(i+1)).toString())){
						newCMLSpeList.add(spectrum);
						validS = true;
					}
				}
				if(validS){
					for(CMLMetadata metadata:ml){
						if(metadata.getDictRef().equals(MZDataConstants.MS_LEVEL) && metadata.getContent().equals("1")){
							scanIDP = spectrum.getId();
						}
					}
				}
			}
			/////////////////////////////////////////////////////////
		    ParentIon ionFinal = extractParentIonN(newCMLSpeList,scanIDP,mzData);
		    
		    DecimalFormat Currency = new DecimalFormat("#0");
	        String formated_1 = Currency.format(ionFinal.getMass());
			String nmPath1 = formated_1;
			listPathInt.put(ionFinal.getCMLPeak(),nmPath1);
		    
		    for(int f2 = 0 ; f2 < ionFinal.getFragments().size(); f2++){
	        	ParentIon ion2 = ionFinal.getFragments().get(f2);
		        String formated_2 = Currency.format(ion2.getMass());
	        	String nmPath2 = formated_1+"||"+formated_2;
        		listPathInt.put(ion2.getCMLPeak(),nmPath2);
	        	
        		for(int f3 = 0 ; f3 < ion2.getFragments().size(); f3++){
	            	ParentIon ion3 = ion2.getFragments().get(f3);
    		        String formated_3 = Currency.format(ion3.getMass());
	            	String nmPath3 = formated_1+"||"+formated_2+"||"+formated_3;
	            	listPathInt.put(ion3.getCMLPeak(),nmPath3);
		        	for(int f4 = 0 ; f4 < ion3.getFragments().size(); f4++){
	                	ParentIon ion4 = ion3.getFragments().get(f4);
	    		        String formated_4 = Currency.format(ion4.getMass());
	                	String nmPath4 = formated_1+"||"+formated_2+"||"+formated_3+"||"+formated_4;
	                	listPathInt.put(ion4.getCMLPeak(),nmPath4);
	    	        	
	                	for(int f5 = 0 ; f5 < ion4.getFragments().size(); f5++){
	                    	ParentIon ion5 = ion4.getFragments().get(f5);
    	    		        String formated_5 = Currency.format(ion5.getMass());
	                    	String nmPath5 = formated_1+"||"+formated_2+"||"+formated_3+"||"+formated_4+"||"+formated_5;
	                    	listPathInt.put(ion5.getCMLPeak(),nmPath5);
	        	        	
	                    	for(int f6 = 0 ; f6 < ion5.getFragments().size(); f6++){
	                        	ParentIon ion6 = ion5.getFragments().get(f6);
	                        	String formated_6 = Currency.format(ion6.getMass());
        	    				String nmPath6 = formated_1+"||"+formated_2+"||"+formated_3+"||"+formated_4+"||"+formated_5+"||"+formated_6;
        	    				listPathInt.put(ion6.getCMLPeak(),nmPath6);
        	    	        	for(int f7 = 0 ; f7 < ion6.getFragments().size(); f7++){
	                            	ParentIon ion7 = ion6.getFragments().get(f7);
	                            	String formated_7 = Currency.format(ion7.getMass());
	            	    		    String nmPath7 = formated_1+"||"+formated_2+"||"+formated_3+"||"+formated_4+"||"+formated_5+"||"+formated_6+"||"+formated_7;
	            	    		    listPathInt.put(ion7.getCMLPeak(),nmPath7);
	            		        	
	                        	}
	                    	}
	                	}
	            	}
	        	}
		    }
		}
        Set<Entry<CMLPeak,String>> set = listPathInt.entrySet();
		for( Iterator<Entry<CMLPeak,String>> ita = set.iterator(); ita.hasNext();){
			Entry<CMLPeak,String> keyl = ita.next();
        	CMLPeak cmlPeak = keyl.getKey();
        	String finalPath = keyl.getValue();

        	if(!hashPath.containsKey(finalPath)){
        		List<CMLPeak> listIA = new ArrayList<CMLPeak>();
        		listIA.add(cmlPeak);
        		hashPath.put(finalPath, listIA);
        	}else{
        		List<CMLPeak> listIA = hashPath.get(finalPath);
        		listIA.add(cmlPeak);
        		hashPath.put(finalPath, listIA);
        	}
		}
	    
	    	
		Map<String, List<CMLPeak>> hashMapS = new TreeMap<String, List<CMLPeak>>(hashPath);
	    return hashMapS;
	}
	private static String efExtractor(IAtomContainer molecule) {

    	String finalPath = "";
	    if(((IMolecule)molecule).getProperty(CDKConstants.FORMULA) instanceof ArrayList){
    		ArrayList<String> mfL = ((ArrayList<String>)((IMolecule)molecule).getProperty(CDKConstants.FORMULA));
    		String path = "";
			for(int i = 0 ; i < mfL.size(); i++){
				path += mfL.get(i).replace(" ", "");
				if(i!=mfL.size()-1)
					path += "@";
    		}
//			path += "||";
			finalPath = path;
		}else if(((IMolecule)molecule).getProperty(CDKConstants.FORMULA) instanceof IMolecularFormula){
			finalPath = (MolecularFormulaManipulator.getString((IMolecularFormula)((IMolecule)molecule).getProperty(CDKConstants.FORMULA))).replace(" ", "");
		}else if(((IMolecule)molecule).getProperty(CDKConstants.FORMULA) instanceof IMolecularFormulaSet){
			IMolecularFormulaSet formulaSet = (IMolecularFormulaSet)((IMolecule)molecule).getProperty(CDKConstants.FORMULA);
			String path = "";
			for(int i = 0; i < formulaSet.size(); i++){
				path += (MolecularFormulaManipulator.getString(formulaSet.getMolecularFormula(i))).replace(" ", "");
				if(i != formulaSet.size()-1)
					path += "@";
			}
//			path += "||";
			finalPath = path;
		}
		return finalPath;
	}
	public static HashMap<String,List<IAtomContainer>> getHashMapPathLoss(IReactionScheme reactionScheme){

		IMoleculeSet molSet = ReactionSchemeManipulator.getAllMolecules(reactionScheme);
		for (IAtomContainer atomContainer : molSet.molecules()) {
	    	if(atomContainer.getID().contains("_Loss"))
	    		molSet.removeAtomContainer(atomContainer);
	    }
		// looking for molecule 1_ . The top
	    List<String> listID = new ArrayList<String>();
		List<IAtomContainer> topMolecules = new ArrayList<IAtomContainer>();
		IReactionSet reactionS = ReactionSchemeManipulator.extractTopReactions(reactionScheme);
	    for(IReaction reaction:reactionS.reactions()){
	    	IAtomContainer molecule = reaction.getReactants().getAtomContainer(0);
	    	if(!listID.contains(molecule.getID())){
	    		listID.add(molecule.getID());
	    		topMolecules.add(molecule);
	    	}
	    }
	    for(IReaction reaction:reactionScheme.reactions()){
	    	IAtomContainer molecule = reaction.getReactants().getAtomContainer(0);
	    	if(!listID.contains(molecule.getID())){
	    		listID.add(molecule.getID());
	    		topMolecules.add(molecule);
	    	}
	    }
	    HashMap<String,List<IAtomContainer>> hashPath = new HashMap<String,List<IAtomContainer>>();
	    for(IAtomContainer molecule : topMolecules){
	    	
	    	for (IAtomContainer atomContainer : molSet.molecules()) {
		    	if(atomContainer.getID().contains("Loss"))
		    		continue;
		    	if(atomContainer.equals(molecule))
		    		continue;
		    	
				String path = getPathLoss((IMolecule)molecule, (IMolecule)atomContainer,reactionScheme);
				if(path!=null)
					if(!hashPath.containsKey(path)){
			    		List<IAtomContainer> listIA = new ArrayList<IAtomContainer>();
			    		listIA.add(atomContainer);
			    		hashPath.put(path, listIA);
			    	}else{
			    		List<IAtomContainer> listIA = hashPath.get(path);
			    		listIA.add(atomContainer);
			    		hashPath.put(path, listIA);
			    	}
		    }
	    	
	    }
	    return hashPath;
	}


//	public static MZData getMZData(List<String> listPath, MZData mzData_Or) {
//		int charge = 0;
//		
//		Polarity polarity = getPolarity(mzData_Or);
//		
//		if(polarity.equals(Polarity.positive))
//			charge = 1;
//		else if(polarity.equals(Polarity.negative))
//			charge = -1;
//		
//		int countSp = 1;
//		int countSi = 1;
//		
//		MZData mzData = new MZData();
//		
//        CMLSpectrumList cmlSpecList = new CMLSpectrumList();
//		CMLMetadataList metadataListL = new CMLMetadataList();
//		/////////////////////////////////
//		CMLMetadata metadataL = new CMLMetadata();
//		metadataL.setDictRef(MZDataConstants.NUM_GROUPS);
//		metadataL.setContent(Integer.toString(1));
//		metadataListL.appendChild(metadataL);
//		/////////////////////////////////
//		cmlSpecList.addMetadataList(metadataListL);
//        HashMap<String, String> hashS = new HashMap<String, String>();
//        
//		for(String sOO : listPath){
//			
//			boolean init = false;
//			
//			if(sOO.lastIndexOf("||") == -1){ // full scan
//				
//				CMLSpectrum cmlSpect = new CMLSpectrum();
//				CMLPeakList cmlPeaks = new CMLPeakList();
//    			CMLMetadataList metadataList = new CMLMetadataList();
//    			cmlSpect.addMetadataList(metadataList);
//				CMLConditionList conditionList = new CMLConditionList();
//				cmlSpect.addConditionList(conditionList);
//    			cmlSpect.setId("spec"+Integer.toString(countSp++));
//    			cmlSpect.setAttribute("type", "MS");
//    			hashS.put(sOO, Integer.toString(countSp-1));
//    			/////////////////////////////////
//    			CMLMetadata metadata = new CMLMetadata();
//    			metadata.setDictRef(MZDataConstants.MS_LEVEL);
//    			metadata.setContent(Integer.toString(1));
//    			metadataList.appendChild(metadata);
//    			/////////////////////////////////
//				metadata = new CMLMetadata();
//				metadata.setDictRef(MZDataConstants.GROUP_PEAK_MSN);
//				metadata.setContent(""+1);
//				metadataList.appendChild(metadata);
//    			/////////////////////////////////
//				metadata = new CMLMetadata();
//				metadata.setDictRef(MZDataConstants.SCAN_NUM);
//				metadata.setContent(Integer.toString(countSp-1));
//				metadataList.appendChild(metadata);
//				/////////////////////////////////
//				CMLScalar condition = new CMLScalar();		
//				condition.setDictRef(MZDataConstants.PRECURSOR_MZ);
//				String precurMZ = "0";
//				condition.setValue(precurMZ);
//				conditionList.appendChild(condition);
//				/////////////////////////////////
//				condition = new CMLScalar();		
//				condition.setDictRef(MZDataConstants.POLARITY);
//				condition.setValue(polarity.toString());
//				conditionList.appendChild(condition);
//				/////////////////////////////////
//    			
//				IMolecularFormula formula = MolecularFormulaManipulator.getMolecularFormula(sOO, builder);
//        		formula.setCharge(charge);
//        		double mass = MolecularFormulaManipulator.getTotalExactMass(formula);
//    			Double xValue = mass;
//            	Double yValue = 100.00;
//            	CMLPeak peak = new CMLPeak();
//    			peak.setId("peak"+Integer.toString(countSi++));
//    			peak.setXValue(xValue.toString());
//    			peak.setYValue(yValue.toString());
//    			cmlPeaks.addPeak(peak); 
//    	        cmlSpect.addPeakList(cmlPeaks);
//    			cmlSpecList.addSpectrum(cmlSpect);
//			}
//			CMLPeakList cmlPeaks = new CMLPeakList();
//			
//    		for(String sO : listPath){
//    			if(sO.lastIndexOf("||") == -1){ // full scan
//    				continue;
//    			}
//	        	String sP = sO.substring(0, sO.lastIndexOf("||"));
//	        	if(sP.equals(sOO)){
//
//        			hashS.put(sO, Integer.toString(countSp));
//	        		if(!init){ // initiate spectrum 
//
//	        	        CMLSpectrum cmlSpect = new CMLSpectrum();
//	        			CMLMetadataList metadataList = new CMLMetadataList();
//	        			cmlSpect.addMetadataList(metadataList);
//	    				CMLConditionList conditionList = new CMLConditionList();
//	    				cmlSpect.addConditionList(conditionList);
//	        			cmlSpect.setId("spec"+Integer.toString(countSp++));
//	        			cmlSpect.setAttribute("type", "MS");
//	        			/////////////////////////////////
//	        			CMLMetadata metadata = new CMLMetadata();
//	        			metadata.setDictRef(MZDataConstants.MS_LEVEL);
//	        			metadata.setContent(Integer.toString(sO.split("\\|\\|").length));
//	        			metadataList.appendChild(metadata);
//	        			/////////////////////////////////
//						metadata = new CMLMetadata();
//						metadata.setDictRef(MZDataConstants.GROUP_PEAK_MSN);
//						metadata.setContent(""+1);
//						metadataList.appendChild(metadata);
//	        			/////////////////////////////////
//	    				metadata = new CMLMetadata();
//	    				metadata.setDictRef(MZDataConstants.SCAN_NUM);
//	    				metadata.setContent(Integer.toString(countSp-1));
//	    				metadataList.appendChild(metadata);
//	        			/////////////////////////////////
//	    				CMLScalar condition = new CMLScalar();		
//	    				condition.setDictRef(MZDataConstants.PRECURSOR_MZ);
//	    				String[] ll = sOO.split("\\|\\|");
//		        		String sL = ll[ll.length-1];
//		        		IMolecularFormula formula = MolecularFormulaManipulator.getMolecularFormula(sL, builder);
//		        		formula.setCharge(charge);
//		        		double mass = MolecularFormulaManipulator.getTotalExactMass(formula);
//		    			DecimalFormat df = new DecimalFormat("#.###");
//		    			condition.setValue(df.format(mass));
//	    				conditionList.appendChild(condition);
//						/////////////////////////////////
//						condition = new CMLScalar();		
//						condition.setDictRef(MZDataConstants.PRECURSOR_SCAN);
//						condition.setValue((String)hashS.get(sOO));
//						conditionList.appendChild(condition);
//						/////////////////////////////////
//						condition = new CMLScalar();		
//						condition.setDictRef(MZDataConstants.POLARITY);
//						condition.setValue(polarity.toString());
//						conditionList.appendChild(condition);
//						/////////////////////////////////
//	        			
//	        	        cmlSpect.addPeakList(cmlPeaks);
//	        			cmlSpecList.addSpectrum(cmlSpect);
//	        			
//	        			init = true;
//	        		}
//	        		
////	        		System.out.println("Si:"+sO+">"+sP+">");
//	        		String[] ll = sO.split("\\|\\|");
//	        		String sL = ll[ll.length-1];
//	        		IMolecularFormula formula = MolecularFormulaManipulator.getMolecularFormula(sL, builder);
//	        		formula.setCharge(charge);
//	        		double mass = MolecularFormulaManipulator.getTotalExactMass(formula);
//	    			Double xValue = mass;
//	            	Double yValue = 100.00;
//	            	CMLPeak peak = new CMLPeak();
//	    			peak.setId("peak"+Integer.toString(countSi++));
//	    			peak.setXValue(xValue.toString());
//	    			peak.setYValue(yValue.toString());
//	    			cmlPeaks.addPeak(peak); 
//	        	}
//	        }
//		}
//       
//		
//		mzData.setListSpectra(cmlSpecList);
////		mzData.setListReactions(reactionScheme);
//		return mzData;
//	}
	
	private static Polarity getPolarity(MZData mzData) {
		CMLConditionList conditionList = mzData.getListSpectra().getSpectrumElements().get(0).getConditionListElements().get(0);
		Iterator conditionIterator = conditionList.getScalarElements().iterator();
		while (conditionIterator.hasNext()) {
			CMLScalar condition = (CMLScalar) conditionIterator.next();
			if(condition.getDictRef().equals(MZDataConstants.POLARITY) ){
				String polarityTMP = condition.getValue();
				if(polarityTMP.equals("positive"))
					return Polarity.positive;
				else if(polarityTMP.equals("negative"))
					return Polarity.negative;
			}
		}
		return Polarity.nothing;
	}

	public static MZData getMZData(Map<String, List<IAtomContainer>> hashMap, List<MZData> mzDataList) {
		int charge = 0; // seting charge molecules
		MZData mzData_old = mzDataList.get(0);
		
		String intrumentRef = "";
		String protocol = "";
		List<CMLMetadata> metadataListI = mzData_old.getListSpectra().getMetadataListElements().get(0).getMetadataDescendants();
		for(CMLMetadata metadata:metadataListI){
			if(metadata.getDictRef().equals(MZDataConstants.INSTRUMENT))
				intrumentRef = metadata.getContent();
			else if(metadata.getDictRef().equals(MZDataConstants.PROTOCOL))
				protocol = metadata.getContent();
		}
		Polarity polarity = getPolarity(mzData_old);
		
		if(polarity.equals(Polarity.positive))
			charge = 1;
		else if(polarity.equals(Polarity.negative))
			charge = -1;

		String activationMethod = getActivationMethod(mzData_old);
		String window = getWindow(mzData_old);
		
		int countSp = 1;
		int countSi = 1;
		
		MZData mzData = new MZData();

		ReactionScheme reactionSet = new ReactionScheme();
		HashMap<String, ReactionScheme> schemeList = new HashMap<String, ReactionScheme>();
		int countScheme = 1;
		int countReaction = 1;
		
		
        CMLSpectrumList cmlSpecList = new CMLSpectrumList();
		CMLMetadataList metadataListL = new CMLMetadataList();
		/////////////////////////////////
		CMLMetadata metadataL = new CMLMetadata();
		metadataL.setDictRef(MZDataConstants.NUM_GROUPS);
		metadataL.setContent(Integer.toString(1));
		metadataListL.appendChild(metadataL);
		/////////////////////////////////
		if(!protocol.equals("")){
			metadataL = new CMLMetadata();
			metadataL.setDictRef(MZDataConstants.PROTOCOL);
			metadataL.setContent(protocol);
			metadataListL.appendChild(metadataL);
		}
		/////////////////////////////////
		if(!intrumentRef.equals("")){
			metadataL = new CMLMetadata();
			metadataL.setDictRef(MZDataConstants.INSTRUMENT);
			metadataL.setContent(intrumentRef);
			metadataListL.appendChild(metadataL);
		}
		/////////////////////////////////
		cmlSpecList.addMetadataList(metadataListL);
        HashMap<String, String> hashS = new HashMap<String, String>();
        HashMap<String, IMolecule> hashAC = new HashMap<String, IMolecule>();
        HashMap<String, IMolecule> hashAC_L = new HashMap<String, IMolecule>();//loss
        
        Set<Entry<String, List<IAtomContainer>>> set = hashMap.entrySet();
        
        // List molecules
		for( Iterator<Entry<String, List<IAtomContainer>>> it = set.iterator(); it.hasNext();){
			Entry<String, List<IAtomContainer>> key = it.next();
			List<IAtomContainer> l_l = key.getValue();
			String sOO = key.getKey();
//			System.out.println(sOO +" "+l_l);
			
			boolean init = false;

			//reaction
			IMolecule reactant;
			IMolecularFormula formula;
			if(sOO.length() == 0)
				continue;
			if(sOO.lastIndexOf("||") == -1){ 
				formula = MolecularFormulaManipulator.getMolecularFormula(sOO, builder);
			}else{
				formula = MolecularFormulaManipulator.getMolecularFormula(sOO.substring(sOO.lastIndexOf("||"), sOO.length()), builder);
			}
				
    		if(!hashAC.containsKey(sOO)){
    			reactant = (IMolecule) l_l.get(0);
//				reactant.setID(l_l.get(0).getID());
//				formula.setCharge(charge);
				reactant.setProperty(CDKConstants.FORMULA, formula);
//				reactant.setProperty("cdk:partialCharge", Double.toString(charge));
				hashAC.put(sOO, reactant);
			}else{
				reactant = (IMolecule) hashAC.get(sOO);
			}
    		
			if(sOO.lastIndexOf("||") == -1){ // full scan
				
				// Spectrum 
				CMLSpectrum cmlSpect = new CMLSpectrum();
				CMLPeakList cmlPeaks = new CMLPeakList();
    			CMLMetadataList metadataList = new CMLMetadataList();
    			cmlSpect.addMetadataList(metadataList);
				CMLConditionList conditionList = new CMLConditionList();
				cmlSpect.addConditionList(conditionList);
    			cmlSpect.setAttribute("type", "MS");
    			cmlSpect.setId(Integer.toString(countSp++));
    			hashS.put(sOO, Integer.toString(countSp-1));
    			/////////////////////////////////
    			CMLMetadata metadata = new CMLMetadata();
    			metadata.setDictRef(MZDataConstants.MS_LEVEL);
    			metadata.setContent(Integer.toString(1));
    			metadataList.appendChild(metadata);
    			/////////////////////////////////
				metadata = new CMLMetadata();
				metadata.setDictRef(MZDataConstants.GROUP_PEAK_MSN);
				metadata.setContent(""+1);
				metadataList.appendChild(metadata);
    			/////////////////////////////////
				metadata = new CMLMetadata();
				metadata.setDictRef(MZDataConstants.SCAN_NUM);
				metadata.setContent(Integer.toString(countSp-1));
				metadataList.appendChild(metadata);
				/////////////////////////////////
				CMLScalar condition = new CMLScalar();		
				condition.setDictRef(MZDataConstants.PRECURSOR_MZ);
				String precurMZ = "0";
				condition.setValue(precurMZ);
				conditionList.appendChild(condition);
				/////////////////////////////////
				condition = new CMLScalar();		
				condition.setDictRef(MZDataConstants.POLARITY);
				condition.setValue(polarity.toString());
				conditionList.appendChild(condition);
				/////////////////////////////////
				condition = new CMLScalar();		
				condition.setDictRef(MZDataConstants.ACTIVATION_METHOD);
				condition.setValue(activationMethod);
				conditionList.appendChild(condition);
				/////////////////////////////////
				condition = new CMLScalar();		
				condition.setDictRef(MZDataConstants.WINDOW);
				condition.setValue(window);
				conditionList.appendChild(condition);
				/////////////////////////////////
				Double[] xyValue = getXYValueAverage(l_l,mzDataList);
    			Double xValue = xyValue[0];
            	Double yValue = xyValue[1];
            	CMLPeak peak = new CMLPeak();
    			peak.setId("peak"+Integer.toString(countSi++));
    			peak.setXValue(xValue.toString());
    			peak.setYValue(yValue.toString());
				peak.setMoleculeRefs(l_l.get(0).getID());
    			cmlPeaks.addPeak(peak); 
    	        cmlSpect.addPeakList(cmlPeaks);
    			cmlSpecList.addSpectrum(cmlSpect);
    			
    			// reaction
    			ReactionScheme scheme = new ReactionScheme();
    			schemeList.put(sOO, scheme);
    			scheme.setID("rs"+countScheme++);
    			reactionSet.add(scheme);

				reactant.setProperty(MZDataConstants.GROUP_PEAK_MSN, "1");
				reactant.setProperty(MZDataConstants.MASS, xValue.toString());
				reactant.setProperty(MZDataConstants.INTENSITIY, yValue.toString());
    		}
			
			CMLPeakList cmlPeaks = new CMLPeakList();
			
			
			Set<Entry<String, List<IAtomContainer>>> setl = hashMap.entrySet();
			for( Iterator<Entry<String, List<IAtomContainer>>> itl = setl.iterator(); itl.hasNext();){
				Entry<String, List<IAtomContainer>> keyl = itl.next();
				List<IAtomContainer> l_ll = keyl.getValue();
				String sO = keyl.getKey();
				if(sO.lastIndexOf("||") == -1){ // full scan
					// get subfragments. If not contained means it is the only one.
					boolean found = false;
					for( Iterator<Entry<String, List<IAtomContainer>>> itlIn = setl.iterator(); itlIn.hasNext();){
						Entry<String, List<IAtomContainer>> keylIn = itlIn.next();
						String sOIn = keylIn.getKey();
						if(sOIn.lastIndexOf("||") != -1){
							String sP = sOIn.substring(0, sOIn.lastIndexOf("||"));
			        		if(sP.equals(sO)){
			        			found = true;
			        			break; // contains fragmetns
			        		}
						}
					}
					if(!found){
		        		// reaction
		    			IReaction reaction = new Reaction();
		    	    	reaction.setID("react"+countReaction++);
		    	    	reaction.addReactant(reactant);
		    	    	if(containProducts(hashMap.entrySet(),sO)){
	    					ReactionScheme scheme1 = new ReactionScheme();
		        			schemeList.put(sO, scheme1);
		    				scheme1.setID("rs"+countScheme++);
		    	    		scheme1.addReaction(reaction);
		    	    		if(schemeList.containsKey(sOO))
		    	    			schemeList.get(sOO).add(scheme1);
		    	    	}else{
		    	    		if(schemeList.containsKey(sOO))
		    	    			schemeList.get(sOO).addReaction(reaction);
		    	    	}
					}
					continue;
    			}
    			
	        	String sP = sO.substring(0, sO.lastIndexOf("||"));
        		if(sP.equals(sOO)){

	        		
	        		if(!init){ // initiate spectrum 

	        	        CMLSpectrum cmlSpect = new CMLSpectrum();
	        			CMLMetadataList metadataList = new CMLMetadataList();
	        			cmlSpect.addMetadataList(metadataList);
	    				CMLConditionList conditionList = new CMLConditionList();
	    				cmlSpect.addConditionList(conditionList);
	        			cmlSpect.setAttribute("type", "MS");
	        			cmlSpect.setId(Integer.toString(countSp++));
	        			/////////////////////////////////
	        			CMLMetadata metadata = new CMLMetadata();
	        			metadata.setDictRef(MZDataConstants.MS_LEVEL);
	        			metadata.setContent(Integer.toString(sO.split("\\|\\|").length));
	        			metadataList.appendChild(metadata);
	        			/////////////////////////////////
						metadata = new CMLMetadata();
						metadata.setDictRef(MZDataConstants.GROUP_PEAK_MSN);
						metadata.setContent(""+1);
						metadataList.appendChild(metadata);
	        			/////////////////////////////////
	    				metadata = new CMLMetadata();
	    				metadata.setDictRef(MZDataConstants.SCAN_NUM);
	    				metadata.setContent(Integer.toString(countSp-1));
	    				metadataList.appendChild(metadata);
	        			/////////////////////////////////
	    				CMLScalar condition = new CMLScalar();		
	    				condition.setDictRef(MZDataConstants.PRECURSOR_MZ);
	    				String[] ll = sOO.split("\\|\\|");
		        		String sL = ll[ll.length-1];
		        		IMolecularFormula formulaP = MolecularFormulaManipulator.getMolecularFormula(sL, builder);
		        		formulaP.setCharge(charge);
		        		double mass = MolecularFormulaManipulator.getTotalExactMass(formulaP);
		    			DecimalFormat df = new DecimalFormat("#.###");
		    			condition.setValue(df.format(mass));
	    				conditionList.appendChild(condition);
						/////////////////////////////////
						condition = new CMLScalar();		
						condition.setDictRef(MZDataConstants.PRECURSOR_SCAN);
						condition.setValue((String)hashS.get(sOO));
						conditionList.appendChild(condition);
						/////////////////////////////////
						condition = new CMLScalar();		
						condition.setDictRef(MZDataConstants.POLARITY);
						condition.setValue(polarity.toString());
						conditionList.appendChild(condition);
						/////////////////////////////////
						condition = new CMLScalar();		
						condition.setDictRef(MZDataConstants.ACTIVATION_METHOD);
						condition.setValue(activationMethod);
						conditionList.appendChild(condition);
						/////////////////////////////////
						condition = new CMLScalar();		
						condition.setDictRef(MZDataConstants.WINDOW);
						condition.setValue(window);
						conditionList.appendChild(condition);
						/////////////////////////////////
	        			
	        	        cmlSpect.addPeakList(cmlPeaks);
	        			cmlSpecList.addSpectrum(cmlSpect);
	        			
	        			init = true;
	        		}
        			hashS.put(sO, Integer.toString(countSp-1));
	        		// spectrum 
	        		String[] ll = sO.split("\\|\\|");
	        		String sL = ll[ll.length-1];
	        		
	        		IMolecularFormula formulaP = MolecularFormulaManipulator.getMolecularFormula(sL.split("@")[0], builder);
	        		formulaP.setCharge(charge);

					Double[] xyValue = getXYValueAverage(l_ll,mzDataList);
	    			Double xValue = xyValue[0];
	            	Double yValue = xyValue[1];
	            	CMLPeak peak = new CMLPeak();
	    			peak.setId("peak"+Integer.toString(countSi++));
	    			peak.setXValue(xValue.toString());
	    			peak.setYValue(yValue.toString());
					peak.setMoleculeRefs(l_ll.get(0).getID());
	    			cmlPeaks.addPeak(peak);

	        		// reaction
	    			IReaction reaction = new Reaction();
	    	    	reaction.setID("react"+countReaction++);
	    	    	reaction.addReactant(reactant);

	    			IMolecule product;
	    			if(!hashAC.containsKey(sO)){
	    				product = (IMolecule) l_ll.get(0);
		    			product.setID(l_ll.get(0).getID());
		    			
		    			if(sL.split("@").length == 1)
			    			product.setProperty(CDKConstants.FORMULA, formulaP);
		    			else{
		    				IMolecularFormulaSet formulaSet = builder.newMolecularFormulaSet();
		    				for(int s = 0 ; s < sL.split("@").length; s++){
		    					IMolecularFormula formulaPS = MolecularFormulaManipulator.getMolecularFormula(sL.split("@")[s], builder);
		    					formulaPS.setCharge(charge);
		    	        		formulaSet.addMolecularFormula(formulaPS);
		    				}
		    				// TODO how to solve that
//			    			product.setProperty(CDKConstants.FORMULA, formulaSet);
		    			}
		    			product.setProperty("cdk:partialCharge", Double.toString(charge));
		    			product.setProperty(MZDataConstants.GROUP_PEAK_MSN, "1");
		    			product.setProperty(MZDataConstants.MASS, xValue.toString());
		    			product.setProperty(MZDataConstants.INTENSITIY, yValue.toString());
		    			product.setProperty(MZDataConstants.PRECURSOR_ID,reactant.getID());

		    			hashAC.put(sO, product);
		    			
						
	    			}else{
	    				product = (Molecule) hashAC.get(sL);
	    			}
	    			
	    			IMoleculeSet productSet = builder.newMoleculeSet();
	    			//Add loss product
	    			IAtomContainer productLoss = getProductLoss(product,mzData_old);
	    			productLoss.setProperty(MZDataConstants.GROUP_PEAK_MSN, "1");
	    			productLoss.setProperty(MZDataConstants.PRECURSOR_ID,reactant.getID());
	    			if(productLoss != null){
	    				if(((IMolecule)productLoss).getProperty(CDKConstants.FORMULA) instanceof ArrayList){
	    			    	ArrayList<String> mfL = ((ArrayList<String>)productLoss.getProperty(CDKConstants.FORMULA));
		    				for(int i = 0 ; i < mfL.size(); i++){
		    					String fm = mfL.get(i).replace(" ", "");
		    					productLoss.setProperty(CDKConstants.FORMULA, MolecularFormulaManipulator.getMajorIsotopeMolecularFormula(fm, builder));
			    				break; // TODO
		    	    		}
	    				}else if(((IMolecule)productLoss).getProperty(CDKConstants.FORMULA) instanceof IMolecularFormula){
	    					String fm = MolecularFormulaManipulator.getString((IMolecularFormula) ((IMolecule)productLoss).getProperty(CDKConstants.FORMULA));
	    					productLoss.setProperty(CDKConstants.FORMULA, MolecularFormulaManipulator.getMajorIsotopeMolecularFormula(fm, builder));
	    				}
		    			productSet.addAtomContainer(productLoss);
	    			}
	    			productSet.addAtomContainer(product);
	    			reaction.setProducts(productSet);
	    			
	    			if(containProducts(hashMap.entrySet(),sO)){
    					ReactionScheme scheme1 = new ReactionScheme();
	        			schemeList.put(sO, scheme1);
	    				scheme1.setID("rs"+countScheme++);
	    	    		scheme1.addReaction(reaction);
	    	    		if(schemeList.containsKey(sOO))
	    	    			schemeList.get(sOO).add(scheme1);
	    	    	}else{
	    	    		if(schemeList.containsKey(sOO))
	    	    			schemeList.get(sOO).addReaction(reaction);
	    	    	}
	        	}
	        }
		}
		mzData.setListSpectra(cmlSpecList);
		mzData.setListReactions(reactionSet);
		return mzData;
	}

	public static MZData getMZDataNominal(Map<String, List<CMLPeak>> hashMap, MZData mzData_old) {
		int charge = 0; // seting charge molecules
		
		String intrumentRef = "";
		String protocol = "";
		List<CMLMetadata> metadataListI = mzData_old.getListSpectra().getMetadataListElements().get(0).getMetadataDescendants();
		for(CMLMetadata metadata:metadataListI){
			if(metadata.getDictRef().equals(MZDataConstants.INSTRUMENT))
				intrumentRef = metadata.getContent();
			else if(metadata.getDictRef().equals(MZDataConstants.PROTOCOL))
				protocol = metadata.getContent();
		}
		Polarity polarity = getPolarity(mzData_old);
		
		if(polarity.equals(Polarity.positive))
			charge = 1;
		else if(polarity.equals(Polarity.negative))
			charge = -1;

		String activationMethod = getActivationMethod(mzData_old);
		String window = getWindow(mzData_old);
		
		int countSp = 1;
		int countSi = 1;
		
		MZData mzData = new MZData();

		ReactionScheme reactionSet = new ReactionScheme();
		HashMap<String, ReactionScheme> schemeList = new HashMap<String, ReactionScheme>();
		int countScheme = 1;
		int countReaction = 1;
		
		
        CMLSpectrumList cmlSpecList = new CMLSpectrumList();
		CMLMetadataList metadataListL = new CMLMetadataList();
		/////////////////////////////////
		CMLMetadata metadataL = new CMLMetadata();
		metadataL.setDictRef(MZDataConstants.NUM_GROUPS);
		metadataL.setContent(Integer.toString(1));
		metadataListL.appendChild(metadataL);
		/////////////////////////////////
		if(!protocol.equals("")){
			metadataL = new CMLMetadata();
			metadataL.setDictRef(MZDataConstants.PROTOCOL);
			metadataL.setContent(protocol);
			metadataListL.appendChild(metadataL);
		}
		/////////////////////////////////
		if(!intrumentRef.equals("")){
			metadataL = new CMLMetadata();
			metadataL.setDictRef(MZDataConstants.INSTRUMENT);
			metadataL.setContent(intrumentRef);
			metadataListL.appendChild(metadataL);
		}
		/////////////////////////////////
		cmlSpecList.addMetadataList(metadataListL);
        HashMap<String, String> hashS = new HashMap<String, String>();
//        HashMap<String, IMolecule> hashAC = new HashMap<String, IMolecule>();
//        HashMap<String, IMolecule> hashAC_L = new HashMap<String, IMolecule>();//loss
        
        Set<Entry<String, List<CMLPeak>>> set = hashMap.entrySet();
        
        // List molecules
		for( Iterator<Entry<String, List<CMLPeak>>> it = set.iterator(); it.hasNext();){
			Entry<String, List<CMLPeak>> key = it.next();
			List<CMLPeak> l_l = key.getValue();
			String sOO = key.getKey();
			
			boolean init = false;

//			//reaction
//			IMolecule reactant;
//			IMolecularFormula formula;
//			if(sOO.length() == 0)
//				continue;
//			if(sOO.lastIndexOf("||") == -1){ 
//				formula = MolecularFormulaManipulator.getMolecularFormula(sOO, builder);
//			}else{
//				formula = MolecularFormulaManipulator.getMolecularFormula(sOO.substring(sOO.lastIndexOf("||"), sOO.length()), builder);
//			}
//				
//    		if(!hashAC.containsKey(sOO)){
//    			reactant = (IMolecule) l_l.get(0);
////				reactant.setID(l_l.get(0).getID());
////				formula.setCharge(charge);
//				reactant.setProperty(CDKConstants.FORMULA, formula);
////				reactant.setProperty("cdk:partialCharge", Double.toString(charge));
//				hashAC.put(sOO, reactant);
//			}else{
//				reactant = (IMolecule) hashAC.get(sOO);
//			}
			if(sOO.lastIndexOf("||") == -1){ // full scan
				
				// Spectrum 
				CMLSpectrum cmlSpect = new CMLSpectrum();
				CMLPeakList cmlPeaks = new CMLPeakList();
    			CMLMetadataList metadataList = new CMLMetadataList();
    			cmlSpect.addMetadataList(metadataList);
				CMLConditionList conditionList = new CMLConditionList();
				cmlSpect.addConditionList(conditionList);
    			cmlSpect.setAttribute("type", "MS");
    			cmlSpect.setId(Integer.toString(countSp++));
    			hashS.put(sOO, Integer.toString(countSp-1));
    			/////////////////////////////////
    			CMLMetadata metadata = new CMLMetadata();
    			metadata.setDictRef(MZDataConstants.MS_LEVEL);
    			metadata.setContent(Integer.toString(1));
    			metadataList.appendChild(metadata);
    			/////////////////////////////////
				metadata = new CMLMetadata();
				metadata.setDictRef(MZDataConstants.GROUP_PEAK_MSN);
				metadata.setContent(""+1);
				metadataList.appendChild(metadata);
    			/////////////////////////////////
				metadata = new CMLMetadata();
				metadata.setDictRef(MZDataConstants.SCAN_NUM);
				metadata.setContent(Integer.toString(countSp-1));
				metadataList.appendChild(metadata);
				/////////////////////////////////
				CMLScalar condition = new CMLScalar();		
				condition.setDictRef(MZDataConstants.PRECURSOR_MZ);
				String precurMZ = "0";
				condition.setValue(precurMZ);
				conditionList.appendChild(condition);
				/////////////////////////////////
				condition = new CMLScalar();		
				condition.setDictRef(MZDataConstants.POLARITY);
				condition.setValue(polarity.toString());
				conditionList.appendChild(condition);
				/////////////////////////////////
				condition = new CMLScalar();		
				condition.setDictRef(MZDataConstants.ACTIVATION_METHOD);
				condition.setValue(activationMethod);
				conditionList.appendChild(condition);
				/////////////////////////////////
				condition = new CMLScalar();		
				condition.setDictRef(MZDataConstants.WINDOW);
				condition.setValue(window);
				conditionList.appendChild(condition);
				/////////////////////////////////
				Double[] xyValue = getXYValueAverageP(l_l);
    			Double xValue = xyValue[0];
            	Double yValue = xyValue[1];
            	CMLPeak peak = new CMLPeak();
    			peak.setId("peak"+Integer.toString(countSi++));
    			peak.setXValue(xValue.toString());
    			peak.setYValue(yValue.toString());
//				peak.setMoleculeRefs(l_l.get(0).getID());
    			cmlPeaks.addPeak(peak); 
    	        cmlSpect.addPeakList(cmlPeaks);
    			cmlSpecList.addSpectrum(cmlSpect);
    			
//    			// reaction
//    			ReactionScheme scheme = new ReactionScheme();
//    			schemeList.put(sOO, scheme);
//    			reactionSet.add(scheme);
//    			scheme.setID("rs"+countScheme++);
    			
    		}
			
			CMLPeakList cmlPeaks = new CMLPeakList();
			
			
			Set<Entry<String, List<CMLPeak>>> setl = hashMap.entrySet();
			for( Iterator<Entry<String, List<CMLPeak>>> itl = setl.iterator(); itl.hasNext();){
				Entry<String, List<CMLPeak>> keyl = itl.next();
				List<CMLPeak> l_ll = keyl.getValue();
				String sO = keyl.getKey();
				if(sO.lastIndexOf("||") == -1){ // full scan
    				continue;
    			}
    			
	        	String sP = sO.substring(0, sO.lastIndexOf("||"));
        		if(sP.equals(sOO)){

	        		
	        		if(!init){ // initiate spectrum 

	        	        CMLSpectrum cmlSpect = new CMLSpectrum();
	        			CMLMetadataList metadataList = new CMLMetadataList();
	        			cmlSpect.addMetadataList(metadataList);
	    				CMLConditionList conditionList = new CMLConditionList();
	    				cmlSpect.addConditionList(conditionList);
	        			cmlSpect.setAttribute("type", "MS");
	        			cmlSpect.setId(Integer.toString(countSp++));
	        			/////////////////////////////////
	        			CMLMetadata metadata = new CMLMetadata();
	        			metadata.setDictRef(MZDataConstants.MS_LEVEL);
	        			metadata.setContent(Integer.toString(sO.split("\\|\\|").length));
	        			metadataList.appendChild(metadata);
	        			/////////////////////////////////
						metadata = new CMLMetadata();
						metadata.setDictRef(MZDataConstants.GROUP_PEAK_MSN);
						metadata.setContent(""+1);
						metadataList.appendChild(metadata);
	        			/////////////////////////////////
	    				metadata = new CMLMetadata();
	    				metadata.setDictRef(MZDataConstants.SCAN_NUM);
	    				metadata.setContent(Integer.toString(countSp-1));
	    				metadataList.appendChild(metadata);
	        			/////////////////////////////////
	    				CMLScalar condition = new CMLScalar();		
	    				condition.setDictRef(MZDataConstants.PRECURSOR_MZ);
	    				String[] ll = sOO.split("\\|\\|");
		        		String sL = ll[ll.length-1];
//		        		IMolecularFormula formulaP = MolecularFormulaManipulator.getMolecularFormula(sL, builder);
//		        		formulaP.setCharge(charge);
//		        		double mass = MolecularFormulaManipulator.getTotalExactMass(formulaP);
//		    			DecimalFormat df = new DecimalFormat("#.###");
		    			condition.setValue(sL);
	    				conditionList.appendChild(condition);
						/////////////////////////////////
						condition = new CMLScalar();		
						condition.setDictRef(MZDataConstants.PRECURSOR_SCAN);
						condition.setValue((String)hashS.get(sOO));
						conditionList.appendChild(condition);
						/////////////////////////////////
						condition = new CMLScalar();		
						condition.setDictRef(MZDataConstants.POLARITY);
						condition.setValue(polarity.toString());
						conditionList.appendChild(condition);
						/////////////////////////////////
						condition = new CMLScalar();		
						condition.setDictRef(MZDataConstants.ACTIVATION_METHOD);
						condition.setValue(activationMethod);
						conditionList.appendChild(condition);
						/////////////////////////////////
						condition = new CMLScalar();		
						condition.setDictRef(MZDataConstants.WINDOW);
						condition.setValue(window);
						conditionList.appendChild(condition);
						/////////////////////////////////
	        			
	        	        cmlSpect.addPeakList(cmlPeaks);
	        			cmlSpecList.addSpectrum(cmlSpect);
	        			
	        			init = true;
	        		}
        			hashS.put(sO, Integer.toString(countSp-1));
	        		// spectrum 
//	        		String[] ll = sO.split("\\|\\|");
//	        		String sL = ll[ll.length-1];
	        		
//	        		IMolecularFormula formulaP = MolecularFormulaManipulator.getMolecularFormula(sL.split("@")[0], builder);
//	        		formulaP.setCharge(charge);

					Double[] xyValue = getXYValueAverageP(l_ll);
	    			Double xValue = xyValue[0];
	            	Double yValue = xyValue[1];
	            	CMLPeak peak = new CMLPeak();
	    			peak.setId("peak"+Integer.toString(countSi++));
	    			peak.setXValue(xValue.toString());
	    			peak.setYValue(yValue.toString());
//					peak.setMoleculeRefs(l_ll.get(0).getID());
	    			cmlPeaks.addPeak(peak);

//	        		// reaction
//	    			IReaction reaction = new Reaction();
//	    	    	reaction.setID("react"+countReaction++);
//	    	    	reaction.addReactant(reactant);
//
//	    			IMolecule product;
//	    			if(!hashAC.containsKey(sO)){
//	    				product = (IMolecule) l_ll.get(0);
//		    			product.setID(l_ll.get(0).getID());
//		    			
//		    			if(sL.split("@").length == 1)
//			    			product.setProperty(CDKConstants.FORMULA, formulaP);
//		    			else{
//		    				IMolecularFormulaSet formulaSet = builder.newMolecularFormulaSet();
//		    				for(int s = 0 ; s < sL.split("@").length; s++){
//		    					IMolecularFormula formulaPS = MolecularFormulaManipulator.getMolecularFormula(sL.split("@")[s], builder);
//		    					formulaPS.setCharge(charge);
//		    	        		formulaSet.addMolecularFormula(formulaPS);
//		    				}
//		    				// TODO how to solve that
////			    			product.setProperty(CDKConstants.FORMULA, formulaSet);
//		    			}
//		    			product.setProperty("cdk:partialCharge", Double.toString(charge));
//
//		    			hashAC.put(sO, product);
//		    			
//						
//	    			}else{
//	    				product = (Molecule) hashAC.get(sL);
//	    			}
	    			
//	    			IMoleculeSet productSet = builder.newMoleculeSet();
//	    			//Add loss product
//	    			IAtomContainer productLoss = getProductLoss(product,mzData_old);
//	    			if(productLoss != null){
//	    				if(((IMolecule)productLoss).getProperty(CDKConstants.FORMULA) instanceof ArrayList){
//	    			    	ArrayList<String> mfL = ((ArrayList<String>)productLoss.getProperty(CDKConstants.FORMULA));
//		    				for(int i = 0 ; i < mfL.size(); i++){
//		    					String fm = mfL.get(i).replace(" ", "");
//		    					productLoss.setProperty(CDKConstants.FORMULA, MolecularFormulaManipulator.getMajorIsotopeMolecularFormula(fm, builder));
//			    				break; // TODO
//		    	    		}
//	    				}else if(((IMolecule)productLoss).getProperty(CDKConstants.FORMULA) instanceof IMolecularFormula){
//	    					String fm = MolecularFormulaManipulator.getString((IMolecularFormula) ((IMolecule)productLoss).getProperty(CDKConstants.FORMULA));
//	    					productLoss.setProperty(CDKConstants.FORMULA, MolecularFormulaManipulator.getMajorIsotopeMolecularFormula(fm, builder));
//	    				}
//		    			productSet.addAtomContainer(productLoss);
//	    			}
//	    			productSet.addAtomContainer(product);
//	    			reaction.setProducts(productSet);
//	    			
//	    			if(containProducts(hashMap.entrySet(),sO)){
//    					ReactionScheme scheme1 = new ReactionScheme();
//	        			schemeList.put(sO, scheme1);
//	    				scheme1.setID("rs"+countScheme++);
//	    	    		scheme1.addReaction(reaction);
//	    	    		if(schemeList.containsKey(sOO))
//	    	    			schemeList.get(sOO).add(scheme1);
//	    	    	}else{
//	    	    		if(schemeList.containsKey(sOO))
//	    	    			schemeList.get(sOO).addReaction(reaction);
//	    	    	}
	        	}
	        }
		}
		mzData.setListSpectra(cmlSpecList);
		mzData.setListReactions(reactionSet);
		return mzData;
	}
	private static IAtomContainer getProductLoss(IMolecule product,	MZData mzDataOld) {
		
    	IMoleculeSet acL = ReactionSchemeManipulator.getAllMolecules(mzDataOld.getListReactions());
    	for(IAtomContainer ac :acL.molecules()){
    		if(((String) ac.getProperty("cdk:partialCharge")).equals("0.0"))
    			if(product.getID().equals(ac.getID().replace("Loss", ""))){
//    				if(sL.split("@").length == 1)
//    					product.setProperty(CDKConstants.FORMULA, formulaP);
//    				else{
//    					IMolecularFormulaSet formulaSet = builder.newMolecularFormulaSet();
//    					for(int s = 0 ; s < sL.split("@").length; s++){
//    						IMolecularFormula formulaPS = MolecularFormulaManipulator.getMolecularFormula(sL.split("@")[s], builder);
//    						formulaPS.setCharge(charge);
//    		        		formulaSet.addMolecularFormula(formulaPS);
//    					}
//    				}
    				
    				return ac;
    			}
    	}
		return null;
	}
	private static Double[] getXYValueAverage(List<IAtomContainer> lAC, List<MZData> mzDataList) {
		double sumX = 0;
		double sumY = 0;
		List<String> listID = new ArrayList<String>();
		for(IAtomContainer ac : lAC)
			listID.add(ac.getProperty("mzData")+"@"+ac.getID());
//		listID.add(ac.getID());
		
		for(MZData mzData:mzDataList){
			CMLSpectrumList cmlSpectList = mzData.getListSpectra();
			CMLElements<CMLSpectrum> specElem = cmlSpectList.getSpectrumElements();
			for(CMLSpectrum spectrum:specElem){
				for(CMLPeak peak:CMLSpectrum.getDescendantPeaks(spectrum)){
					if(peak.getMoleculeRefs() != null){
						if(listID.contains(mzData.getID()+"@"+peak.getMoleculeRefs()[0])){
							sumX += peak.getXValue();
							Object test = peak.getMoleculeRefs()[0];
							sumY += peak.getYValue();
						}
					}
				}
			}
		}
		Double[] result = {sumX/lAC.size(),sumY/lAC.size()};
		return result;
	}

	private static Double[] getXYValueAverageP(List<CMLPeak> lp) {
		double sumX = 0;
		double sumY = 0;
//		List<String> listID = new ArrayList<String>();
		for(CMLPeak peak:lp){
			sumX += peak.getXValue();
			sumY += peak.getYValue();
		}
		Double[] result = {sumX/lp.size(),sumY/lp.size()};
		return result;
	}
	private static String getActivationMethod(MZData mzData) {
		CMLConditionList conditionList = mzData.getListSpectra().getSpectrumElements().get(0).getConditionListElements().get(0);
		Iterator conditionIterator = conditionList.getScalarElements().iterator();
		while (conditionIterator.hasNext()) {
			CMLScalar condition = (CMLScalar) conditionIterator.next();
			if(condition.getDictRef().equals(MZDataConstants.ACTIVATION_METHOD) ){
				String activationMethod = condition.getValue();
				return activationMethod;
			}
		}
		return null;
	}

	private static String getWindow(MZData mzData) {
		CMLConditionList conditionList = mzData.getListSpectra().getSpectrumElements().get(0).getConditionListElements().get(0);
		Iterator conditionIterator = conditionList.getScalarElements().iterator();
		while (conditionIterator.hasNext()) {
			CMLScalar condition = (CMLScalar) conditionIterator.next();
			if(condition.getDictRef().equals(MZDataConstants.WINDOW) ){
				String activationMethod = condition.getValue();
				return activationMethod;
			}
		}
		return null;
	}
	/**
	 * looking if contains products
	 * 
	 * @param seta
	 * @param sOO
	 * @return
	 */
	private static boolean containProducts(Set<Entry<String, List<IAtomContainer>>> seta, String sOO) {
		for( Iterator<Entry<String, List<IAtomContainer>>> ita = seta.iterator(); ita.hasNext();){
			Entry<String, List<IAtomContainer>> keyl = ita.next();
			String sO = keyl.getKey();
			if(sO.lastIndexOf("||") == -1)
				continue;
			String sP = sO.substring(0, sO.lastIndexOf("||"));
        	if(sP.equals(sOO)){
        		return true;
        	}
		}
		return false;
	}
//	/**
//	 * Get a MZData from a matrix object simple containing the information.
//	 * 0:ID, 1:IDparent, 2:level, 3:rt 4:m/z, 5:int, 6:file, 7:group
//	 * 
//	 * @param data   The matrix (ID, precuID, mass, int)
//	 * @param mzData The MZData
//	 * @return       The MZData object fill with the data
//	 */
//	public static MZData getMZData2(double[][] data, MZData mzData) {
//		return getMZData2(data, mzData, null, "", "",1);
//	}
//	/**
//	 * Get a MZData from a matrix object containing the information.
//	 * 0:ID, 1:IDparent, 2:level, 3:rt 4:m/z, 5:int, 6:file, 7:group
//	 * 
//	 * @param data   The matrix
//	 * @param mzData The MZData
//	 * @param pol    The polarity of the data
//	 * @return       The MZData object fill with the data
//	 */
//	public static MZData getMZData2(double[][] data, MZData mzData, int pol) {
//		return getMZData2(data, mzData, null, "", "",pol);
//	}
//	/**
//	 * Get a MZData from a matrix object containing the information.
//	 * 0:ID, 1:IDparent, 2:level, 3:rt 4:m/z, 5:int, 6:file, 7:group
//	 * 
//	 * @param data   The matrix
//	 * @param mzData The MZData
//	 * @param intrumentRef
//	 * @param protocol
//	 * @param pol    The polarity of the data
//	 * @return       The MZData object fill with the data
//	 */
//	public static MZData getMZData2(double[][] data, MZData mzData, String intrumentRef, String protocol, int pol) {
//		return getMZData2(data, mzData, null, intrumentRef, protocol, pol);
//	}
	/**
	 * Get a MZData from a matrix object containing the information. The matrix contains
	 * 0:ID, 1:IDparent, 2:level, 3:rt 4:m/z, 5:int, 6:file, 7:group
	 * 
	 * @param data   The matrix
	 * @param mzData The MZData
	 * @param pathFile
	 * @param intrumentRef
	 * @param protocol
	 * @param pol    The polarity of the data
	 * @return       The MZData object fill with the data
	 */
	public static MZData getMZData2(double[][] data, MZData mzData, String pathFile, String intrumentRef, String protocol, int pol) {
		
		CMLSpectrumList specList = new CMLSpectrumList();
		CMLMetadataList metadataListL = new CMLMetadataList();
		specList.addMetadataList(metadataListL);
		
		int countScans = 0;
		int numGroups = 0;
		try {
			DocumentBuilderFactory docBuilderFactory = DocumentBuilderFactory.newInstance();
	        DocumentBuilder docBuilder = docBuilderFactory.newDocumentBuilder();
	        NodeList listOfScans =  null;
	        Document doc = null;
	    	if(pathFile != null){
		        doc = docBuilder.parse(pathFile);
		        // normalize text representation
		        doc.getDocumentElement().normalize();
		        listOfScans = doc.getElementsByTagName("scan");
        	}
			// restructure the links
			for(int i = 0; i < data.length; i++){
				double precID = data[i][1];
				if(precID == 0.0)
					continue;
				// extraction of the rt 
				double rt = -1;
				double precMass = 0.0;
				for(int j = 0; j < data.length; j++){
					if(data[j][0] == precID){
						precMass = data[j][4];
						rt = data[j][3];
						break;
					}
				}
				// look for those peaks with the same rt and a window of CID (0.5) for the highest peak
				double intensity = 0;
				for(int j = 0; j < data.length; j++){
					if(data[j][3] == rt){
						double max = precMass+0.5;
						double min = precMass-0.5;
						double mass = data[j][4];
						if(mass > max || mass < min)
							continue;
						if(intensity < data[j][5]){
							intensity = data[j][5];
							precID = data[j][0];
						}
					}
				}
				data[i][1] = precID;
			}
			DecimalFormat df = new DecimalFormat("#");
	        double rt_old = -1;
			CMLPeakList cmlPeaks = null;
			String polarity = "unknown";
			HashMap<String, String> mapScan = new HashMap<String, String>();
			for(int i = 0; i < data.length; i++){
				double rt = data[i][3];
				if(rt != rt_old){
					ScanMZXML scanMZXML = getScanMZXML(rt,listOfScans);
////					c.eval("scan_pos <- which(xr@scantime == "+(String) data[i][3]+")");
////					c.eval("scan_pos_msn <- which(xr@msnRt == "+(String) data[i][3]+")");
					rt_old = rt;
					CMLSpectrum cmlSpect = new CMLSpectrum();
					cmlSpect.setAttribute("type", "MS");
					CMLConditionList conditionList = new CMLConditionList();
					cmlSpect.addConditionList(conditionList);
					CMLMetadataList metadataList = new CMLMetadataList();
					cmlSpect.addMetadataList(metadataList);
					countScans++;
					cmlSpect.setId(Integer.toString(countScans));
					specList.addSpectrum(cmlSpect);
					/////////////////////////////////
					CMLMetadata metadata = new CMLMetadata();
					metadata.setDictRef(MZDataConstants.SCAN_NUM);
					if(pathFile != null)
				        metadata.setContent(scanMZXML.getNum());
					else
					    metadata.setContent(Integer.toString(countScans));
					        
					metadataList.appendChild(metadata);
					/////////////////////////////////
					metadata = new CMLMetadata();
					metadata.setDictRef(MZDataConstants.MS_LEVEL);
					Double nLevel = data[i][2];  
					
					if(nLevel.equals(1.0))
						numGroups++;
					metadata.setContent(df.format(nLevel));
					metadataList.appendChild(metadata);
					/////////////////////////////////
					metadata = new CMLMetadata();
					metadata.setDictRef(MZDataConstants.RT);
					metadata.setContent((new Double(data[i][3])).toString());
					metadataList.appendChild(metadata);
					mapScan.put((new Double(data[i][3])).toString(), "" + (countScans));
					/////////////////////////////////
					metadata = new CMLMetadata();
					metadata.setDictRef(MZDataConstants.GROUP_PEAK_MSN);
					metadata.setContent(df.format(data[i][7]));
					metadataList.appendChild(metadata);
					/////////////////////////////////
					CMLScalar condition = new CMLScalar();		
					condition.setDictRef(MZDataConstants.PRECURSOR_MZ);
					double precurMZ = 0.0;
					double precurRT = 0.0;
					double precurID = data[i][1];
					for(int j = 0; j < data.length; j++){
						if(precurID == data[j][0]){
							precurMZ = data[j][4];
							precurRT = data[j][3];
							break;
						}
					}
					if(precurMZ != 0.0){
						condition.setValue(precurMZ);
						conditionList.appendChild(condition);
					}
					/////////////////////////////////
					condition = new CMLScalar();		
					condition.setDictRef(MZDataConstants.PRECURSOR_SCAN);
					if(precurMZ != 0.0){
//					if(c.eval("(length(scan_pos_msn) > 0) != 0").asString().equals("true")){
						condition.setValue(mapScan.get((new Double(precurRT)).toString()));
						conditionList.appendChild(condition);
					}
					/////////////////////////////////
					condition = new CMLScalar();		
					condition.setDictRef(MZDataConstants.WINDOW);
					condition.setValue("1");
					conditionList.appendChild(condition);
					/////////////////////////////////
					condition = new CMLScalar();		
					condition.setDictRef(MZDataConstants.COLLISION_ENERGY);
					if(precurMZ != 0.0){
//					String ce = "0";
////					if(c.eval("(length(scan_pos_msn) > 0) != 0").asString().equals("true")){
////						ce = c.eval("xr@msnCollisionEnergy[scan_pos_msn]").asString();
						if(pathFile != null){
							condition.setValue(scanMZXML.getcollisionEnergy());
							conditionList.appendChild(condition);
						}
					}
					/////////////////////////////////
					condition = new CMLScalar();		
					condition.setDictRef(MZDataConstants.ACTIVATION_METHOD);
					condition.setValue("CID");
					conditionList.appendChild(condition);
//					/////////////////////////////////
					condition = new CMLScalar();		
					condition.setDictRef(MZDataConstants.POLARITY);
					if(precurMZ == 0.0){
						if(pathFile != null){
						    if(scanMZXML.getpolarity().equals("+"))
								polarity = "positive";
							else if(scanMZXML.getpolarity().equals("-"))
								polarity = "negative";
							else
								polarity = "unknown";
						}else{
							if(pol == 1)
								polarity = "positive";
							else if(pol == -1)
								polarity = "negative";
							else
								polarity = "unknown";
						}
					}
//					if(c.eval("(length(scan_pos) > 0) != 0").asString().equals("true"))
//						polarity = c.eval("xr@polarity[scan_pos]").asString();
					condition.setValue(polarity);
					conditionList.appendChild(condition);
//					/////////////////////////////////
					cmlPeaks = new CMLPeakList();

					// if it is full scan
					if(nLevel == 1){
//					if(c.eval("(length(scan_pos) > 0) != 0").asString().equals("true")){
////						c.eval("scan = getScan(xr,scan_pos)");
//////						c.eval("peaks <- specPeaks(scan, sn = 0.7, mzgap = 0.01)");
////						c.eval("peaks <- scan"); // TMP
////						c.eval("ll <- order(peaks[,1])");
////						c.eval("peaks <- peaks[ll,]");
						String id = df.format(data[i][0]);
////						if(c.eval("length(ll)==1").asString().equals("true")){
							CMLPeak peak = new CMLPeak();
//							String[] values = c.eval("peaks").asStrings();
							peak.setId(df.format(data[i][0]));
							peak.setXValue((new Double(data[i][4])).toString());
							peak.setYValue((new Double(data[i][5])).toString());
							cmlPeaks.addPeak(peak);
////						
////						}else{
////							int numRow2  = c.eval("length(peaks[,1])").asInteger();
////							for(int k = 1; k <= numRow2; k++){
////								CMLPeak peak = new CMLPeak();
////								String[] values = c.eval("peaks["+k+",]").asStrings();
////								peak.setId(id+"."+k);
////								peak.setXValue((String) values[0]);
////								
////								peak.setYValue((String) values[1]);
////								cmlPeaks.addPeak(peak);
////							}
////						}
						cmlSpect.addPeakList(cmlPeaks);
						
					}else{
						CMLPeak peak = new CMLPeak();
						peak.setId(df.format(data[i][0]));
						peak.setXValue((new Double(data[i][4])).toString());
						peak.setYValue((new Double(data[i][5])).toString());
						cmlPeaks.addPeak(peak);
						cmlSpect.addPeakList(cmlPeaks);
					}
				}else{
					CMLPeak peak = new CMLPeak();
					peak.setId(df.format(data[i][0]));
					peak.setXValue((new Double(data[i][4])).toString());
					peak.setYValue((new Double(data[i][5])).toString());
					cmlPeaks.addPeak(peak);
					
				}
//				break;
			}

			/////////////////////////////////
			CMLMetadata metadataL = new CMLMetadata();
			metadataL.setDictRef(MZDataConstants.SCANS);
			metadataL.setContent(Integer.toString(countScans));
			metadataListL.appendChild(metadataL);
			/////////////////////////////////
			metadataL = new CMLMetadata();
			metadataL.setDictRef(MZDataConstants.NUM_GROUPS);
			metadataL.setContent(Integer.toString(numGroups));
			metadataListL.appendChild(metadataL);
			/////////////////////////////////
			if(!protocol.equals("")){
				metadataL = new CMLMetadata();
				metadataL.setDictRef(MZDataConstants.PROTOCOL);
				metadataL.setContent(protocol);
				metadataListL.appendChild(metadataL);
			}
			/////////////////////////////////
			if(!intrumentRef.equals("")){
				metadataL = new CMLMetadata();
				metadataL.setDictRef(MZDataConstants.INSTRUMENT);
				metadataL.setContent(intrumentRef);
				metadataListL.appendChild(metadataL);
			}
			/////////////////////////////////parentFile RAW
			if(pathFile != null){
			    metadataL = new CMLMetadata();
				metadataL.setDictRef(MZDataConstants.PARENT_FILE);
				NodeList listOfModels = doc.getElementsByTagName("parentFile");
		        Element modelElement = (Element)listOfModels.item(0);
		        if(modelElement != null){
		        	String fileName = modelElement.getAttribute("fileName");
	//	        	if(fileName.endsWith(".RAW"))
	//	        		fileName = fileName.replace(".RAW", ".mzXML");
	//	        	else if(fileName.endsWith(".raw"))
	//	        		fileName = fileName.replace(".raw", ".mzXML");
		        	metadataL.setContent(fileName);
		        	metadataL.setConvention("RAWData");
					metadataListL.appendChild(metadataL);
		        }
				/////////////////////////////////parentFile mzXML
				metadataL = new CMLMetadata();
				metadataL.setDictRef(MZDataConstants.PARENT_FILE);
		        metadataL.setContent(pathFile);
		        metadataL.setConvention("mzXMLData");
				metadataListL.appendChild(metadataL);
			}
//		} catch (RserveException e) {
//			e.printStackTrace();
//		} catch (REXPMismatchException e) {
//			e.printStackTrace();
//		} catch (ParserConfigurationException e) {
//			e.printStackTrace();
//		} catch (SAXException e) {
//			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (SAXException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ParserConfigurationException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		mzData.setListSpectra(specList);
		Map<Object,Object> properties = new HashMap<Object, Object>();
		properties.put(MZDataConstants.NUM_GROUPS, Integer.toString(numGroups));
		properties.put(MZDataConstants.SCANS, Integer.toString(countScans));
		PackageVersion v = new PackageVersion();
    	if(v.existManifest()){
    		String version = v.getAttribute(PackageVersion.ATTRIBUTE_IMPLEMENTATIONVERSION);
    		properties.put(MZDataConstants.VERSION_SAMS, version);
    	}
    	mzData.setProperties(properties);
		return mzData;
	}
	/**
	 * Extract number of scans based on the different number
	 * of retention times
	 * 
	 * @param data
	 * @return
	 */
	private static int extractListScans(double[][] data) {
		List<String> rtList = new ArrayList<String>(); 
		for(int j = 0; j < data.length; j++){
			String rt = Double.toString(data[j][3]);
			if(!rtList.contains(rt)){
				rtList.add(rt);
			}
		}
		
		return rtList.size();
	}
	public static ScanMZXML getScanMZXML(double rt, NodeList listOfModels2) {
		
		String rtN = (new Double(rt)).toString();
		if(rtN.endsWith(".0"))
			rtN = rtN.replace(".0", "");

		if(listOfModels2 != null && listOfModels2.getLength() > 0) {
			for(int i = 0 ; i < listOfModels2.getLength();i++) 
				if(listOfModels2.item(i) != null){
					Element modelElement3 = (Element)listOfModels2.item(i);
			        String rtt = modelElement3.getAttribute("retentionTime");
			        if(rtt.contains(rtN)){
			        	return new ScanMZXML(modelElement3.getAttribute("num"), modelElement3.getAttribute("polarity"), modelElement3.getAttribute("collisionEnergy"));
			        }
				}
		}
		return null;
	}
	/**
	 * Get a MZData from a matrix object simple containing the information.
	 * 0:ID, 1:IDparent, 2:level, 3:rt 4:m/z, 5:int, 6:file, 7:group
	 * 
	 * @param matrix
	 * @param mzData
	 * @return
	 */
	public static MZData getMZData(double[][] matrix, MZData mzData) {
		int countSp = 1;
		int countSi = 1;
		
		List<Integer> levelL = extractLevels(matrix);
		
        CMLSpectrumList cmlSpecList = new CMLSpectrumList();
		CMLMetadataList metadataListL = new CMLMetadataList();
		/////////////////////////////////
		CMLMetadata metadataL = new CMLMetadata();
		metadataL.setDictRef(MZDataConstants.NUM_GROUPS);
		metadataL.setContent(Integer.toString(1));
		metadataListL.appendChild(metadataL);
		/////////////////////////////////
		cmlSpecList.addMetadataList(metadataListL);
        HashMap<String, String> hashS = new HashMap<String, String>();
                
        for(int i = 0 ; i < matrix.length; i++){

			boolean init = false;
			
			if(i == 0){ // full scan
				
				CMLSpectrum cmlSpect = new CMLSpectrum();
				CMLPeakList cmlPeaks = new CMLPeakList();
    			CMLMetadataList metadataList = new CMLMetadataList();
    			cmlSpect.addMetadataList(metadataList);
				CMLConditionList conditionList = new CMLConditionList();
				cmlSpect.addConditionList(conditionList);
    			cmlSpect.setAttribute("type", "MS");
    			cmlSpect.setId(Integer.toString(countSp++));
    			/////////////////////////////////
    			CMLMetadata metadata = new CMLMetadata();
    			metadata.setDictRef(MZDataConstants.MS_LEVEL);
    			metadata.setContent(Integer.toString(1));
    			metadataList.appendChild(metadata);
    			/////////////////////////////////
				metadata = new CMLMetadata();
				metadata.setDictRef(MZDataConstants.GROUP_PEAK_MSN);
				metadata.setContent(""+1);
				metadataList.appendChild(metadata);
    			/////////////////////////////////
				metadata = new CMLMetadata();
				metadata.setDictRef(MZDataConstants.SCAN_NUM);
				metadata.setContent(Integer.toString(countSp-1));
				metadataList.appendChild(metadata);
				/////////////////////////////////
				CMLScalar condition = new CMLScalar();		
				condition.setDictRef(MZDataConstants.PRECURSOR_MZ);
				String precurMZ = "0";
				hashS.put(Integer.toString(i), Integer.toString(countSp-1));
				condition.setValue(precurMZ);
				conditionList.appendChild(condition);
				/////////////////////////////////
				condition = new CMLScalar();		
				condition.setDictRef(MZDataConstants.POLARITY);
				Polarity polarity = Polarity.nothing;
				if(matrix[i][4] == 1.0)
					polarity = Polarity.positive;
				else if(matrix[i][4] == -1.0)
					polarity = Polarity.negative;
				condition.setValue(polarity.toString());
				conditionList.appendChild(condition);
				/////////////////////////////////
    			
    			Double xValue = matrix[i][2];
            	Double yValue = matrix[i][3];
            	CMLPeak peak = new CMLPeak();
    			peak.setId(Integer.toString(countSi++));
    			peak.setXValue(xValue.toString());
    			peak.setYValue(yValue.toString());
    			cmlPeaks.addPeak(peak); 
    	        cmlSpect.addPeakList(cmlPeaks);
    			cmlSpecList.addSpectrum(cmlSpect);
			}
			CMLPeakList cmlPeaks = new CMLPeakList();
			
			for(int j = 0 ; j < matrix.length; j++){
				if(j == 0){ // full scan
    				continue;
    			}
    			
				if(matrix[i][0] == matrix[j][1]){
					
					if(!init){ // initiate spectrum 

	        	        CMLSpectrum cmlSpect = new CMLSpectrum();
	        			CMLMetadataList metadataList = new CMLMetadataList();
	        			cmlSpect.addMetadataList(metadataList);
	    				CMLConditionList conditionList = new CMLConditionList();
	    				cmlSpect.addConditionList(conditionList);
	        			cmlSpect.setAttribute("type", "MS");
	        			cmlSpect.setId(Integer.toString(countSp++));
	        			/////////////////////////////////
	        			CMLMetadata metadata = new CMLMetadata();
	        			metadata.setDictRef(MZDataConstants.MS_LEVEL);
//	        			metadata.setContent(Integer.toString(sO.split("\\|\\|").length));
	        			metadata.setContent(Integer.toString(levelL.get(j)+1));
	        			metadataList.appendChild(metadata);
	        			/////////////////////////////////
						metadata = new CMLMetadata();
						metadata.setDictRef(MZDataConstants.GROUP_PEAK_MSN);
						metadata.setContent(""+1);
						metadataList.appendChild(metadata);
	        			/////////////////////////////////
	    				metadata = new CMLMetadata();
	    				metadata.setDictRef(MZDataConstants.SCAN_NUM);
	    				metadata.setContent(Integer.toString(countSp-1));
	    				metadataList.appendChild(metadata);
	        			/////////////////////////////////
	    				CMLScalar condition = new CMLScalar();		
	    				condition.setDictRef(MZDataConstants.PRECURSOR_MZ);
		    			DecimalFormat df = new DecimalFormat("#.###");
		    			double mass = matrix[i][2];
		    			condition.setValue(df.format(mass));
	    				conditionList.appendChild(condition);
						/////////////////////////////////
	    				condition = new CMLScalar();		
						condition.setDictRef(MZDataConstants.PRECURSOR_SCAN);
						condition.setValue((String)hashS.get(Integer.toString(i)));
	        			conditionList.appendChild(condition);
						/////////////////////////////////
						condition = new CMLScalar();		
						condition.setDictRef(MZDataConstants.POLARITY);
						Polarity polarity = Polarity.nothing;
		    			if(matrix[j][4] ==  1){
		        			polarity = Polarity.positive;
		    			}else if(matrix[j][4] ==  -1){
		        			polarity = Polarity.negative;
	        			}
						condition.setValue(polarity.toString());
						conditionList.appendChild(condition);
						/////////////////////////////////
	        			
	        	        cmlSpect.addPeakList(cmlPeaks);
	        			cmlSpecList.addSpectrum(cmlSpect);
	        			
	        			init = true;
	        		}
	        		Double xValue = matrix[j][2];
	            	Double yValue = matrix[j][3];
	            	CMLPeak peak = new CMLPeak();
	            	hashS.put(Integer.toString(j), Integer.toString(countSp-1));
					peak.setId(Integer.toString(countSi++));
	    			peak.setXValue(xValue.toString());
	    			peak.setYValue(yValue.toString());
	    			cmlPeaks.addPeak(peak);
	        	}
	        }
        }
		/////////////////////////////////
		metadataL = new CMLMetadata();
		metadataL.setDictRef(MZDataConstants.SCANS);
		metadataL.setContent(Integer.toString(countSp-1));
		metadataListL.appendChild(metadataL);
		/////////////////////////////////
		
        mzData.setListSpectra(cmlSpecList);
		return mzData;
	}
	/**
	 * Extract the levels from the matrix containing 
	 * the relations between the ions.
	 * 
	 * @param matrix containing the relation
	 * @return A list with the levels
	 */
	public static List<Integer> extractLevels(double[][] matrix) {
		List<Integer> levelsL = new ArrayList<Integer>();
		for(int i = 0 ; i < matrix.length ; i++){
			double id = matrix[i][1];
			List<Double> precurL = extractPrecursors(id,matrix);
			levelsL.add(precurL.size());
		}
		return levelsL;
	}
	/**
	 * Extract all precursor ID's 
	 * 
	 * @param precur ID of the ion to investigate
	 * @param matrix
	 * @return
	 */
	public static List<Double> extractPrecursors(double id,double[][] matrix) {

		List<Double> precurLS = new ArrayList<Double>();
		for(int i = 0 ; i < matrix.length ; i++){
			if(matrix[i][0] == id){
				precurLS.add(matrix[i][0]);
				List<Double> precurL = extractPrecursors(matrix[i][1],matrix);
				precurLS.addAll(precurL);
				break;		
			}
		}
		return precurLS;
	}
	/**
	 * Print the mzData object. It prints only when it is 
	 * containing some molecule objects
	 * 
	 * @param mzData The MZData object
	 * @return  The ouput as a string
	 */
	public static String toString(MZData mzData) {
		return toString(mzData,true);
	}
	
	/**
	 * Print the mzData object. It prints only when it is 
	 * containing some molecule objects
	 * 
	 * @param mzData The MZData object
	 * @param treeStructure True, if will be in tree structure
	 * @return  The ouput as a string
	 */
	public static String toString(MZData mzData,boolean treeStructure) {
		String str = "";
		CMLSpectrumList cmlSpectList = mzData.getListSpectra();
		List<CMLMetadata> metadataList = cmlSpectList.getMetadataListElements().get(0).getMetadataDescendants();
		int numGroups = 0;
		for(CMLMetadata metadata:metadataList){
			if(metadata.getDictRef().equals(MZDataConstants.NUM_GROUPS))
				numGroups = Integer.parseInt(metadata.getContent());
		}
		CMLElements<CMLSpectrum> specElem = cmlSpectList.getSpectrumElements();
		IMoleculeSet molecules = null;
		if(mzData.getListReactions() != null)
			molecules = ReactionSchemeManipulator.getAllMolecules(mzData.getListReactions());
		for(int i = 0 ; i < numGroups; i++){
			List<CMLSpectrum> newCMLSpeList = new ArrayList<CMLSpectrum>();
			String scanIDP = "";
			for(CMLSpectrum spectrum:specElem){
				List<CMLMetadata> ml = spectrum.getMetadataListElements().get(0).getMetadataDescendants();
				boolean validS = false;
				for(CMLMetadata metadata:ml){
					if(metadata.getDictRef().equals(MZDataConstants.GROUP_PEAK_MSN) && metadata.getContent().equals((new Integer(i+1)).toString())){
						newCMLSpeList.add(spectrum);
						validS = true;
					}
				}
				if(validS){
					for(CMLMetadata metadata:ml){
						if(metadata.getDictRef().equals(MZDataConstants.MS_LEVEL) && metadata.getContent().equals("1")){
							scanIDP = spectrum.getId();
						}
					}
				}
			}
			/////////////////////////////////////////////////////////
			ParentIon parentIon = extractParentIon(molecules,newCMLSpeList,scanIDP,mzData);
			if(parentIon == null)
	        	continue;
	        CMLPeak cmlPeak = getCMLPeakFromIDP(newCMLSpeList, scanIDP);
	        IsotopePattern isoP = getIsotopePatternFromIDP(newCMLSpeList,cmlPeak);
	        str += printOutPutCompress(parentIon,isoP,treeStructure);
		}
		return str;
	}
		/**
		 * Print the output in a tree way. 
		 * 
		 * @param ionFinal
		 * @param isoP
		 * @param treeStructure
		 * @return
		 */
		private static String printOutPutCompress(ParentIon ionFinal, IsotopePattern isoP, boolean treeStructure) {
			String output = "";
			DecimalFormat Currency = new DecimalFormat("#0");
	        String formated_1 = Currency.format(ionFinal.getMass());
			if(treeStructure)
				output = printingLineTree(ionFinal,isoP);
			else {
        		String nmPath1 = formated_1;
				output = printingLineTable(ionFinal,"0",nmPath1);
			}
		    
		    for(int f2 = 0 ; f2 < ionFinal.getFragments().size(); f2++){
	        	ParentIon ion2 = ionFinal.getFragments().get(f2);
		        String formated_2 = Currency.format(ion2.getMass());
	        	if(treeStructure)
					output += printingLineTree(ion2,null);
	        	else {
            		String nmPath2 = formated_1+"||"+formated_2;
					output += printingLineTable(ion2,ionFinal.getID(),nmPath2);
	        	}
	        	for(int f3 = 0 ; f3 < ion2.getFragments().size(); f3++){
	            	ParentIon ion3 = ion2.getFragments().get(f3);
    		        String formated_3 = Currency.format(ion3.getMass());
	            	if(treeStructure)
	    				output += printingLineTree(ion3, null);
	            	else {
	            		String nmPath3 = formated_1+"||"+formated_2+"||"+formated_3;
	    				output += printingLineTable(ion3,ion2.getID(),nmPath3);
	            	}
		        	for(int f4 = 0 ; f4 < ion3.getFragments().size(); f4++){
	                	ParentIon ion4 = ion3.getFragments().get(f4);
	    		        String formated_4 = Currency.format(ion4.getMass());
	                	if(treeStructure)
		    				output += printingLineTree(ion4, null);
		            	else {
    	            		String nmPath4 = formated_1+"||"+formated_2+"||"+formated_3+"||"+formated_4;
		    				output += printingLineTable(ion4,ion3.getID(),nmPath4);
		            	}
		            	for(int f5 = 0 ; f5 < ion4.getFragments().size(); f5++){
	                    	ParentIon ion5 = ion4.getFragments().get(f5);
    	    		        String formated_5 = Currency.format(ion5.getMass());
	                    	if(treeStructure)
	    	    				output += printingLineTree(ion5, null);
	    	            	else {
        	            		String nmPath5 = formated_1+"||"+formated_2+"||"+formated_3+"||"+formated_4+"||"+formated_5;
	    	    				output += printingLineTable(ion5,ion4.getID(),nmPath5);
	    	            	}
	    	            	for(int f6 = 0 ; f6 < ion5.getFragments().size(); f6++){
	                        	ParentIon ion6 = ion5.getFragments().get(f6);
	                        	String formated_6 = Currency.format(ion6.getMass());
        	    				if(treeStructure)
	        	    				output += printingLineTree(ion6, null);
	        	            	else {
	        	            		String nmPath6 = formated_1+"||"+formated_2+"||"+formated_3+"||"+formated_4+"||"+formated_5+"||"+formated_6;
            	    				output += printingLineTable(ion6,ion5.getID(),nmPath6);
	        	            	}
	        	            	for(int f7 = 0 ; f7 < ion6.getFragments().size(); f7++){
	                            	ParentIon ion7 = ion6.getFragments().get(f7);
	                            	if(treeStructure)
	            	    				output += printingLineTree(ion7, null);
	            	            	else {
	            	    		        String formated_7 = Currency.format(ion7.getMass());
	            	    		        String nmPath7 = formated_1+"||"+formated_2+"||"+formated_3+"||"+formated_4+"||"+formated_5+"||"+formated_6+"||"+formated_7;
	            	    				output += printingLineTable(ion7,ion6.getID(),nmPath7);
	            	            	}
	                        	}
	                    	}
	                	}
	            	}
	        	}
		    }
		   return output;
		}

//		private static String printingLine(ParentIon ionFinal, boolean treeStructure) {
//			if(treeStructure)
//				return printingLineTree(ionFinal, null);
//			else
//				return printingLineTable(ionFinal);
//			
//		}
		/**
		 * print Line with a tree form 
		 * @param ionFinal The ParentIon
		 * @param isoP 
		 */
		private static String printingLineTree(ParentIon ionFinal, IsotopePattern isoP) {
			
			String output = "";
			String form = "";
			String spaceIni = "";
	    	double massExp = ionFinal.getMass();
			
			if(ionFinal.getLevel() != 1)
				for(int i = 0; i < ionFinal.getLevel()-1; i++)
					spaceIni += "   ";
			if(ionFinal.getFormulaSet() != null)
			    for(int i = 0; i < ionFinal.getFormulaSet().size(); i++){
			    	
			    	IMolecularFormula formula = ionFinal.getFormulaSet().getMolecularFormula(i);
			    	String formulaS = MolecularFormulaManipulator.getString(formula);
			    	IMolecularFormula nformula = MolecularFormulaManipulator.getMolecularFormula(formulaS, builder);
			    	double charge = Double.valueOf((String) ionFinal.getMolecule().getProperty("cdk:partialCharge"));
					nformula.setCharge((int)charge);
			    	double diff = (massExp-MolecularFormulaManipulator.getTotalExactMass(nformula));
			    	
			    	DecimalFormat df1 = new DecimalFormat("####.0");
			    	String result = df1.format(diff*1000000/massExp);
			    	form += formulaS+"<ACC("+result+")";
			    	if(isoP != null){
						IsotopePatternSimilarity is = new IsotopePatternSimilarity();
			    		IsotopePatternGenerator isotopeGe = new IsotopePatternGenerator(minAbISO);
		    			IsotopePattern patternIsoPredicted = isotopeGe.getIsotopes(nformula);
			    		is.seTolerance(toleranceISO);
			    		double tempScore = is.compare(isoP, patternIsoPredicted);
			    		form += "/ISO("+tempScore+")";
			    	}
			    	
			    	if(i != ionFinal.getFormulaSet().size()-1)
			    		form += ", ";
			    }
			else form += " ? ";
		    output += spaceIni+ionFinal.getID()+":"+massExp+"("+ionFinal.getIntensity_Abs()+")"+" ["+form+"] ";
		    IMolecularFormulaSet lossSet = ionFinal.getFormulaLossSet();
//		    System.out.println(lossSet.size());
		    if(ionFinal.getID().equals("1_")){
		    	output += "\n";
		    }else{
		    	output += " > ";
		    	if(lossSet == null){
		    		output += "\n";
		    	}else if(lossSet.size() == 1){
		    		output += MolecularFormulaManipulator.getString(lossSet.getMolecularFormula(0))+"\n";
		    	}else{
				    	for(IMolecularFormula formulaLoss : lossSet.molecularFormulas()){
				    		output += MolecularFormulaManipulator.getString(formulaLoss) +", ";
				    	}
				    	output += "\n";
				    }
		    }
			return output;
		}
		

		/**
		 * print Line with a tree form 
		 * @param ionFinal The ParentIon
		 * @param isoP 
		 */
		private static String printingLineTable(ParentIon ionFinal,String precursorID, String nomPath) {
			String output = "";
			String form = "";
	    	double massExp = ionFinal.getMass();
	    	int level = ionFinal.getLevel();
			double charge = 0;
			if(ionFinal.getFormulaSet() != null){
			    for(int i = 0; i < ionFinal.getFormulaSet().size(); i++){
			    	IMolecularFormula formula = ionFinal.getFormulaSet().getMolecularFormula(i);
			    	String formulaS = MolecularFormulaManipulator.getString(formula);
			    	IMolecularFormula nformula = MolecularFormulaManipulator.getMolecularFormula(formulaS, builder);

					charge = Double.valueOf((String) ionFinal.getMolecule().getProperty("cdk:partialCharge"));
					nformula.setCharge((int)charge);
			    	double diff = (massExp-MolecularFormulaManipulator.getTotalExactMass(nformula));
			    	
			    	DecimalFormat df1 = new DecimalFormat("####.0");
			    	String result = df1.format(diff*1000000/massExp);
			    	form += formulaS+","+result;
			    	
			    	break;
			    }
			}else form += "?,?";
			
			IMolecularFormulaSet lossSet = ionFinal.getFormulaLossSet();
			String formNL = "?";
		    if(ionFinal.getID().equals("1_")){
		    	formNL = ",";
		    }else{
		    	if(lossSet == null){
		    		formNL = "?";
		    	}else if(lossSet.size() == 1){
		    		formNL = MolecularFormulaManipulator.getString(lossSet.getMolecularFormula(0));
		    	}else{
		    		formNL = "";
			    	for(IMolecularFormula formulaLoss : lossSet.molecularFormulas()){
			    		formNL += MolecularFormulaManipulator.getString(formulaLoss)+"/";
			    	}
				}
		    }
		    String inchi = "?";
			if(!form.contains("?")){
				inchi = "InChi=1S/"+form.split(",")[0];
				if((int)charge == 1)
					inchi = inchi+"/p+1";
				else if((int)charge == -1)
					inchi = inchi+"/p-1";
					
		    }
			
			String mfP = "?";
			if(ionFinal.getFormulaPath() != null)
				mfP = ionFinal.getFormulaPath();
//			if(form.contains("?"))
//				form = "?";
			String inchiL = "?";
			if(formNL.contains("?") || formNL.equals(""))
				formNL = "?";
			else
				inchiL = "InChi=1S/"+formNL+"/";

			String color = "0";
			if(ionFinal.getMolecule() != null)
				if(ionFinal.getMolecule().getProperty(MZDataConstants.COLOR) != null){
					color = (String)ionFinal.getMolecule().getProperty(MZDataConstants.COLOR);
				}
			
			String colorLo = "0";
			if(ionFinal.getMoleculeLoss() != null){
				if(ionFinal.getMoleculeLoss().getProperty(MZDataConstants.COLOR) != null){
					colorLo = (String)ionFinal.getMoleculeLoss().getProperty(MZDataConstants.COLOR);
				}
				
			}
			
			double massLoss = 0.0;
			if(ionFinal.getParent() != null)
				massLoss = ionFinal.getParent().getMass()-massExp;
//			System.out.println(ionFinal.getGroup()+",:"+ionFinal.getID()+",:"+precursorID+",:"+massExp+",:"+ionFinal.getIntensity_Abs()+",:"+ionFinal.getIntensity_Nor()+",:"+mfP+",:"+nomPath+",:"+level+",:"+form+",:"+formNL+",:"+inchi+",:"+inchiL+",:"+massLoss+",:"+color+",:"+colorLo+">");
		    output += ionFinal.getGroup()+","+ionFinal.getID()+","+precursorID+","+massExp+","+ionFinal.getIntensity_Abs()+","+ionFinal.getIntensity_Nor()+","+mfP+","+nomPath+","+level+","+form+","+formNL+","+inchi+","+inchiL+","+massLoss+","+color+","+colorLo+","+ionFinal.getIntensity_Nor2()+"\n";
		    
		    return output;
		}

	private static ParentIon extractParentIon(IMoleculeSet molecules, List<CMLSpectrum> newCMLSpeList,String scanIDP, MZData mzData) {
		CMLPeak cmlPeak = getCMLPeakFromIDP(newCMLSpeList,scanIDP);
		if(cmlPeak == null)
			return null;
		ParentIon parentIon1 = new ParentIon(cmlPeak.getXValue(),cmlPeak.getYValue(),cmlPeak.getId(), 1);
		
		CMLSpectrum spectrum1 = (CMLSpectrum) cmlPeak.getParent().getParent();
		List<CMLMetadata> cmlL = spectrum1.getMetadataListElements().get(0).getMetadataDescendants();
		Integer numGroup = null;
		for(CMLMetadata metadata : cmlL)
			if(metadata.getDictRef().equals(MZDataConstants.GROUP_PEAK_MSN))
				numGroup = Integer.parseInt(metadata.getContent());
		
		parentIon1.setGroup(numGroup);
		
		parentIon1.setIntensity_Nor(1.0);
		parentIon1.setIntensity_Nor2(1.0);
		if(molecules == null)
			return null;
		parentIon1.addFormulaSet(getFormulas(molecules,cmlPeak));
		IMolecule molecule1 = (IMolecule) getMolecule(molecules,cmlPeak);
		String path = "";
		if(molecule1 != null){
			if(((IMolecule)molecule1).getProperty(CDKConstants.FORMULA) instanceof ArrayList){
				path = ((ArrayList<String>)((IMolecule)molecule1).getProperty(CDKConstants.FORMULA)).get(0).replace(" ", "");
			}else if(((IMolecule)molecule1).getProperty(CDKConstants.FORMULA) instanceof IMolecularFormula){
				path = (MolecularFormulaManipulator.getString((IMolecularFormula)((IMolecule)molecule1).getProperty(CDKConstants.FORMULA))).replace(" ", "");
			}
		}
		parentIon1.addFormulaPath(path);
		parentIon1.setScan(scanIDP);
		parentIon1.addID_XCMS(cmlPeak.getId());
		parentIon1.setMolecule(molecule1);
		HashMap<String,ParentIon> hashMap = new HashMap<String,ParentIon>();
		hashMap.put(""+0, parentIon1);
		List<CMLPeak> peakList = new ArrayList<CMLPeak>();
		for(CMLSpectrum spectrum:newCMLSpeList){
			for(CMLPeak peak:spectrum.getPeakListElements().get(0).getPeakDescendants()){
				peakList.add(peak);
			}
		}
		int count = 1;
//		System.out.println("peakList: "+peakList.size());
		for(int i = 0 ; i < peakList.size(); i++){
			ParentIon parentIon = hashMap.get(""+i);
			if(parentIon == null)
				break;
//			System.out.print(parentIon.getID());
//			System.out.print("("+parentIon.getMass()+")");
			for(CMLSpectrum spectrum2:newCMLSpeList){

				int nLevel = -1;
				String scanNum = "";
				scanNum = spectrum2.getId();
				List<CMLMetadata> ml = spectrum2.getMetadataListElements().get(0).getMetadataDescendants();
				for(CMLMetadata metadata:ml)
					if(metadata.getDictRef().equals(MZDataConstants.MS_LEVEL) )
						nLevel = Integer.parseInt(metadata.getContent());
//					else if(metadata.getDictRef().equals(MZDataConstants.SCAN_NUM) )
//						scanNum =metadata.getContent();
				
				CMLConditionList conditionList = spectrum2.getConditionListElements().get(0);
				boolean precuB = false;
				boolean precuS = false;
				Iterator conditionIterator = conditionList.getScalarElements().iterator();
				while (conditionIterator.hasNext()) {
					CMLScalar condition = (CMLScalar) conditionIterator.next();
					if(condition.getDictRef().equals(MZDataConstants.PRECURSOR_MZ) ){
						double value = Double.parseDouble(condition.getValue());
						if(nLevel == 2)
							value = cmlPeak.getXValue();
						double range = 0.01;// case during nominal masss
						if(molecules.getAtomContainerCount() == 0)
							range = 0.5;
						if(((value + range) > parentIon.getMass()) && ((value - range) < parentIon.getMass()))
							precuB = true;
					}
					if(condition.getDictRef().equals(MZDataConstants.PRECURSOR_SCAN) ){
						String scan = condition.getValue();
						if(scan.equals(parentIon.getScan())){
							precuS = true;
						}
					}
				}
				if(precuB && precuS){
					double highestMass = 0.0;
					for(CMLPeak peak2:spectrum2.getPeakListElements().get(0).getPeakDescendants()){ // extraction of the normalization
						if(peak2.getYValue() > highestMass)
							highestMass = peak2.getYValue();
					}

					List<CMLMetadata> cmlL2 = spectrum2.getMetadataListElements().get(0).getMetadataDescendants();
					Integer numGroup2 = null;
					for(CMLMetadata metadata : cmlL2)
						if(metadata.getDictRef().equals(MZDataConstants.GROUP_PEAK_MSN))
							numGroup2 = Integer.parseInt(metadata.getContent());
					
					for(CMLPeak peak2:spectrum2.getPeakListElements().get(0).getPeakDescendants()){
						ParentIon newIon = new ParentIon(correctingMass(peak2.getXValue(),Polarity.nothing),peak2.getYValue(), ""+peak2.getId(), nLevel);
						newIon.addID_XCMS(peak2.getId());
						newIon.setScan(scanNum);
//						System.out.print(" -> "+spectrum2.getId()+":["+count+"/"+newIon.getID_XCMS()+"]("+newIon.getMass()+")");
						
						newIon.setGroup(numGroup2);

						newIon.addFormulaSet(getFormulas(molecules,peak2));
						newIon.addFormulaLossSet(getFormulasLoss(molecules,peak2));
						IMolecule molecule2 = (IMolecule) getMolecule(molecules,peak2);
						if(molecule2 != null){
							String pathMF = getPath(molecule1,molecule2,mzData.getListReactions());
							newIon.addFormulaPath(pathMF);
						}
						newIon.setMolecule(molecule2);
						newIon.setMoleculeLoss(getMoleculeLoss(molecule2,mzData.getListReactions()));
						parentIon.addFragment(newIon);
			        	newIon.addParent(parentIon);
			        	newIon.setIntensity_Nor(peak2.getYValue()/highestMass);
			        	newIon.setIntensity_Nor2(peak2.getYValue()/parentIon.getIntensity_Abs());
			        	hashMap.put(""+count, newIon);
			        	count ++;
					}
				}
				
			}
//			System.out.println();
		}
		return parentIon1;
	}
	private static IMolecule getMoleculeLoss(IAtomContainer ac,IReactionScheme reactionScheme) {
		boolean flag = false;
		IMolecule moleculeLoss = null;
		IReactionSet reactions = ReactionSchemeManipulator.getAllReactions(reactionScheme);
		for(IReaction reaction:reactions.reactions()){
			IMoleculeSet molecules = ReactionManipulator.getAllProducts(reaction);
			for(IAtomContainer mo:molecules.molecules()){
				if(mo.equals(ac)){
					flag = true;
					break;
				}
			}
			if(flag){
				for(IAtomContainer mo:molecules.molecules()){
					if(!mo.equals(ac)){
						moleculeLoss = (IMolecule) mo;
						break;
					}
				}
				break;
			}
		}
		return moleculeLoss;
	}
	private static ParentIon extractParentIonN(List<CMLSpectrum> newCMLSpeList,String scanIDP, MZData mzData) {
		CMLPeak cmlPeak = getCMLPeakFromIDP(newCMLSpeList,scanIDP);
		if(cmlPeak == null)
			return null;
		ParentIon parentIon1 = new ParentIon(cmlPeak.getXValue(),cmlPeak.getYValue(),cmlPeak.getId(), 1);
		
		CMLSpectrum spectrum1 = (CMLSpectrum) cmlPeak.getParent().getParent();
		List<CMLMetadata> cmlL = spectrum1.getMetadataListElements().get(0).getMetadataDescendants();
		Integer numGroup = null;
		for(CMLMetadata metadata : cmlL)
			if(metadata.getDictRef().equals(MZDataConstants.GROUP_PEAK_MSN))
				numGroup = Integer.parseInt(metadata.getContent());
		
		parentIon1.setGroup(numGroup);
		
		parentIon1.setIntensity_Nor(1.0);
		parentIon1.setIntensity_Nor2(1.0);
		parentIon1.setScan(scanIDP);
		parentIon1.addID_XCMS(cmlPeak.getId());
		parentIon1.addCMLPeak(cmlPeak);
		HashMap<String,ParentIon> hashMap = new HashMap<String,ParentIon>();
		hashMap.put(""+0, parentIon1);
		List<CMLPeak> peakList = new ArrayList<CMLPeak>();
		for(CMLSpectrum spectrum:newCMLSpeList){
			for(CMLPeak peak:spectrum.getPeakListElements().get(0).getPeakDescendants()){
				peakList.add(peak);
			}
		}
		int count = 1;
//		System.out.println("peakList: "+peakList.size());
		for(int i = 0 ; i < peakList.size(); i++){
			ParentIon parentIon = hashMap.get(""+i);
			if(parentIon == null)
				break;
//			System.out.print(parentIon.getID());
//			System.out.print("("+parentIon.getMass()+")");
			for(CMLSpectrum spectrum2:newCMLSpeList){

				int nLevel = -1;
				String scanNum = "";
				scanNum = spectrum2.getId();
				List<CMLMetadata> ml = spectrum2.getMetadataListElements().get(0).getMetadataDescendants();
				for(CMLMetadata metadata:ml)
					if(metadata.getDictRef().equals(MZDataConstants.MS_LEVEL) )
						nLevel = Integer.parseInt(metadata.getContent());
//					else if(metadata.getDictRef().equals(MZDataConstants.SCAN_NUM) )
//						scanNum =metadata.getContent();
				
				CMLConditionList conditionList = spectrum2.getConditionListElements().get(0);
				boolean precuB = false;
				boolean precuS = false;
				Iterator conditionIterator = conditionList.getScalarElements().iterator();
				while (conditionIterator.hasNext()) {
					CMLScalar condition = (CMLScalar) conditionIterator.next();
					if(condition.getDictRef().equals(MZDataConstants.PRECURSOR_MZ) ){
						double value = Double.parseDouble(condition.getValue());
						if(nLevel == 2)
							value = cmlPeak.getXValue();
						
						double range = 0.01;// case during nominal masss
//						if(molecules.getAtomContainerCount() == 0)
							range = 0.5;
						if(((value + range) > parentIon.getMass()) && ((value - range) < parentIon.getMass()))
							precuB = true;
					}
					if(condition.getDictRef().equals(MZDataConstants.PRECURSOR_SCAN) ){
						String scan = condition.getValue();
						if(scan.equals(parentIon.getScan())){
							precuS = true;
						}
					}
				}
				if(precuB && precuS){
					double highestMass = 0.0;
					for(CMLPeak peak2:spectrum2.getPeakListElements().get(0).getPeakDescendants()){ // extraction of the normalization
						if(peak2.getYValue() > highestMass)
							highestMass = peak2.getYValue();
					}

					List<CMLMetadata> cmlL2 = spectrum2.getMetadataListElements().get(0).getMetadataDescendants();
					Integer numGroup2 = null;
					for(CMLMetadata metadata : cmlL2)
						if(metadata.getDictRef().equals(MZDataConstants.GROUP_PEAK_MSN))
							numGroup2 = Integer.parseInt(metadata.getContent());
					
					for(CMLPeak peak2:spectrum2.getPeakListElements().get(0).getPeakDescendants()){
						ParentIon newIon = new ParentIon(correctingMass(peak2.getXValue(),Polarity.nothing),peak2.getYValue(), ""+peak2.getId(), nLevel);
						newIon.addID_XCMS(peak2.getId());
						newIon.setScan(scanNum);
//						System.out.print(" -> "+spectrum2.getId()+":["+count+"/"+newIon.getID_XCMS()+"]("+newIon.getMass()+")");
						
						newIon.setGroup(numGroup2);
						newIon.addCMLPeak(peak2);

						parentIon.addFragment(newIon);
			        	newIon.addParent(parentIon);
			        	newIon.setIntensity_Nor(peak2.getYValue()/highestMass);
			        	newIon.setIntensity_Nor2(peak2.getYValue()/parentIon.getIntensity_Abs());
			        	hashMap.put(""+count, newIon);
			        	count ++;
					}
				}
				
			}
//			System.out.println();
		}
		return parentIon1;
	}
	private static IMolecularFormulaSet getFormulas(IMoleculeSet molecules, CMLPeak cmlPeak) {
		if(cmlPeak.getMoleculeRefs() == null)
			return null;
		String ref = cmlPeak.getMoleculeRefs()[0];
		for(IAtomContainer molecule : molecules.molecules()){
			if(molecule.getID().equals(ref)){
				IMolecularFormulaSet formSet = extractMFSet(molecule);
				return formSet;
			}
		}
		return null;
	}
	private static IAtomContainer getMolecule(IMoleculeSet molecules, CMLPeak cmlPeak) {
		if(cmlPeak.getMoleculeRefs() == null)
			return null;
		String ref = cmlPeak.getMoleculeRefs()[0];
		for(IAtomContainer molecule : molecules.molecules()){
			if(molecule.getID().equals(ref)){
				return molecule;
			}
		}
		return null;
	}
	private static IMolecularFormulaSet getFormulasLoss(IMoleculeSet molecules, CMLPeak cmlPeak) {
		if(cmlPeak.getMoleculeRefs() == null)
			return null;
		String ref = cmlPeak.getMoleculeRefs()[0];
		for(IAtomContainer molecule : molecules.molecules()){
			if(molecule.getID().equals(ref+"Loss")){
				IMolecularFormulaSet formSet = extractMFSet(molecule);
				return formSet;
			}
		}
		return null;
	}
	private static IMolecularFormulaSet extractMFSet(IAtomContainer molecule) {
		double charge = Double.valueOf((String) molecule.getProperty("cdk:partialCharge"));
		IMolecularFormulaSet formSet = molecule.getBuilder().newMolecularFormulaSet();
		if(molecule.getProperty(CDKConstants.FORMULA) instanceof ArrayList){
    		ArrayList<String> mfL = ((ArrayList<String>)((IMolecule)molecule).getProperty(CDKConstants.FORMULA));
    		for(String formula: mfL){
    			IMolecularFormula form = MolecularFormulaManipulator.getMolecularFormula(formula, molecule.getBuilder());
    			form.setCharge((int)charge);
    			formSet.addMolecularFormula(form);
    		}
		}else if(((IMolecule)molecule).getProperty(CDKConstants.FORMULA) instanceof IMolecularFormula){
			IMolecularFormula formula = (IMolecularFormula)((IMolecule)molecule).getProperty(CDKConstants.FORMULA);
			formula.setCharge((int)charge);
			formSet.addMolecularFormula(formula);
		}else if(((IMolecule)molecule).getProperty(CDKConstants.FORMULA) instanceof IMolecularFormulaSet){
			IMolecularFormulaSet formulaSet = (IMolecularFormulaSet)((IMolecule)molecule).getProperty(CDKConstants.FORMULA);
			for(IMolecularFormula form : formulaSet.molecularFormulas()){
				form.setCharge((int)charge);
				formSet.addMolecularFormula(form);
			}
		}
		return formSet;
	}
	/**
	 * Correct the Mass according the  Polarity. 
	 * Polarity.positive designed for positive mode: To compare
	 * with the Theoretical EF you need to add the mass of one electron.
	 * For negative mode the other way around. 
	 * 
	 * @param hashMap  HashMap with the masses 
	 * @param polarity The Polarity
	 * @return         The correct value mass
	 */
	private static double correctingMass(double mass, 
			Polarity polarity) {
		/** Charge of the electron. */
		double chargeE = 0.00054857990927;
		
		if(polarity.equals(Polarity.positive))
			return (mass+chargeE);
		else if(polarity.equals(Polarity.negative))
			return (mass-chargeE);
		else
			return mass;

	}
	/**
	 * Extract the cmlPeak containing a moleculeRefs
	 * 
	 * @param newCMLSpeList A list of CMLSpectrum objects
	 * @param scanIDP  ID of the spectrum
	 * @return The CMLPeak
	 */
	private static CMLPeak getCMLPeakFromIDP(List<CMLSpectrum> newCMLSpeList,String scanIDP) {
		// TODO : windonw at the moment set to 1
		for(CMLSpectrum spectrum:newCMLSpeList){
			List<CMLMetadata> ml = spectrum.getMetadataListElements().get(0).getMetadataDescendants();
			for(CMLMetadata metadata:ml){
				if(metadata.getDictRef().equals(MZDataConstants.MS_LEVEL) && metadata.getContent().equals("1")){
					if(spectrum.getId().equals(scanIDP)){
						for(CMLPeak peak:CMLSpectrum.getDescendantPeaks(spectrum)){
							if(peak.getMoleculeRefs() != null){
								return peak;
							}
						}
					}
				}
			}
		}
		// when it is null take the highest peak of MS1
		CMLPeak cmlPeak = null;
		for(CMLSpectrum spectrum:newCMLSpeList){
			List<CMLMetadata> ml = spectrum.getMetadataListElements().get(0).getMetadataDescendants();
			for(CMLMetadata metadata:ml){
				if(metadata.getDictRef().equals(MZDataConstants.MS_LEVEL) && metadata.getContent().equals("1")){
					if(spectrum.getId().equals(scanIDP)){
						for(CMLPeak peak:CMLSpectrum.getDescendantPeaks(spectrum)){
							if(cmlPeak != null){
								if(peak.getYValue() > cmlPeak.getYValue())
									cmlPeak = peak;
							}else
								cmlPeak = peak;
						}
					}
				}
			}
		}
		return cmlPeak;
	}

	private static IsotopePattern getIsotopePatternFromIDP(List<CMLSpectrum> cmlSpectList, CMLPeak peakID) {
		if(peakID == null)
			return null;
		CMLSpectrum spectrum = CMLSpectrum.getSpectrum(peakID);
		List<CMLPeak> peakList = new ArrayList<CMLPeak>();
		for(CMLPeak peak:spectrum.getPeakListElements().get(0).getPeakDescendants()){
			if((peak.getXValue() > peakID.getXValue())  && (peakID.getXValue() + 3.5 > peak.getXValue()))
				if(peakID.getYValue() > peak.getYValue()){
					peakList.add(peak);
				}
		}
		
		IsotopePattern spExp = new IsotopePattern();

		double massM = peakID.getXValue();
		double inteM = peakID.getYValue();
		spExp.setMonoIsotope(new IsotopeContainer(massM, inteM/inteM));
		for(CMLPeak peaki:peakList){
			spExp.addIsotope(new IsotopeContainer(peaki.getXValue(), peaki.getYValue()/inteM));
		}
		return spExp;
	}
	/**
	 * I take an employee element and read the values in, create
	 * an Employee object and return it
	 */
	private ScanMZXML getScanMZXML(Element scanEl) {

		//for each <employee> element get text or int values of
		//name ,id, age and name
		String num = scanEl.getAttribute("num");
		String pol = getTextValue(scanEl,"polarity");
		String collisionEnergy = getTextValue(scanEl,"collisionEnergy");

		//Create a new Employee with the value read from the xml nodes
		ScanMZXML e = new ScanMZXML(num,pol,collisionEnergy);

		return e;
	}


	/**
	 * I take a xml element and the tag name, look for the tag and get
	 * the text content
	 * i.e for <employee><name>John</name></employee> xml snippet if
	 * the Element points to employee node and tagName is 'name' I will return John
	 */
	private String getTextValue(Element ele, String tagName) {
		String textVal = null;
		NodeList nl = ele.getElementsByTagName(tagName);
		if(nl != null && nl.getLength() > 0) {
			Element el = (Element)nl.item(0);
			textVal = el.getFirstChild().getNodeValue();
		}

		return textVal;
	}
	
	/**
	 * Return a MZData object with the theoretical mass of the elemental composition.
	 * 
	 * @param mzData The MZData to change
	 * @return       The MZData with the theoretical masses
	 */
	public static MZData getTheoretical(MZData mzData){
		try {
			MZData newMZData = mzData.clone();
			CMLSpectrumList cmlSpectList = newMZData.getListSpectra();
			IReactionScheme cmlReacList = newMZData.getListReactions();
			IMoleculeSet molecules = ReactionSchemeManipulator.getAllMolecules(cmlReacList);
			double charge = Double.valueOf((String) molecules.getMolecule(0).getProperty("cdk:partialCharge"));
			CMLElements<CMLSpectrum> specElem = cmlSpectList.getSpectrumElements();
			for(CMLSpectrum spectrum:specElem){
				for(CMLPeak peak:CMLSpectrum.getDescendantPeaks(spectrum)){
					peak.setXValue(10.0);
					if(peak.getMoleculeRefs() != null){
						String ref = peak.getMoleculeRefs()[0];
						for(IAtomContainer molecule : molecules.molecules()){
							if(molecule.getID().equals(ref)){
								IMolecularFormulaSet formSet = extractMFSet(molecule);
								if(formSet.size() != 0){
									IMolecularFormula fm = formSet.getMolecularFormula(0);
									IMolecularFormula formula = MolecularFormulaManipulator.getMajorIsotopeMolecularFormula(MolecularFormulaManipulator.getString(fm),builder);
								    formula.setCharge((int) charge);
								    double mass2 = MolecularFormulaManipulator.getTotalExactMass(formula);
								    peak.setXValue(mass2);
								}
							}
						}
					}
				}
			}
			return newMZData;

		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
	}
	/**
	 * Colored a MZData based on what nodes are present from other MZData.
	 *  
	 * @param mzDataA  The 1 MZData object
	 * @param mzDataB  The 2 MZData object
	 * @return         The MZData of the object A including the color property
	 */
	public static MZData getColored(MZData mzDataA,MZData mzDataB) {
		FingerMZData fpMZD = new FingerMZData();
		List<String> fp_AL0 = fpMZD.getFingerprintL0(mzDataA);
//		List<String> fp_AL1 = fpMZD.getFingerprintL1(mzDataA);
		List<String> fp_AL2 = fpMZD.getFingerprintL2(mzDataA);
		List<String> fp_AL3 = fpMZD.getFingerprintL3(mzDataA);
		List<String> fp_AL4 = fpMZD.getFingerprintL4(mzDataA);
		List<String> fp_AL5 = fpMZD.getFingerprintL5(mzDataA);
		Collections.sort(fp_AL0);
		

		List<String> fp_AlL0 = fpMZD.getFingerprintL0_loss(mzDataA);
//		List<String> fp_AlL1 = fpMZD.getFingerprintL1_Loss(mzDataA);
		List<String> fp_AlL2 = fpMZD.getFingerprintL2_Loss(mzDataA);
		List<String> fp_AlL3 = fpMZD.getFingerprintL3_Loss(mzDataA);
		List<String> fp_AlL4 = fpMZD.getFingerprintL4_Loss(mzDataA);
		List<String> fp_AlL5 = fpMZD.getFingerprintL5_Loss(mzDataA);
		Collections.sort(fp_AlL0);
		
		List<MZData> mzDataList = new ArrayList<MZData>();
		mzDataList.add(mzDataA);
		
		Map<String, List<IAtomContainer>> hashMap = getHashMapPath(mzDataList);
		Map<String, List<IAtomContainer>> hashMapL = getHashMapPathLoss(mzDataList.get(0).getListReactions());
		
		List<String> fp_BL0 = fpMZD.getFingerprintL0(mzDataB);
		List<String> fp_BL1 = fpMZD.getFingerprintL1(mzDataB);
		Collections.sort(fp_BL0);
		
		List<String> fp_BlL0 = fpMZD.getFingerprintL0_loss(mzDataB);
		List<String> fp_BlL1 = fpMZD.getFingerprintL1_Loss(mzDataB);
		Collections.sort(fp_BlL0);
		
		if(fp_BL0.size() == 0)
			return mzDataA;
		
		for(String fp_A : fp_AL0){
			String countC = "0";
			String formula = fp_A.substring(fp_A.lastIndexOf("|")+1, fp_A.length());

			if(fp_BL1.contains(formula)){
				countC = getllP(formula,fp_BL0, fp_A, fp_AL2, fp_AL3, fp_AL4, fp_AL5);
    		}
			IAtomContainer ac = hashMap.get(fp_A).get(0);
			ac.setProperty(MZDataConstants.COLOR, countC);
			
		}
		
		for(String fp_A : fp_AlL0){
			String countC = "0";
			String formula = fp_A.substring(fp_A.lastIndexOf("|")+1, fp_A.length());
			if(fp_BlL1.contains(formula)){
				countC = getllP(formula,fp_BlL0, fp_A, fp_AlL2, fp_AlL3, fp_AlL4, fp_AlL5);
    		}
			IAtomContainer ac = hashMapL.get(fp_A).get(0);
			IAtomContainer acl = getMoleculeLoss(ac, mzDataA.getListReactions());
			acl.setProperty(MZDataConstants.COLOR, countC);
//			System.out.println(ac.getID()+">"+acl.getID()+" "+formula+" "+fp_A+" "+countC);
		}
		return mzDataA;
	}
	
	
	private static String getllP(String mf, List<String>fpF2, String fp, List<String> fpL12, List<String> fpL13, List<String> fpL14, List<String> fpL15) {

		boolean fp2B = false;
		boolean fp3B = false;
		boolean fp4B = false;
		boolean fp5B = false;
		int lim = fp.split("\\|\\|").length;
//		System.out.println("lim:"+lim);
//		System.out.println("fpF2:"+fpF2);// all paths
		for(String nmptt : fpF2){
//			if(!nmptt.endsWith(fp))
//				continue;
//			System.out.println("nmptt:"+nmptt);
			List<String> ff = getFingerprintsL(nmptt,2);
			List<String> fff = getFingerprintsL(fp,2);
			if(lim > 1)
				for(String f : ff){
					if(fff.contains(f)){
						String mf2 = f.substring(f.lastIndexOf("@")+1, f.length());
						if(mf2.equals(mf)){
//							System.out.println("f2:"+f);
//							System.out.println("mf2:"+mf2);
							fp2B = true;
							break;
						}
					}
				}
			ff = getFingerprintsL(nmptt,3);
			fff = getFingerprintsL(fp,3);
			if(fp2B && lim > 2)
				for(String f : ff){//B fp3
					if(fff.contains(f)){
						String mf2 = f.substring(f.lastIndexOf("@")+1, f.length());
						if(mf2.equals(mf)){
//							System.out.println("f3:"+f);
//							System.out.println("mf3:"+mf2);
//							System.out.println("fpL13:"+fpL13);
							fp3B = true;
							break;
						}
					}
				}
			ff = getFingerprintsL(nmptt,4);
			fff = getFingerprintsL(fp,4);
			if(fp3B && lim > 3)
				for(String f : ff){
					if(fff.contains(f)){
						String mf2 = f.substring(f.lastIndexOf("@")+1, f.length());
						if(mf2.equals(mf)){
							fp4B = true;
							break;
						}
					}
				}
			ff = getFingerprintsL(nmptt,5);
			fff = getFingerprintsL(fp,5);
			if(fp4B && lim > 4)
				for(String f : ff){
					if(fff.contains(f)){
						String mf2 = f.substring(f.lastIndexOf("@")+1, f.length());
						if(mf2.equals(mf)){
							fp5B = true;
							break;
						}
					}
				}
		}
//		System.out.println("fp2B:"+fp2B+", fp3B:"+fp3B+", fp4B:"+fp4B+", fp5B:"+fp5B);
		String llP = "1";
		if(fp2B == true){
			llP = "2";
			if(fp3B == true){
				llP = "3";
				if(fp4B == true){
					llP = "4";
					if(fp5B == true){
						llP = "5";
					}
				}
			}
		}
		return llP;
	}


	private static List<String> getFingerprintsL(String path, int sign) {
		List<String> fingerList = new ArrayList<String>();
		String[] fm = path.split("\\|\\|");
			for(int i = 0 ; i < fm.length; i++){
				if(fm.length >= i+sign){
					String finger = "";
					for(int j = i ; j < i+sign; j++){
						if(i == j)
							finger += fm[j];
						else
							finger += "@"+fm[j];
						
					}
					if(!fingerList.contains(finger)){
						fingerList.add(finger);
					}
				}
			}
		return fingerList;
	}
	public static MZData SetInChIs(MZData mzData, HashMap<String, String> inchiMAP) {
		IMoleculeSet molecules = ReactionSchemeManipulator.getAllMolecules(mzData.getListReactions());
		
		for(IAtomContainer molecule : molecules.molecules()){
			if(inchiMAP.containsKey(molecule.getID())){
	    		molecule.setProperty(CDKConstants.INCHI, inchiMAP.get(molecule.getID()));
	    	}
		}
		return mzData;
	}

	public static double[][] getMZDataMatrixFromFile(String filePath) throws IOException {
		String line;
		int count = 0 ;
		int colnum = 0;
		BufferedReader br = new BufferedReader(new FileReader(filePath));
		while ((line = br.readLine()) != null) {
	    	count ++;
	    	String[] splits = line.split(",");
	    	colnum = splits.length;
	    }
	    br.close();
	    
	    if(colnum == 5){
		    double[][] matrix = new double[count][5];
			//0:ID, 1:IDparent, 2:m/z, 3:int,4:charge
			count = 0 ;
		    br = new BufferedReader(new FileReader(filePath));
			while ((line = br.readLine()) != null) {
		    	String[] splits = line.split(",");
		    	matrix[count][0] = Double.parseDouble(splits[0]);  //ID
		    	matrix[count][1] = Double.parseDouble(splits[1]);  //IDparent  
		    	matrix[count][2] = Double.parseDouble(splits[2]);  //mass
		    	matrix[count][3] = Double.parseDouble(splits[3]);  //int
		    	matrix[count][4] = Double.parseDouble(splits[4]);  //charge
		    	count ++;
			}
		    br.close();
		    
			return matrix;
	    }if(colnum == 9){
	    	double[][] matrix = new double[count][8];
			//0:ID, 1:IDparent, 2:m/z, 3:int,4:charge
			count = 0 ;
		    br = new BufferedReader(new FileReader(filePath));
			while ((line = br.readLine()) != null) {
		    	String[] splits = line.split(",");
		    	matrix[count][0] = Double.parseDouble(splits[0]);  //ID
		    	matrix[count][1] = Double.parseDouble(splits[1]);  //IDparent  
		    	matrix[count][2] = Double.parseDouble(splits[2]);  //level 
		    	matrix[count][3] = Double.parseDouble(splits[3]);  //rt
		    	matrix[count][4] = Double.parseDouble(splits[4]);  //mass
		    	matrix[count][5] = Double.parseDouble(splits[5]);  //int
		    	matrix[count][6] = Double.parseDouble(splits[6]);  //file
		    	matrix[count][7] = Double.parseDouble(splits[7]);  //group
		    	count ++;
			}
		    br.close();
		    
			return matrix;
	    }else
	    	return null;
	}
	/**
	 * Split a mzData containing more of 1 groups in a List of mzData
	 * 
	 * @param mzData
	 * @return
	 */
	public static List<MZData> splitMZData(MZData mzData) {
		List<MZData> mzDataList = new ArrayList<MZData>();
		
		Map<Object, Object> properties = mzData.getProperties();
		int numGroups = Integer.parseInt((String)properties.get(MZDataConstants.NUM_GROUPS));
		if(numGroups == 1){
			mzDataList.add(mzData);
		}else if(numGroups > 1){
			
			for(int i = 1; i < numGroups+1; i++){
				CMLElements<CMLSpectrum> specElem = mzData.getListSpectra().getSpectrumElements();
				MZData mzs = new MZData();
		        CMLSpectrumList cmlSpecList = new CMLSpectrumList();
		        mzs.setListSpectra(cmlSpecList);
				CMLMetadataList metadataList = new CMLMetadataList();
		        cmlSpecList.addMetadataList(metadataList);
				
				// collecting all spectra from the group
				List<String> molsID = new ArrayList<String>();
				for(CMLSpectrum spectrum: specElem){
					List<CMLMetadata> metadataDescendants = spectrum.getMetadataListElements().get(0).getMetadataDescendants();
					for(CMLMetadata metadata : metadataDescendants){
						if(metadata.getDictRef().equals(MZDataConstants.GROUP_PEAK_MSN)){
							int groupScan = Integer.parseInt((String)metadata.getContent());
							if(groupScan == i){
								List<CMLMetadata> ml = spectrum.getMetadataListElements().get(0).getMetadataDescendants();
								for(CMLMetadata md:ml)
									if(md.getDictRef().equals(MZDataConstants.GROUP_PEAK_MSN)){
										md.setContent("1");
										break;
									}
									
								cmlSpecList.addSpectrum(spectrum);
								// extract the molecules
								List<CMLPeak> peakList = CMLSpectrum.getDescendantPeaks(spectrum);
								for(CMLPeak peak : peakList){
									if(peak.getMoleculeRefs() != null){
										molsID.add(peak.getMoleculeRefs()[0]);
									}
								}
							}
						}
					}
				}
				// Add metadata
				int numSpectra = cmlSpecList.getChildCount();
				List<CMLMetadata> cmlMDList = mzData.getListSpectra().getMetadataListElements().get(0).getMetadataDescendants();
				for(CMLMetadata metadata:cmlMDList){
		        	if(metadata.getDictRef().equals(MZDataConstants.NUM_GROUPS))
		        		metadata.setContent("1");
	        		else if(metadata.getDictRef().equals(MZDataConstants.SCANS))
		        		metadata.setContent(Integer.toString(numSpectra));
        			
		        	metadataList.appendChild(metadata);
				}
				Map<Object, Object> properties1 = mzData.getProperties();
				mzs.setProperties(properties1);
				mzs.setProperty(MZDataConstants.NUM_GROUPS, "1");
				mzs.setProperty(MZDataConstants.SCANS,Integer.toString(numSpectra));
		        
				// collecting all the molecules from the group
		        IReactionScheme reaSche = mzData.getListReactions();
    	    	IReactionSet reaSet = ReactionSchemeManipulator.getAllReactions(reaSche);
    			List<String> reacRemove = new ArrayList<String>();
    			for(IReaction reaction : reaSet.reactions()){
    				IMoleculeSet molecules = ReactionManipulator.getAllMolecules(reaction);
					for(int j = 0; j < molecules.getMoleculeCount() ; j++){
    					String idTo = molecules.getAtomContainer(j).getID();
        				if(!molsID.contains(idTo) && !idTo.contains("Loss")){
							reacRemove.add(reaction.getID());
        					break;
        				}
    				}
    			}
    			if(reacRemove.size() > 0){
    				IReactionSet reaSet2 = reaSet.getBuilder().newReactionSet();
					for(IReaction reaction : reaSet.reactions()){
						if(!reacRemove.contains(reaction.getID()))
								reaSet2.addReaction(reaction);
					}

					IReactionScheme reaSche2 = ReactionSchemeManipulator.createReactionScheme(reaSet2);
			        mzs.setListReactions(reaSche2);
    			}else{
			        mzs.setListReactions(reaSche);
    			}
		        mzDataList.add(mzs);
		        
			}
		}
		
		return mzDataList;
	}
	/**
	 * Get the retention time from a single mzData containing only a fragmentation tree.
	 * 
	 * @param mzData
	 * @return
	 */
	public static double getRt(MZData mzData) {
		double rt = 0.0;
		Map<Object, Object> properties = mzData.getProperties();
		int numGroups = Integer.parseInt((String)properties.get(MZDataConstants.NUM_GROUPS));
		if(numGroups == 1){
			double rtSmall = 10000000.0;
			double rtBig = 0.0;
			CMLElements<CMLSpectrum> specElem = mzData.getListSpectra().getSpectrumElements();
			for(CMLSpectrum spectrum: specElem){
				List<CMLMetadata> metadataDescendants = spectrum.getMetadataListElements().get(0).getMetadataDescendants();
				for(CMLMetadata metadata : metadataDescendants){
					if(metadata.getDictRef().equals(MZDataConstants.RT)){
						double rtPunt = Double.parseDouble((String)metadata.getContent());
						if(rtPunt < rtSmall)
							rtSmall = rtPunt;
						if(rtPunt > rtBig)
							rtBig = rtPunt;
						
					}
				}
			}
			rt = rtBig - rtSmall;
		}
		
		return rt;
	}
}
