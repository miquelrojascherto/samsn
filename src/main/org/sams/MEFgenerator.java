package org.sams;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.ReactionScheme;
import org.openscience.cdk.config.IsotopeFactory;
import org.openscience.cdk.formula.IsotopeContainer;
import org.openscience.cdk.formula.IsotopePattern;
import org.openscience.cdk.formula.MolecularFormulaRange;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IIsotope;
import org.openscience.cdk.interfaces.IMolecularFormulaSet;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IMoleculeSet;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.IReactionScheme;
import org.openscience.cdk.tools.LoggingTool;
import org.openscience.cdk.tools.manipulator.ReactionSchemeManipulator;
import org.sams.manipulator.PrintTools;
import org.sams.spect.Accuracy;
import org.sams.spect.ParentIon;
import org.xmlcml.cml.base.CMLElements;
import org.xmlcml.cml.element.CMLConditionList;
import org.xmlcml.cml.element.CMLMetadata;
import org.xmlcml.cml.element.CMLPeak;
import org.xmlcml.cml.element.CMLScalar;
import org.xmlcml.cml.element.CMLSpectrum;
import org.xmlcml.cml.element.CMLSpectrumList;

/**
 * Class that process spectral tree data.
 *  
 * @author Miguel Rojas-Cherto
 */
public class MEFgenerator {
	/** Bilder to initiate CDK objects */
	static DefaultChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
	
	private MZData mzData;
	public static enum Polarity {positive,negative,nothing}
	private List<Integer> ppm = Arrays.asList(3,3,3,3,3,3,3);
	private Polarity polarity;
	private int countMol = 1;
	private int countScheme = 1;
	private int countReaction = 1;
	private MolecularFormulaRange range;
	private String[] rules = {"nitrogenR","RDBER"};
	private Double massParent;
	private double rangeParent = 0.01;
	private int toleranceISO;
	private double minAbISO;
	private double scoreISO;
	private double rtStart = -1.0;
	private double rtEnd = -1.0;
	private int level = 7;

	private LoggingTool logger;
	

	public MEFgenerator(MZData mzData){
		this.mzData = mzData;
		range = creatingDefaultRange();
	}
	
	/**
	 * Set the mass accuracy in ppm as a List of integers
	 * 
	 * @param ppmList  A List of integers
	 */
	public void setMassAccuracy(List<Integer>  ppmList){
		this.ppm = ppmList;
	}
	
	/**
	 * Set the range of the elemental formulas
	 * 
	 * @param range The MolecularFormulaRange containing the range
	 */
	public void setRangeMolecularFormula(MolecularFormulaRange range) {
		this.range = range;
	}
	
	/**
	 * Fix the mass of the top parent ion to look for.
	 * 
	 * @param mass  The mass of the top ion
	 * @param range The range to look for
	 */
	public void addMassParent(double mass, double range){
		this.massParent = mass;
		this.rangeParent  = range;
	}

	/**
	 * Set the rules
	 * @param rules  An array with the rules
	 */
	public void setRules(String[] rules) {
		this.rules = rules;
	}
	/**
	 * Set the level deep maximal to analyze the fragmentation tree
	 * 
	 * @param level  The level deep
	 */
	public void setLevelMax(int level){
		this.level = level;
	}

	/**
	 * @param toleranceISO Tolerance in ppm scanning the signals
	 * @param scoreISO     Score comparing isotope patterns 
	 * @param minAb        Minimal abundance of the isotopes to be added 
	 * 				in the combinatorial search
	 */
	public void setIsopePatternParameters(double minAbISO, int toleranceISO, double scoreISO) {
		this.toleranceISO = toleranceISO;
		this.minAbISO = minAbISO;
		this.scoreISO = scoreISO;
	}

	/**
	 * Set the retention time range to look for
	 * 
	 * @param rtStart The start
	 * @param rtEnd   The end
	 */
	public void setRT(double rtStart, double rtEnd) {
		this.rtStart = rtStart;
		this.rtEnd = rtEnd;
	}
	/**
	 * Initiate the processing
	 * 
	 * @return The MZData object containing all the information.
	 */
	public MZData init(){
		MZData newMZData = new MZData();
		CMLSpectrumList cmlSpectList = mzData.getListSpectra();
		Map<Object, Object> properties = mzData.getProperties();
//		CMLMetadataList cp_metadataList = cmlSpectList.get;
		
		//TODO substitute by loggingTool
//		logger = new LoggingTool(this);
//		PrintTools printTools = new PrintTools("/tmp/out-tree-process.txt");
		
		
		// Extraction of the number of groups. It defines the number of spectral trees.
		int numGroups = 0;
		List<CMLMetadata> metadataList = cmlSpectList.getMetadataListElements().get(0).getMetadataDescendants();

        for(CMLMetadata metadata:metadataList){
			if(metadata.getDictRef().equals(MZDataConstants.NUM_GROUPS))
				numGroups = Integer.parseInt(metadata.getContent());
		}
        
		//////////////////////////////////////////////////////ACCURACY
		CMLMetadata metadataL = new CMLMetadata();
		metadataL.setDictRef(MZDataConstants.ACCURACY);
		String accS = "";
		for(int i = 0; i < ppm.size() ; i++){
			accS += ppm.get(i);
			if(i != ppm.size()-1)
				accS += ",";
		}
		metadataL.setContent(accS);
		cmlSpectList.getMetadataListElements().get(0).addMetadata(metadataL);
		properties.put(MZDataConstants.ACCURACY, accS);
		//////////////////////////////////////////////////////RULES
		metadataL = new CMLMetadata();
		metadataL.setDictRef(MZDataConstants.RULES);
		String rulesS = "";
		for(int i = 0; i < rules.length ; i++){
			rulesS += rules[i];
			if(i != rules.length-1)
				rulesS += ",";
		}
		metadataL.setContent(rulesS);
		cmlSpectList.getMetadataListElements().get(0).addMetadata(metadataL);
		properties.put(MZDataConstants.RULES, rulesS);
		//////////////////////////////////////////////////////ELEMENTS
		metadataL = new CMLMetadata();
		metadataL.setDictRef(MZDataConstants.ELEMENTS);
		String rangeL = getRangeLine(range);
		metadataL.setContent(rangeL);
		cmlSpectList.getMetadataListElements().get(0).addMetadata(metadataL);
		properties.put(MZDataConstants.ELEMENTS, rangeL);
		
		
		// Each spectral tree will be set in an own List of CMLSpectrum
		List<List<CMLSpectrum>> cmlSpectListList = new ArrayList<List<CMLSpectrum>>();
		// Create a list with all parent ion top 
		List<ParentIon> listPI = new ArrayList<ParentIon>();
		CMLElements<CMLSpectrum> specElem = cmlSpectList.getSpectrumElements();
		for(int i = 0 ; i < numGroups; i++){
			List<CMLSpectrum> newCMLSpeList = new ArrayList<CMLSpectrum>();
			cmlSpectListList.add(newCMLSpeList);
			
			String peakIDP = null; // Peak ID of the top 
			String scanNum = "";   // scan number belonging this peak ID
			boolean inRT = false;
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
					// when we are interested in specify the rt to look for
					if((rtStart == -1.0) && (rtEnd == -1.0))
						inRT = true;
					else
						for(CMLMetadata metadata:ml){
							if(metadata.getDictRef().equals(MZDataConstants.RT) ){
								double valueRT = Double.parseDouble(metadata.getContent());
								if((valueRT >= rtStart*60) && (valueRT <= rtEnd*60)){
									inRT = true ;
								}
							}
						}
				}
				if(validS){
					// determine the peakID, polarity and scanNum
					for(CMLMetadata metadata:ml){
						if(metadata.getDictRef().equals(MZDataConstants.MS_LEVEL) && metadata.getContent().equals("2")){
							CMLConditionList conditionList = spectrum.getConditionListElements().get(0);
							Iterator<CMLScalar> conditionIterator = conditionList.getScalarElements().iterator();
							while (conditionIterator.hasNext()) {
								CMLScalar condition = conditionIterator.next();
								if(condition.getDictRef().equals(MZDataConstants.PRECURSOR_MZ) ){
									peakIDP = condition.getValue();
//									System.out.println("Precursor selected: "+peakIDP);
								}else if(condition.getDictRef().equals(MZDataConstants.POLARITY) ){
									String polarityTMP = condition.getValue();
									if(polarityTMP.equals("positive"))
										polarity = Polarity.positive;
									else if(polarityTMP.equals("negative"))
										polarity = Polarity.negative;
									else 
										polarity = Polarity.nothing;
								}else if(condition.getDictRef().equals(MZDataConstants.PRECURSOR_SCAN) ){
									scanNum = condition.getValue();
								}
							}
						}
					}
				}
			}
			/////////////////////////////////////////////////////////
			if(peakIDP == null)
				continue;
			
			if(!inRT)
				continue;
			
			// create the parentIon from the parent ion Top
	        ParentIon parentIon = extractParentIon(newCMLSpeList,peakIDP,scanNum);
	        
	        if(parentIon == null)
	        	continue;
			
//	        List<Object> param = new ArrayList<Object>();
			Accuracy acc = new Accuracy();
			acc.setAccuracyPPM(parentIon.getMass(), ppm.get(0));
//	        param.add(acc);
//	        param.add(range);
//	        param.add(polarity);
//	        param.add(rules);
	        
//	        GeneralRun run = new GeneralRun();
//	        run.setParameters(param);
	        
	        IsotopePattern isoP = null;
	        for(int r = 0; r < rules.length; r++){
        		if(rules[r].equals("isotopePatternR")){
        			isoP = getIsotopePatternFromIDP(newCMLSpeList,getCMLPeakFromIDP(newCMLSpeList, peakIDP));
        		}
	        }
	        int cyclesMax = 20;
			Generator2TreeMF generator = new Generator2TreeMF(parentIon,range,acc.getAccuracyDa(),level,cyclesMax,polarity,isoP,rules,minAbISO,toleranceISO,scoreISO);

	        
//			ParentIon ionFinal = run.initiate(parentIon, level, isoP,minAbISO,toleranceISO,scoreISO);
	        
	        // printing in /tmp/out-tree-process.txt the logfile with all the process
//			printTools.printThis(generator.getOutprint());
//			printTools.printOutPutCompress(parentIon,polarity);
			
	        addRefMol(newCMLSpeList,parentIon);
	        listPI.add(parentIon);

		}
        IReactionScheme scheme = convert2ReactionScheme(listPI);
//        printTools.close();
        
        // add properties spectrum to the molecules
        relateBothObjects(cmlSpectList,scheme);

		newMZData.setListSpectra(cmlSpectList);
        newMZData.setListReactions(scheme);
        newMZData.setProperties(properties);

		return newMZData;
	}

	private void relateBothObjects(CMLSpectrumList cmlSpectList, IReactionScheme scheme) {
		IMoleculeSet molSet = ReactionSchemeManipulator.getAllMolecules(scheme);
		for(IAtomContainer molecule: molSet.molecules()){
			String idMol = molecule.getID();
			CMLElements<CMLSpectrum> specElem = mzData.getListSpectra().getSpectrumElements();
			for(CMLSpectrum spectrum: specElem){
				List<CMLPeak> peakList = CMLSpectrum.getDescendantPeaks(spectrum);
				for(CMLPeak peak : peakList){
					if(peak.getMoleculeRefs() != null){
						if(peak.getMoleculeRefs()[0].equals(idMol)){
							String mass = Double.toString(peak.getXValue());
							String intensity = Double.toString(peak.getYValue());
							
							String numGroups = "0";
							List<CMLMetadata> metadataDescendants = spectrum.getMetadataListElements().get(0).getMetadataDescendants();
							for(CMLMetadata metadata : metadataDescendants){
								if(metadata.getDictRef().equals("nmc:groupPeakMSn"))
									numGroups = metadata.getContent();
							}
							molecule.setProperty(MZDataConstants.GROUP_PEAK_MSN, numGroups);
							molecule.setProperty(MZDataConstants.MASS, mass);
							molecule.setProperty(MZDataConstants.INTENSITIY, intensity);
						}
					}
				}
			}
		}
		

	}

	private String getRangeLine(MolecularFormulaRange range) {
		String rangeL = "";
		int countISO = 0;
		for (IIsotope isotope : range.isotopes()) {
            rangeL += isotope.getSymbol()+""+range.getIsotopeCountMin(isotope)+".."+range.getIsotopeCountMax(isotope);
            if(countISO < range.getIsotopeCount()-1)
            	rangeL += ",";
            countISO ++;
        }
		return rangeL;
	}

	/**
	 * Add the reference between signal in a spectrum and the molecule of the reaction
	 * 
	 * @param newCMLSpeList The CMLSpectrum object
	 * @param ionFinal      The ParentIon top
	 */
	private void addRefMol(List<CMLSpectrum> newCMLSpeList, ParentIon ionFinal) {

    	setRef(ionFinal,newCMLSpeList);
		for(int f2 = 0 ; f2 < ionFinal.getFragments().size(); f2++){
        	ParentIon ion2 = ionFinal.getFragments().get(f2);
        	setRef(ion2,newCMLSpeList);
    		for(int f3 = 0 ; f3 < ion2.getFragments().size(); f3++){
            	ParentIon ion3 = ion2.getFragments().get(f3);
            	setRef(ion3,newCMLSpeList);
            	for(int f4 = 0 ; f4 < ion3.getFragments().size(); f4++){
                	ParentIon ion4 = ion3.getFragments().get(f4);
                	setRef(ion4,newCMLSpeList);
                	for(int f5 = 0 ; f5 < ion4.getFragments().size(); f5++){
                    	ParentIon ion5 = ion4.getFragments().get(f5);
                    	setRef(ion5,newCMLSpeList);
                    	for(int f6 = 0 ; f6 < ion5.getFragments().size(); f6++){
                        	ParentIon ion6 = ion5.getFragments().get(f6);
                        	setRef(ion6,newCMLSpeList);
                        	for(int f7 = 0 ; f7 < ion6.getFragments().size(); f7++){
                            	ParentIon ion7 = ion6.getFragments().get(f7);
                            	setRef(ion7,newCMLSpeList);
                            }
                    	}
                	}
            	}
        	}
	    }
		
	}

	private void setRef(ParentIon ion, List<CMLSpectrum> newCMLSpeList) {
		if(ion.getFormulaSet() == null)
			return;
		for(CMLSpectrum spectrum:newCMLSpeList){
			for(CMLPeak peak:spectrum.getPeakListElements().get(0).getPeakDescendants()){
				if(peak.getId().equals(ion.getID_XCMS())){
					peak.setMoleculeRefs(ion.getID());
				}
			}
		}
	}

	/**
	 * Create a ParentIon given a peakID.
	 * 
	 * @param newCMLSpeList A list of CMLSpectrum
	 * @param peakIDP       The peakID as a string
	 * @return              The ParentIon object
	 */
	private ParentIon extractParentIon(List<CMLSpectrum> newCMLSpeList, String peakIDP, String parentScan) {
		CMLPeak cmlPeak = getCMLPeakFromIDP(newCMLSpeList,peakIDP);
		if(cmlPeak == null)
			return null;
		ParentIon parentIon1 = new ParentIon(correctingMass(cmlPeak.getXValue(),polarity),cmlPeak.getYValue(),""+countMol, 1);
		parentIon1.setScan(parentScan);
		countMol++;
		parentIon1.addID_XCMS(cmlPeak.getId());
		///////////////////////////////////////////////////////////////////////////
		Accuracy acc = new Accuracy();
		acc.setAccuracyPPM(parentIon1.getMass(), ppm.get(0));
		double accur = acc.getAccuracyDa();
		///////////////////////////////////////////////////////////////////////////
		parentIon1.setAccuracy(accur);
		HashMap<String,ParentIon> hashMap = new HashMap<String,ParentIon>();
		hashMap.put(""+0, parentIon1);
		List<CMLPeak> peakList = new ArrayList<CMLPeak>();
		for(CMLSpectrum spectrum:newCMLSpeList){
			for(CMLPeak peak:spectrum.getPeakListElements().get(0).getPeakDescendants()){
				peakList.add(peak);
			}
		}
		int count = 1;
		for(int i = 0 ; i < peakList.size(); i++){
			ParentIon parentIon = hashMap.get(""+i);
			if(parentIon == null)
				continue;
//			System.out.print(parentIon.getID_XCMS()+"("+parentIon.getMass()+")");
			for(CMLSpectrum spectrum2:newCMLSpeList){

				int nLevel = -1;
				String scanNum = "";
				scanNum = spectrum2.getId();
				List<CMLMetadata> ml = spectrum2.getMetadataListElements().get(0).getMetadataDescendants();
				for(CMLMetadata metadata:ml)
					if(metadata.getDictRef().equals(MZDataConstants.MS_LEVEL) )
						nLevel = Integer.parseInt(metadata.getContent());
//					else if(metadata.getDictRef().equals("nmc:scanNum") )
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
						if(((value + rangeParent) > parentIon.getMass()) && ((value - rangeParent) < parentIon.getMass()))
							precuB = true;
					}
					if(condition.getDictRef().equals(MZDataConstants.PRECURSOR_SCAN) ){
						String scan = condition.getValue();
						if(scan.equals(parentIon.getScan()))
							precuS = true;
					}
				}
				
				if(precuB && precuS)
					for(CMLPeak peak2:spectrum2.getPeakListElements().get(0).getPeakDescendants()){
						ParentIon newIon = new ParentIon(correctingMass(peak2.getXValue(),polarity),peak2.getYValue(), ""+countMol, nLevel);
						countMol++;
						newIon.addID_XCMS(peak2.getId());
						newIon.setScan(scanNum);
						if(newIon.getMass() > parentIon.getMass() - 0.1) // remove those peaks higher
							continue;
//						System.out.print(" -> "+spectrum2.getId()+":"+newIon.getID_XCMS()+"("+newIon.getMass()+")");

						Accuracy acc2 = new Accuracy();
						if(nLevel == 2)
							acc2.setAccuracyPPM(newIon.getMass(), ppm.get(1));
						else if(nLevel == 3)
							acc2.setAccuracyPPM(newIon.getMass(), ppm.get(2));
						else if(nLevel == 4)
							acc2.setAccuracyPPM(newIon.getMass(), ppm.get(3));
						else if(nLevel == 5)
							acc2.setAccuracyPPM(newIon.getMass(), ppm.get(4));
						else if(nLevel == 6)
							acc2.setAccuracyPPM(newIon.getMass(), ppm.get(5));
						else if(nLevel == 7)
							acc2.setAccuracyPPM(newIon.getMass(), ppm.get(6));
						
						double accur2 = acc2.getAccuracyDa();
						newIon.setAccuracy(accur2);
			        	parentIon.addFragment(newIon);
			        	newIon.addParent(parentIon);
			        	hashMap.put(""+count, newIon);
			        	count ++;
					}
				
			}
//			System.out.println();
		}
//		System.out.println("hashMap.size(): "+hashMap.size());
		return parentIon1;
	}

	/**
	 * Extract the CMLPeak given a peakID of the Parent Ion top
	 * 
	 * @param cmlSpectList A List of CMLSpectrum
	 * @param peakIDP      The peak ID
	 * @return             The CMLPeak object
	 */
	private CMLPeak getCMLPeakFromIDP(List<CMLSpectrum> cmlSpectList, String peakIDP) {
		Double peakValue = Double.parseDouble(peakIDP);
		if(massParent != null)
			peakValue = massParent;
//		System.out.println("     precursor using: "+peakValue+", rang:"+rangeParent+", from peakID:"+peakIDP);
		// TODO : windonw at the moment set to 1
		CMLPeak peakp = null;
		for(CMLSpectrum spectrum:cmlSpectList){
			boolean ms1 = false;
			List<CMLMetadata> ml = spectrum.getMetadataListElements().get(0).getMetadataDescendants();
			for(CMLMetadata metadata:ml)
				if(metadata.getDictRef().equals(MZDataConstants.MS_LEVEL) )
						if(metadata.getContent().equals("1")) // TOP
							ms1 = true;
			if(!ms1)
				continue;
			for(CMLPeak peak:spectrum.getPeakListElements().get(0).getPeakDescendants()){
				if(peak.getXValue() < (peakValue + rangeParent) && (peakValue -rangeParent) < peak.getXValue()){
					if(peakp != null){
						if(peak.getYValue() > peakp.getYValue())
							peakp =  peak;
					}else
						peakp =  peak;
				}
			}
		}
		return peakp;
	}

	/**
	 * Get a IsotopePattern object given a peak ID
	 * 
	 * @param cmlSpectList  A List with CMLSpectrum
	 * @param peakID        The peak ID
	 * @return              The IsotopePattern object
	 */
	private IsotopePattern getIsotopePatternFromIDP(List<CMLSpectrum> cmlSpectList, CMLPeak peakID) {
		List<CMLPeak> peakList = new ArrayList<CMLPeak>();
		for(CMLSpectrum spectrum:cmlSpectList){
			boolean ms1 = false;
			List<CMLMetadata> ml = spectrum.getMetadataListElements().get(0).getMetadataDescendants();
			for(CMLMetadata metadata:ml)
				if(metadata.getDictRef().equals(MZDataConstants.MS_LEVEL) )
						if(metadata.getContent().equals("1"))
							ms1 = true;
			if(!ms1)
				continue;
			for(CMLPeak peak:spectrum.getPeakListElements().get(0).getPeakDescendants()){
				if((peak.getXValue() > peakID.getXValue())  && (peakID.getXValue() + 3.5 > peak.getXValue()))
					if(peakID.getYValue() > peak.getYValue()){
						peakList.add(peak);
					}
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
	 * Generate default Elemental formula range containing C,H,N and O with a relation
	 * 50-1,100-1,30-0,30-0.
	 * 
	 * @return The MolecularFormulaRange object
	 */
	public static  MolecularFormulaRange creatingDefaultRange() {
		String[] stringElements = {"C","H","N","O"};
		int[] elementsMax = {50,100,30,30};
		int[] elementsMin = {1,1,0,0,0};
		return creatingRange(stringElements, elementsMax, elementsMin);
	}
	/**
	 * Get the MolecularFormulaRange object to generate.
	 * 
	 * @param stringElements
	 * @param elementsMax
	 * @param elementsMin
	 * @return                 The MolecularFormulaRange object
	 */
	public static  MolecularFormulaRange creatingRange(String[] stringElements,int[] elementsMax,int[] elementsMin) {
        
        MolecularFormulaRange mfRangei = new MolecularFormulaRange();
		try {
			IsotopeFactory ifac = IsotopeFactory.getInstance(DefaultChemObjectBuilder.getInstance());
	        for(int sE = 0 ; sE < stringElements.length ; sE++){
	        	if(elementsMax[sE] != 0)
	        		mfRangei.addIsotope(ifac.getMajorIsotope(stringElements[sE]),elementsMin[sE],elementsMax[sE]);
	        }
		} catch (IOException e) {
			e.printStackTrace();
		}
		return mfRangei;
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
	public static double correctingMass(double mass, Polarity polarity) {
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
	 * Correct the Mass inversly according the  Polarity. 
	 * Polarity.positive designed for positive mode: To compare
	 * with the Theoretical EF you need to add the mass of one electron.
	 * For negative mode the other way around. 
	 * 
	 * @param hashMap  HashMap with the masses 
	 * @param polarity The Polarity
	 * @return         The correct value mass
	 */
	static double correctingInverseMass(double mass, Polarity polarity) {
		/** Charge of the electron. */
		double chargeE = 0.00054857990927;
		
		if(polarity.equals(Polarity.positive))
			return (mass-chargeE);
		else if(polarity.equals(Polarity.negative))
			return (mass+chargeE);
		else
			return mass;

	}
	/**
	 * Conversion from ParentIon to reactionScheme
	 * 
	 * @param ionFinal
	 */
	private ReactionScheme convert2ReactionScheme(List<ParentIon> listPI) {

		String pol = "0.0";
		if(polarity == Polarity.positive)
			pol = "1.0";
		else if(polarity == Polarity.negative)
			pol = "-1.0";
		
		ReactionScheme reactionSet = new ReactionScheme();
		for(ParentIon ionFinal:listPI){
			ReactionScheme scheme1 = new ReactionScheme();
			reactionSet.add(scheme1);
			scheme1.setID("rs"+countScheme++);
			IMolecule reactant1 = new Molecule();
			reactant1.setID(ionFinal.getID());
			IMolecularFormulaSet ioFS = ionFinal.getFormulaSet();
			
			if(ioFS == null)
				reactant1.setProperty("cdk:partialCharge", pol);
			else if(ioFS.size() == 1){
				reactant1.setProperty(CDKConstants.FORMULA, ioFS.getMolecularFormula(0));
				reactant1.setProperty("cdk:partialCharge", pol);
			}else{
				reactant1.setProperty("cdk:partialCharge", pol);
				reactant1.setProperty(CDKConstants.FORMULA, ionFinal.getFormulaSet());
			}
			// exception if the parent ion doesn't contain any fragment
			if(ionFinal.getFragments().size() == 0){
				IReaction reaction0 = new Reaction();
		    	reaction0.setID("react"+countReaction);
				reaction0.addReactant(reactant1);
				scheme1.addReaction(reaction0);
			}
	        for(int f2 = 0 ; f2 < ionFinal.getFragments().size(); f2++){
	        	ParentIon ion2 = ionFinal.getFragments().get(f2);
	        	IMolecule reactant2 = new Molecule();
	    		reactant2.setID(ion2.getID());
	    		ReactionScheme scheme2 = creatingScheme(scheme1,reactant1,reactant2,ion2,countScheme++,countReaction++);
	    		
	    		if(scheme2 != null)
	        	for(int f3 = 0 ; f3 < ion2.getFragments().size(); f3++){
	            	ParentIon ion3 = ion2.getFragments().get(f3);
	            	IMolecule reactant3 = new Molecule();
	        		reactant3.setID(ion3.getID());
	        		ReactionScheme scheme3 = creatingScheme(scheme2,reactant2,reactant3,ion3,countScheme++,countReaction++);
	        		if(scheme3 != null)
	        		for(int f4 = 0 ; f4 < ion3.getFragments().size(); f4++){
	                	ParentIon ion4 = ion3.getFragments().get(f4);
	                	IMolecule reactant4 = new Molecule();
	            		reactant4.setID(ion4.getID());
	                	ReactionScheme scheme4 = creatingScheme(scheme3,reactant3,reactant4,ion4,countScheme++,countReaction++);
	            		
	                	if(scheme4 != null)
	    	        	for(int f5 = 0 ; f5 < ion4.getFragments().size(); f5++){
	                    	ParentIon ion5 = ion4.getFragments().get(f5);
	                    	IMolecule reactant5 = new Molecule();
	                		reactant5.setID(ion5.getID());
	                    	ReactionScheme scheme5 = creatingScheme(scheme4,reactant4,reactant5,ion5,countScheme++,countReaction++);
	                		
	                    	if(scheme5 != null)
	        	        	for(int f6 = 0 ; f6 < ion5.getFragments().size(); f6++){
	                        	ParentIon ion6 = ion5.getFragments().get(f6);
	                        	IMolecule reactant6 = new Molecule();
	                    		reactant6.setID(ion6.getID());
	                        	ReactionScheme scheme6 = creatingScheme(scheme5,reactant5,reactant6,ion6,countScheme++,countReaction++);
	                    		
	                        	if(scheme6 != null)
	            	        	for(int f7 = 0 ; f7 < ion6.getFragments().size(); f7++){
	                            	ParentIon ion7 = ion6.getFragments().get(f7);
	                            	IMolecule reactant7 = new Molecule();
	                        		reactant7.setID(ion7.getID());
	                            	ReactionScheme scheme7 = creatingScheme(scheme6,reactant6,reactant7,ion7,countScheme++,countReaction++);
	                        		
	                            }
	                    	}
	                	}
	            	}
	        	}
		    }
		}
		return reactionSet;
	}

	private ReactionScheme creatingScheme(ReactionScheme scheme1,
			IMolecule reactant1, IMolecule product2, ParentIon ion2,
			int countScheme, int countReaction) {
		
		IReaction reaction2 = new Reaction();
    	reaction2.setID("react"+countReaction);
		reaction2.addReactant(reactant1);
		
		String pol = "0.0";
		if(polarity == Polarity.positive)
			pol = "1.0";
		else if(polarity == Polarity.negative)
			pol = "-1.0";
		
		if(ion2.getFormulaSet() == null)
			return null;
		
		if(ion2.getFormulaSet().size() == 1){
			product2.setProperty(CDKConstants.FORMULA, ion2.getFormulaSet().getMolecularFormula(0));
			product2.setProperty("cdk:partialCharge", pol);
		}else{
			product2.setProperty("cdk:partialCharge", pol);
			product2.setProperty(CDKConstants.FORMULA, ion2.getFormulaSet());
		}

		product2.setProperty(MZDataConstants.PRECURSOR_ID, reactant1.getID());
		
		IMoleculeSet productSet = builder.newMoleculeSet();
		productSet.addAtomContainer(product2);
		
		IMolecule productLoss = reaction2.getBuilder().newMolecule();
		productLoss.setID(ion2.getID()+"Loss");
		productLoss.setProperty(MZDataConstants.PRECURSOR_ID, reactant1.getID());
		if(ion2.getFormulaLossSet().size() == 1){
			productLoss.setProperty(CDKConstants.FORMULA, ion2.getFormulaLossSet().getMolecularFormula(0));
			productLoss.setProperty("cdk:partialCharge", "0.0");
			productSet.addAtomContainer(productLoss);
		}else{
			productLoss.setProperty("cdk:partialCharge", "0.0");
			productLoss.setProperty(CDKConstants.FORMULA, ion2.getFormulaLossSet());
			productSet.addAtomContainer(productLoss);
		}
		reaction2.setProducts(productSet);
		
		if(ion2.getFragments().size() != 0){
			ReactionScheme scheme = new ReactionScheme();
			scheme.setID("rs"+countScheme++);
    		scheme1.add(scheme);
    		scheme.addReaction(reaction2);
    		return scheme;
    	}else{
    		scheme1.addReaction(reaction2);
    		return scheme1;
    	}
		
	}
}
