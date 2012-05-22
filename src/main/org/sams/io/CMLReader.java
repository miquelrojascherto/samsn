package org.sams.io;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.ReactionScheme;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IMolecularFormulaSet;
import org.openscience.cdk.interfaces.IMoleculeSet;
import org.openscience.cdk.interfaces.IReactionScheme;
import org.openscience.cdk.interfaces.IReactionSet;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.openscience.cdk.tools.manipulator.ReactionSchemeManipulator;
import org.openscience.cdk.tools.manipulator.ReactionSetManipulator;
import org.sams.MZData;
import org.sams.MZDataConstants;
import org.sams.ext.ResetOnCloseInputStream;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;
import org.xmlcml.cml.element.CMLConditionList;
import org.xmlcml.cml.element.CMLMetadata;
import org.xmlcml.cml.element.CMLMetadataList;
import org.xmlcml.cml.element.CMLPeak;
import org.xmlcml.cml.element.CMLPeakList;
import org.xmlcml.cml.element.CMLScalar;
import org.xmlcml.cml.element.CMLSpectrum;
import org.xmlcml.cml.element.CMLSpectrumList;
/**
 * Class to read CML files
 * 
 * @author Miguel Rojas-Cherto
 */
public class CMLReader {
	
	private CMLSpectrumList cmlSpecList;
	private IReactionScheme reactionScheme;
	private HashMap<Object, Object> properties;

	/**
     * Constructor of the CMLReader object.
     * Reads CML from an java.io.InputStream, for example the FileInputStream.
     *
     * @param in InputStream type
     */
    public CMLReader(InputStream in) {
        init(in);
    }
    
    /**
     * Initiate the process extracting the information.
     * 
     * @param in InputStream type
     */
    private void init(InputStream in) {

        cmlSpecList = new CMLSpectrumList();
        reactionScheme = new ReactionScheme();
		properties = new HashMap<Object, Object>();
        
        try {
	    	DocumentBuilderFactory docBuilderFactory = DocumentBuilderFactory.newInstance();
	        DocumentBuilder docBuilder = docBuilderFactory.newDocumentBuilder();
	        ResetOnCloseInputStream decoratedIs = new ResetOnCloseInputStream(in);
	        Document doc = docBuilder.parse(decoratedIs);
	        // normalize text representation
            doc.getDocumentElement().normalize();

			////////////////////////////////////////
			CMLMetadataList metadataListList = new CMLMetadataList();
			NodeList mdl = doc.getElementsByTagName("metadataList");
        	for (int i = 0; i < mdl.getLength();i++) {
        		Element node = (Element) mdl.item(i);
        		NodeList ss = node.getElementsByTagName("metadata");
        		for(int j = 0 ; j < ss.getLength(); j++){
	            	Element scal = (Element) ss.item(j);
	            	CMLMetadata metadata = new CMLMetadata();
	            	metadata.setDictRef(scal.getAttribute("dictRef"));
	            	String content = scal.getAttribute("content");
					metadata.setContent(content);
					metadataListList.appendChild(metadata); 
					if(!scal.getAttribute("convention").equals(""))
						metadata.setConvention(scal.getAttribute("convention"));
					properties.put(scal.getAttribute("dictRef"), scal.getAttribute("content"));
				}
            	break;
        	}
			cmlSpecList.addMetadataList(metadataListList);
			////////////////////////////////////////
			HashMap<String, IChemObject> hashMol = new HashMap<String, IChemObject>();
			
    		NodeList listSpect = doc.getElementsByTagName("spectrum");
            for(int i = 0 ; i < listSpect.getLength(); i++){
            	Element spEl = (Element)listSpect.item(i);
            	
            	CMLSpectrum cmlSpect = new CMLSpectrum();
				cmlSpect.setAttribute("type", "MS");
				cmlSpect.setId(spEl.getAttribute("id"));
				////////////////////////////////////////
				CMLConditionList conditionList = new CMLConditionList();
				cmlSpect.addConditionList(conditionList);
				CMLMetadataList metadataList = new CMLMetadataList();
				cmlSpect.addMetadataList(metadataList);
				////////////////////////////////////////
				NodeList scalList = spEl.getElementsByTagName("scalar");
				for(int j = 0 ; j < scalList.getLength(); j++){
                	Element scal = (Element) scalList.item(j);
					CMLScalar condition = new CMLScalar();		
					condition.setDictRef(scal.getAttribute("dictRef"));
					condition.setValue(scal.getTextContent());
					conditionList.addScalar(condition);
				}
				////////////////////////////////////////
				String group = null;
				NodeList metaList = spEl.getElementsByTagName("metadata");
				for(int j = 0 ; j < metaList.getLength(); j++){
                	Element scal = (Element) metaList.item(j);
                	CMLMetadata metadata = new CMLMetadata();
					metadata.setDictRef(scal.getAttribute("dictRef"));
					metadata.setContent(scal.getAttribute("content"));
					metadataList.appendChild(metadata);
					if(metadata.getDictRef().equals("nmc:groupPeakMSn"))
						group = scal.getAttribute("content");
				}
				////////////////////////////////////////
				cmlSpecList.addSpectrum(cmlSpect);
				CMLPeakList cmlPeaks = new CMLPeakList();
				cmlSpect.addPeakList(cmlPeaks);
				NodeList peaklist = spEl.getElementsByTagName("peak");
            	for(int j = 0 ; j < peaklist.getLength(); j++){
                	Element peakEl = (Element)peaklist.item(j);
                	String id_p = peakEl.getAttribute("id");
                	Double xValue = Double.parseDouble(peakEl.getAttribute("xValue"));
                	Double yValue = Double.parseDouble(peakEl.getAttribute("yValue"));
                	String moleculeRefs = peakEl.getAttribute("moleculeRefs");
                	CMLPeak peak = new CMLPeak();
					peak.setId(id_p.toString());
					peak.setXValue(xValue.toString());
					peak.setYValue(yValue.toString());
					if(moleculeRefs != null){
						peak.setMoleculeRefs(moleculeRefs);
						
						IChemObject object = reactionScheme.getBuilder().newChemObject();
						object.setProperty(MZDataConstants.GROUP_PEAK_MSN, group);
						object.setProperty(MZDataConstants.MASS, xValue.toString());
						object.setProperty(MZDataConstants.INTENSITIY, yValue.toString());
						hashMol.put(moleculeRefs, object);
					}
					cmlPeaks.addPeak(peak); 
            	}
            }
            decoratedIs.close();
	    	// read list reactions
            org.openscience.cdk.io.CMLReader reader = new org.openscience.cdk.io.CMLReader(decoratedIs);
    	    IChemFile chemFile = (IChemFile)reader.read(new org.openscience.cdk.ChemFile());
    	    IReactionSet reactionSet = chemFile.getChemSequence(0).getChemModel(0).getReactionSet();
    	    if(reactionSet != null){
    	    	reactionScheme = ReactionSchemeManipulator.createReactionScheme(reactionSet);
    	    	IMoleculeSet allM = ReactionSchemeManipulator.getAllMolecules(reactionScheme);
    			IReactionSet reaSet = ReactionSchemeManipulator.getAllReactions(reactionScheme);
    			for(IAtomContainer mol : allM.molecules()){
    				IChemObject object = hashMol.get(mol.getID());
    				if(object != null){
	    				Map<Object, Object> chemProperties = object.getProperties();
	    				for( Object key: chemProperties.keySet() ){
	    				        Object property = key;
	    				        String value = (String) chemProperties.get( key );
	    				        mol.setProperty(property,value);
	    				}
	    				// setting the precursor molecule id
						IReactionSet reacSetR = ReactionSetManipulator.getRelevantReactionsAsProduct(reaSet, mol);
						if(reacSetR.getReactionCount() > 0){
							String precursorID = reacSetR.getReaction(0).getReactants().getAtomContainer(0).getID();
					        mol.setProperty(MZDataConstants.PRECURSOR_ID,precursorID);
						}
    				}
					
    				//The property formula needs to be convert to IMolecularFormula
        	    	ArrayList<String> formSt = (ArrayList<String>) mol.getProperty(CDKConstants.FORMULA);
					IMolecularFormulaSet mfSet = mol.getBuilder().newMolecularFormulaSet();
    				for(String formula :formSt){
    					mfSet.addMolecularFormula(MolecularFormulaManipulator.getMolecularFormula(formula,chemFile.getBuilder()));
    				}
    				mol.setProperty(CDKConstants.FORMULA, mfSet);
    			}
    	    }
    	    
            decoratedIs.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		} catch (ParserConfigurationException e) {
			e.printStackTrace();
		} catch (SAXException e) {
			e.printStackTrace();
		} catch (CDKException e) {
			e.printStackTrace();
		}
    	
    }

	/**
     * Read a MZData from input.
     *
     * @return the content in a MZData object
     */
	public MZData read(MZData mzData) {
		mzData.setListSpectra(cmlSpecList);
		mzData.setListReactions(reactionScheme);
		mzData.setProperties(properties);
		return mzData;
	}
}
