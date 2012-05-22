package org.sams;

import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.openscience.cdk.interfaces.IReactionScheme;
import org.xmlcml.cml.element.CMLMetadata;
import org.xmlcml.cml.element.CMLSpectrumList;

/**
 * Class containing information about chemical structures and mass 
 * spectrometry.
 * 
 * @author Miguel Rojas-Cherto
 */
public class MZData {

	private String id;
	private CMLSpectrumList listSpect;
	private IReactionScheme listReact;
	/**
	 *  A hashtable for the storage of any kind of properties of this IChemObject.
	 */
	private Map<Object, Object> properties;
	

	/**
	 * Constructor of MZData object
	 */
	public MZData(){
		properties = new HashMap<Object, Object>();
	}

	/**
	 * Set a list of Spectra
	 * 
	 * @param listSpec The CMLSpectrumList object
	 */
	public void setListSpectra(CMLSpectrumList listSpec) {
		this.listSpect = listSpec;
	}

	/**
	 * 
	 * Get the list of Spectra
	 * 
	 * @return The CMLSpectrumList object
	 */
	public CMLSpectrumList getListSpectra(){
		return listSpect;
	}

	/**
	 * Set a list of reactions
	 * 
	 * @param listReact The IReactionscheme object
	 */
	public void setListReactions(IReactionScheme listReact) {
		this.listReact = listReact;
	}

	/**
	 * Get a list of reactions
	 * 
	 * @return The IReactionscheme object
	 */
	public IReactionScheme getListReactions(){
		return listReact;
	}

	/**
	 * Get the Id of the mzData
	 * 
	 * @return The id
	 */
	public String getID() {
		return id;
	}

	/**
	 * Set the Id of the mzData
	 * 
	 * @param idmzData The id
	 */
	public void setID(String idmzData) {
		id = idmzData;
	}
	
	/**
	 * Sets a property for a MZData.
	 *
	 * @param  description  An object description of the property (most likely a
	 *                      unique string)
	 * @param  property     An object with the property itself
	 * @see                 #getProperty
	 */
	public void setProperty(Object description, Object property){
		properties.put(description, property);
		setPropertyIn(description, property);
	}

	/**
	 * Set property in spectrum metadata
	 * 
	 * @param description Object with the description
	 * @param property    Object with the property
	 */
	private void setPropertyIn(Object description, Object property) {
		if(listSpect != null){
			List<CMLMetadata> metadataList = listSpect.getMetadataListElements().get(0).getMetadataDescendants();
	        boolean co = false;
			for(CMLMetadata metadata:metadataList){
				if(metadata.getDictRef().equals(description)){
					if(property instanceof Integer)
						metadata.setContent(Integer.toString((Integer)property));
					else if(property instanceof Double)
						metadata.setContent(Double.toString((Double)property));
					else
						metadata.setContent((String)property);
					co = true;
				}
			}
			if(!co){
				CMLMetadata metadataL = new CMLMetadata();
				metadataL.setDictRef((String)description);
				if(property instanceof Integer)
					metadataL.setContent(Integer.toString((Integer)property));
				else if(property instanceof Double)
					metadataL.setContent(Double.toString((Double)property));
				else
					metadataL.setContent((String)property);
				listSpect.getMetadataListElements().get(0).addMetadata(metadataL);
			}
		}
	}

	/**
	 * Returns a property for the MZData.
	 *
	 * @param  description  An object description of the property (most likely a
	 *                      unique string)
	 * @return              The object containing the property. Returns null if
	 *                      property is not set.
	 * @see                 #setProperty
	 */
	public Object getProperty(Object description){
		return properties.get(description);
	}
	
	/**
	 *  Returns a Map with the MZData's properties.
	 *
	 *@return    The object's properties as an Map
	 *@see       #setProperties
	 */
	public Map<Object,Object> getProperties(){
		return properties;
	}
	
	/**
	 * Sets the properties of this MZData.
	 *
	 * @param  properties  a Map specifying the property values
	 * @see                #getProperties
	 */
	public void setProperties(Map<Object,Object> properties){
		this.properties = properties;
		Set set = properties.entrySet();
		Iterator i = set.iterator();
		while(i.hasNext()){
	      Map.Entry me = (Map.Entry)i.next();
	      setPropertyIn(me.getKey(), me.getValue());
	    }
	}
	
	/**
	 * Clone the MZData object
	 * 
	 * @param mzData Object to clone
	 * @return       The Object cloned
	 */
	public MZData clone()throws CloneNotSupportedException {
		MZData newMZData = new MZData();
		// cloning reactionscheme
		IReactionScheme reactionS = (IReactionScheme) listReact.clone();
		newMZData.setListReactions(reactionS);
		// cloning spectra
		CMLSpectrumList spectList = (CMLSpectrumList) listSpect.copy();
		newMZData.setListSpectra(spectList);
		return newMZData;
	}

}
