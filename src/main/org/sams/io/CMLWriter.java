package org.sams.io;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;
import java.io.UnsupportedEncodingException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import nu.xom.Document;
import nu.xom.Element;
import nu.xom.Serializer;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMoleculeSet;
import org.openscience.cdk.io.cml.CustomSerializer;
import org.openscience.cdk.libio.cml.Convertor;
import org.openscience.cdk.tools.manipulator.ReactionSchemeManipulator;
import org.sams.MZData;
import org.sams.MZDataConstants;
import org.xmlcml.cml.base.CMLElements;
import org.xmlcml.cml.element.CMLMolecule;
import org.xmlcml.cml.element.CMLMoleculeList;

/**
 * Class to write CML files
 * 
 * @author Miguel Rojas-Cherto
 */
public class CMLWriter {

    private OutputStream output;
	private Writer writer;

	public CMLWriter(OutputStream output) {
        this.output = output;
    }
	/**
     * Constructs a new CMLWriter class. Output will be stored in the Writer
     * class given as parameter. The CML code will be valid CML code with a
     * XML header. Only one object can be stored.
     *
     * @param out Writer to redirect the output to.
     */
    public CMLWriter(Writer out) {
        this.writer = out;
        this.output = new OutputStream() {
			public void write(int anInt) throws IOException {
				writer.write(anInt);
			}
        };
    }

    /**
     * Serializes the MZData to CML and redirects it to the output Writer.
     *
     * @param object A MZData object
     */
	public void write(MZData mzData) {
        Element root = new Element("cml");
        root.addNamespaceDeclaration("xsi", "http://www.w3.org/2001/XMLSchema-instance");
        if(mzData.getListSpectra() != null){
        	root.appendChild(mzData.getListSpectra().copy());
        }
        if(mzData.getListReactions() != null){
        	Convertor convertor = new Convertor(false,null);
//        	try {
	        	//removing those properties not chemObject
//        		IReactionScheme reacClone = (IReactionScheme) mzData.getListReactions().clone();
	        	IMoleculeSet allM = ReactionSchemeManipulator.getAllMolecules(mzData.getListReactions());
	        	for(IAtomContainer mol : allM.molecules()){
    				mol.removeProperty(MZDataConstants.GROUP_PEAK_MSN);
    				mol.removeProperty(MZDataConstants.MASS);
    				mol.removeProperty(MZDataConstants.INTENSITIY);
    				mol.removeProperty(MZDataConstants.PRECURSOR_ID);
    			}
				
	        	//order molecules by alphabetic
	        	CMLMoleculeList cmlMolList = convertor.cdkMoleculeSetToCMLList(allM);
	        	CMLElements<CMLMolecule> cmlElList = cmlMolList.getMoleculeElements();
	        	List<String> mapValues = new ArrayList<String>();
	        	for(int i = 0 ; i < cmlElList.size(); i++){
	        		mapValues.add(cmlElList.get(i).getId());
	        	}
	        	Collections.sort(mapValues);	
	        	
	        	CMLMoleculeList cmlList = new CMLMoleculeList();
	            cmlList.setConvention(cmlMolList.getConvention());
	            cmlList.setId(cmlMolList.getId());
	            for(String id2:mapValues){
	        		int iF = 0;
	        		for(int i = 0 ; i < cmlElList.size(); i++){
		        		if(cmlElList.get(i).getId().equals(id2)){
		        			iF = i;
		        			break;
		        		}
		        			
		        	}
	        		cmlList.appendChild(cmlElList.get(iF));
	        	}
	        	
	        	
	        	root.appendChild(cmlList.copy());
            	convertor.setRef(true);
            	root.appendChild(convertor.cdkReactionSchemeToCMLReactionScheme(mzData.getListReactions()).copy());
//			} catch (CloneNotSupportedException e) {
//				e.printStackTrace();
//			}
        }
        Document doc = new Document(root);
        try {
		    Serializer serializer = new CustomSerializer(output, "ISO-8859-1");
		    serializer.setIndent(2);
	        serializer.write(doc);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (UnsupportedEncodingException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}

    /**
     * Flushes the output and closes this object
     */
	public void close() throws IOException {
        output.close();
	}
}
