package org.sams.io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Writer;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Map;
import java.util.UUID;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.openscience.cdk.ReactionScheme;
import org.rosuda.REngine.REXPMismatchException;
import org.rosuda.REngine.Rserve.RConnection;
import org.rosuda.REngine.Rserve.RserveException;
import org.sams.MZData;
import org.sams.MZDataConstants;
import org.sams.main.PackageVersion;
import org.sams.manipulator.MZDataManipulator;
import org.sams.spect.ScanMZXML;
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
 * Class to read mzXML files
 * 
 * @author Miguel Rojas-Cherto
 */
public class MZXMLReader {
	
    private BufferedReader input;
	private CMLSpectrumList listSpec;
	private double mzgap=0.2;
	private double snthresh=2;
	private String intrumentRef = "";
	private String protocol = "";
	private double rInt = 0;
	private int scansN;
	private int groupsN;
	private String version;

	/**
     * Reads MZXML from an java.io.InputStream, for example the FileInputStream.
     *
     * @param input InputStream type in
     */
    public MZXMLReader(InputStream in) {
    	PackageVersion v = new PackageVersion();
    	if(v.existManifest())
    		version = v.getAttribute(PackageVersion.ATTRIBUTE_IMPLEMENTATIONVERSION);
    	input = new BufferedReader(new InputStreamReader(in));
    }
    
    private void init() {
    	
		// TODO change method to read mzXML
    	//create tmp file with the input to be read from XCMS
    	StringBuilder contents = new StringBuilder();
        String line = null; //not declared within while loop
        try {
        	while (( line = input.readLine()) != null){
				contents.append(line);
				contents.append(System.getProperty("line.separator"));
			}
        	String pFile = "/tmp/"+UUID.randomUUID().toString()+".mzXML";
        	File aFile = new File(pFile);
	      	Writer output = new BufferedWriter(new FileWriter(aFile));
	        output.write( contents.toString() );
	        input.close();
	        output.close();
	        
	        listSpec = extractTable(pFile);
		} catch (IOException e) {
			e.printStackTrace();
		}
    }
    public void setMZGap(double mzgap){
    	this.mzgap = mzgap;
    }
    public void setSNThresh(double snthresh){
    	this.snthresh = snthresh;
    }
	private CMLSpectrumList extractTable(String pFile) {
		CMLSpectrumList specList = new CMLSpectrumList();
		CMLMetadataList metadataListL = new CMLMetadataList();
		specList.addMetadataList(metadataListL);
		
		try {
			
			DocumentBuilderFactory docBuilderFactory = DocumentBuilderFactory.newInstance();
	        DocumentBuilder docBuilder = docBuilderFactory.newDocumentBuilder();
	        Document doc = docBuilder.parse(pFile);
	        NodeList listOfScans = doc.getElementsByTagName("scan");

	        // normalize text representation
	        doc.getDocumentElement().normalize();
	        
			RConnection c = new RConnection();
			c.eval("rm(list = ls(all = TRUE))");
			c.eval("library(xcms)");
			c.eval("xset <- xcmsSet('"+pFile+"', method = 'MS1');");
			c.eval("xr <- xcmsRaw('"+pFile+"',includeMSn=TRUE)");
			c.eval("xcmsfrag <- xcmsFragments(xset,mzgap="+mzgap+",snthresh="+snthresh+")");
			
			if(rInt != 0.0){
				// removing those peaks according relative intensity
				c.eval("uniRT <- unique(xcmsfrag@peaks[,'rt'])");
				c.eval("relInt <- rep(1,length(xcmsfrag@peaks[,'rt']))");
				String loop = "for(rt in uniRT){ \n"+
					"poss <- which(xcmsfrag@peaks[,'rt'] == rt)\n"+
					"maxInt <- max(xcmsfrag@peaks[poss,'intensity'])\n"+
					"for(i in 1:length(poss)){\n"+
					"	relInt[poss[i]] <- xcmsfrag@peaks[poss[i],'intensity']/maxInt \n"+
	//				"	cat(poss[i],',',xcmsfrag@peaks[poss[i],'intensity']/maxInt,'=',xcmsfrag@peaks[poss[i],'intensity'],'/',maxInt,'\n') \n"+
					"}	\n"+
				"}";
				c.eval(loop);
				c.eval("rInt <- "+rInt);
	//			c.eval("cat(rInt,'\n')");
				c.eval("rmPI <- which(relInt <= rInt)");
//				c.eval("cat(rmPI)");
				

				// removing peaks out of the retention
				String iff = "if(length(rmPI) > 0) \n"+
					"xcmsfrag@peaks <- xcmsfrag@peaks[-rev(rmPI),]"; 
				c.eval(iff);
				
				// removing peaks not linked
				String rm = "poMS1 <- which(xcmsfrag@peaks[,'msLevel'] == 1);\n"+
				"poRM <- c();\n"+
				"for(i in 1:length(xcmsfrag@peaks[,1])){\n"+
				"	if( length(which(i == poMS1)) == 0){\n"+
				"		pos = which(xcmsfrag@peaks[,'peakID'] == xcmsfrag@peaks[i,'MSnParentPeakID']) \n"+
				"		if(length(pos) == 0){ \n"+
				"			poRM <- c(poRM,i) \n"+
				"		}else{\n"+
				"			posp <- which(xcmsfrag@peaks[,'peakID'] == xcmsfrag@peaks[i,'MSnParentPeakID'])[1]\n"+
				"			pos = which(poRM == posp)\n"+
				"			if(length(pos) != 0){\n"+
				"				poRM <- c(poRM,i)\n"+
				"			}\n"+
				"		}\n"+
				"	}\n"+
				"}\n";
				c.eval(rm);
				String poRM = "if(length(poRM) > 0) \n"+
				"xcmsfrag@peaks <- xcmsfrag@peaks[-rev(poRM),];";
				c.eval(poRM);
			}

			// grouping xcmsfrag
			String grouping =   "poMS1 <- which(xcmsfrag@peaks[,'msLevel'] == 1);\n"+
                                "for(i in poMS1){\n"+
                                "        xcmsfrag@peaks[i,'GroupPeakMSn'] <- i;\n"+
                                "}\n"+
                                "for(i in 1:length(xcmsfrag@peaks[,1])){\n"+
                                "        if( length(which(i == poMS1)) == 0){\n"+
                                "                pos = which(xcmsfrag@peaks[,'peakID'] == xcmsfrag@peaks[i,'MSnParentPeakID']);\n"+
                                "                xcmsfrag@peaks[i,'GroupPeakMSn'] <- xcmsfrag@peaks[pos,'GroupPeakMSn']\n"+
                                "        }\n"+
                                "}";
			c.eval(grouping);
			int numRow  = c.eval("length(xcmsfrag@peaks[,1])").asInteger();
			int numCol  = c.eval("length(xcmsfrag@peaks[1,])").asInteger();
			Object[][] data = new Object[numRow][numCol];
			String[] colnames = c.eval("colnames(xcmsfrag@peaks)").asStrings();
//			for(int i = 0; i < numCol; i++){
//				System.out.print(colnames[i]+" ");
//			}
//			System.out.println();
			for(int i = 1; i <= numRow; i++){
				String[] values = c.eval("xcmsfrag@peaks["+i+",]").asStrings();
				for(int j = 0; j< values.length; j++){
					data[i-1][j] = values[j];
//					System.out.print(values[j]+" ");
				}
//				System.out.println();
			}
			// restructure the links
			for(int i = 0; i < numRow; i++){
				String precID = (String)data[i][1];
				if(precID.equals("0.0"))
					continue;
				// extraction of the rt 
				String rt = "";
				double precMass = 0.0;
				for(int j = 0; j < numRow; j++){
					if(data[j][0].equals(precID)){
						precMass = Double.valueOf((String)data[j][4]);
						rt = (String)data[j][3];
						break;
					}
				}
				// look for those peaks with the same rt and a window of CID (0.5) for the highest peak
				double intensity = 0;
				for(int j = 0; j < numRow; j++){
					if(data[j][3].equals(rt)){
						double max = precMass+0.5;
						double min = precMass-0.5;
						double mass = Double.valueOf((String)data[j][4]);
						if(mass > max || mass < min)
							continue;
						if(intensity < Double.valueOf((String)data[j][5])){
							intensity = Double.valueOf((String)data[j][5]);
							precID = (String)data[j][0];
						}
					}
				}
				data[i][1] = precID;
			}
			DecimalFormat df = new DecimalFormat("#");
	        String rt_old = "";
			CMLPeakList cmlPeaks = null;
			int countScans = 0;
			int numGroups = 0;
			String polarity = "unknown";
			HashMap<String, String> mapScan = new HashMap<String, String>();
			for(int i = 0; i < data.length; i++){
				String rt = (String) data[i][3];
				if(!rt.equals(rt_old)){
					ScanMZXML scanMZXML = MZDataManipulator.getScanMZXML(Double.parseDouble(rt),listOfScans);
					c.eval("scan_pos <- which(xr@scantime == "+(String) data[i][3]+")");
					c.eval("scan_pos_msn <- which(xr@msnRt == "+(String) data[i][3]+")");
					rt_old = rt;
					CMLSpectrum cmlSpect = new CMLSpectrum();
					cmlSpect.setAttribute("type", "MS");
					CMLConditionList conditionList = new CMLConditionList();
					cmlSpect.addConditionList(conditionList);
					CMLMetadataList metadataList = new CMLMetadataList();
					cmlSpect.addMetadataList(metadataList);
//					int valueScan = countScans++;
					countScans++;
					cmlSpect.setId(Integer.toString(countScans));
					specList.addSpectrum(cmlSpect);
					/////////////////////////////////
					CMLMetadata metadata = new CMLMetadata();
					metadata.setDictRef(MZDataConstants.SCAN_NUM);
//					metadata.setContent(Integer.toString(valueScan));
					metadata.setContent(scanMZXML.getNum());
					metadataList.appendChild(metadata);
					/////////////////////////////////
					metadata = new CMLMetadata();
					metadata.setDictRef(MZDataConstants.MS_LEVEL);
					Double nLevel = Double.parseDouble((String) data[i][2]);  
					
					if(nLevel.equals(1.0))
						numGroups++;
					metadata.setContent(df.format(nLevel));
					metadataList.appendChild(metadata);
					/////////////////////////////////
					metadata = new CMLMetadata();
					metadata.setDictRef(MZDataConstants.RT);
					metadata.setContent((String) data[i][3]);
					metadataList.appendChild(metadata);
					mapScan.put((String) data[i][3], "" + (countScans));
					/////////////////////////////////
					metadata = new CMLMetadata();
					metadata.setDictRef(MZDataConstants.GROUP_PEAK_MSN);
					metadata.setContent(df.format(Double.parseDouble((String) data[i][7])));
					metadataList.appendChild(metadata);
					/////////////////////////////////
					CMLScalar condition = new CMLScalar();		
					condition.setDictRef(MZDataConstants.PRECURSOR_MZ);
					String precurMZ = "0";
					String precurRT = "0";
					String precurID = ((String) data[i][1]);
					for(int j = 0; j < data.length; j++){
						if(precurID.equals(((String) data[j][0]))){
							precurMZ = ((String) data[j][4]);
							precurRT = ((String) data[j][3]);
							break;
						}
					}
					condition.setValue(precurMZ);
					conditionList.appendChild(condition);
					/////////////////////////////////
					condition = new CMLScalar();		
					condition.setDictRef(MZDataConstants.PRECURSOR_SCAN);
					if(c.eval("(length(scan_pos_msn) > 0) != 0").asString().equals("true")){
						condition.setValue(mapScan.get(precurRT));
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
					String ce = "0";
					if(c.eval("(length(scan_pos_msn) > 0) != 0").asString().equals("true")){
						ce = c.eval("xr@msnCollisionEnergy[scan_pos_msn]").asString();
						condition.setValue(ce);
						conditionList.appendChild(condition);
					}
					/////////////////////////////////
					condition = new CMLScalar();		
					condition.setDictRef(MZDataConstants.ACTIVATION_METHOD);
					condition.setValue("CID");
					conditionList.appendChild(condition);
					/////////////////////////////////
					condition = new CMLScalar();		
					condition.setDictRef(MZDataConstants.POLARITY);
					if(c.eval("(length(scan_pos) > 0) != 0").asString().equals("true"))
						polarity = c.eval("xr@polarity[scan_pos]").asString();
					condition.setValue(polarity);
					conditionList.appendChild(condition);
					/////////////////////////////////
					cmlPeaks = new CMLPeakList();

					if(c.eval("(length(scan_pos) > 0) != 0").asString().equals("true")){
						c.eval("scan = getScan(xr,scan_pos)");
//						c.eval("peaks <- specPeaks(scan, sn = 0.7, mzgap = 0.01)");
						c.eval("peaks <- scan"); // TMP
						c.eval("ll <- order(peaks[,1])");
						c.eval("peaks <- peaks[ll,]");
						String id = df.format(Double.parseDouble( (String)data[i][0]));
						if(c.eval("length(ll)==1").asString().equals("true")){
							CMLPeak peak = new CMLPeak();
							String[] values = c.eval("peaks").asStrings();
							peak.setId(id);
							peak.setXValue((String) values[0]);
							
							peak.setYValue((String) values[1]);
							cmlPeaks.addPeak(peak);
						
						}else{
							int numRow2  = c.eval("length(peaks[,1])").asInteger();
							for(int k = 1; k <= numRow2; k++){
								CMLPeak peak = new CMLPeak();
								String[] values = c.eval("peaks["+k+",]").asStrings();
								peak.setId(id+"."+k);
								peak.setXValue((String) values[0]);
								
								peak.setYValue((String) values[1]);
								cmlPeaks.addPeak(peak);
							}
						}
						cmlSpect.addPeakList(cmlPeaks);
						
					}else{
						CMLPeak peak = new CMLPeak();
						peak.setId(df.format(Double.parseDouble( (String)data[i][0])));
						peak.setXValue((String) data[i][4]);
						peak.setYValue((String) data[i][5]);
						cmlPeaks.addPeak(peak);
						cmlSpect.addPeakList(cmlPeaks);
					}
				}else{
					CMLPeak peak = new CMLPeak();
					peak.setId(df.format(Double.parseDouble( (String)data[i][0])));
					peak.setXValue((String) data[i][4]);
					peak.setYValue((String) data[i][5]);
					cmlPeaks.addPeak(peak);
					
				}
			}

			/////////////////////////////////
			CMLMetadata metadataL = new CMLMetadata();
			metadataL.setDictRef(MZDataConstants.SCANS);
			metadataL.setContent(Integer.toString(countScans));
			metadataListL.appendChild(metadataL);
			scansN = countScans;
			/////////////////////////////////
			metadataL = new CMLMetadata();
			metadataL.setDictRef(MZDataConstants.NUM_GROUPS);
			metadataL.setContent(Integer.toString(numGroups));
			metadataListL.appendChild(metadataL);
			groupsN = numGroups;
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
			/////////////////////////////////parentFile
			NodeList listOfModels = doc.getElementsByTagName("parentFile");
	        Element modelElement = (Element)listOfModels.item(0);
	        if(modelElement != null){
	        	String fileName = modelElement.getAttribute("fileName");
	        	
				metadataL = new CMLMetadata();
//				metadataL.setDictRef(MZDataConstants.PARENT_FILE);
//				metadataL.setContent(fileName);
//		        metadataL.setConvention("RAWData");
//				metadataListL.appendChild(metadataL);
		        
				//////////////////////////////////////
				if(fileName.endsWith(".RAW"))
	        		fileName = fileName.replace(".RAW", ".mzXML");
	        	else if(fileName.endsWith(".raw"))
	        		fileName = fileName.replace(".raw", ".mzXML");
				
				metadataL.setDictRef(MZDataConstants.PARENT_FILE);
	        	metadataL.setContent(fileName);
	        	metadataL.setConvention("mzXMLData");
				metadataListL.appendChild(metadataL);
				
	        }
	        /////////////////////////////////MZGAP
	        metadataL = new CMLMetadata();
			metadataL.setDictRef(MZDataConstants.MZGAP);
			metadataL.setContent(Double.toString(mzgap));
			metadataListL.appendChild(metadataL);
	        /////////////////////////////////SNTHRESH
	        metadataL = new CMLMetadata();
			metadataL.setDictRef(MZDataConstants.SNTHRESH);
			metadataL.setContent(Double.toString(snthresh));
			metadataListL.appendChild(metadataL);
	        /////////////////////////////////RINT
	        metadataL = new CMLMetadata();
			metadataL.setDictRef(MZDataConstants.RINT);
			metadataL.setContent(Double.toString(rInt));
			metadataListL.appendChild(metadataL);
			
		} catch (RserveException e) {
			e.printStackTrace();
			return specList;
		} catch (REXPMismatchException e) {
			e.printStackTrace();
			return specList;
		} catch (ParserConfigurationException e) {
			e.printStackTrace();
			return specList;
		} catch (SAXException e) {
			e.printStackTrace();
			return specList;
		} catch (IOException e) {
			e.printStackTrace();
			return specList;
		}
		return specList;
	}
	/**
     * Read a MZData from input
     *
     * @return the content in a MZData object
     */
	public MZData read(MZData mzData) {
        init();
        mzData.setListSpectra(listSpec);
		mzData.setListReactions(new ReactionScheme());
		Map<Object,Object> properties = new HashMap<Object, Object>();
		properties.put(MZDataConstants.MZGAP, mzgap);
		properties.put(MZDataConstants.SNTHRESH, snthresh);
		properties.put(MZDataConstants.RINT, rInt);
		properties.put(MZDataConstants.NUM_GROUPS, groupsN);
		properties.put(MZDataConstants.SCANS, scansN);
		if(version != null)
			properties.put(MZDataConstants.VERSION_SAMS, version);
		mzData.setProperties(properties);
		return mzData;
	}

	/**
	 * 
	 * @param intrumentRef
	 */
	public void addInstrumentRef(String intrumentRef) {
		this.intrumentRef = intrumentRef;
	}
	public void addProtocol(String protocol) {
		this.protocol = protocol;
	}

	public void setRInt(double i) {
		this.rInt  = i;
	}
}
