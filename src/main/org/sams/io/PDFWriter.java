package org.sams.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.io.StringWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.rosuda.REngine.Rserve.RConnection;
import org.rosuda.REngine.Rserve.RserveException;
import org.sams.MZData;
import org.sams.MZDataConstants;
import org.sams.manipulator.MZDataManipulator;

/**
 * Class to write PDF files
 * 
 * @author Miguel Rojas-Cherto
 */
public class PDFWriter {

    private String output;
	private String digits = "#";
	private String title = "";
	private String subtitle = "";
	private String type = MZDataConstants.FRAGMENTS;

	/**
     * Constructs a new PDFWriter class. Output will be stored in the Writer
     * class given as parameter. 
     *
     * @param out Path to redirect the output to.
     */
    public PDFWriter(String out) {
    	this.output = out;
    }

    /**
     * Serializes the MZData to PDF and redirects it to the output Writer.
     *
     * @param object A MZData object
     */
	public void write(MZData mzData) {
		write(mzData, "");
	}
    /**
     * Serializes the MZData to PDF and redirects it to the output Writer.
     *
     * @param object A MZData object
     * @param String The InChI identifier
     */
	public void write(MZData mzData, String inchi) {
		String matrix = MZDataManipulator.toString(mzData,false);
//		System.out.println(matrix);
		//"peakID","MSnParentPeakID","msLevel","sy","mass","los","color","colorLo"
		StringWriter outputW = new StringWriter();
		outputW.write(matrix);
		boolean colors = false;
		try {
			
			BufferedReader reader = new BufferedReader(new StringReader(outputW.toString()));
			String str;
			
			DecimalFormat currency = new DecimalFormat(digits );
    	
		
			RConnection c = new RConnection();
			c.eval("rm(list = ls(all = TRUE))");
			c.eval("library('Rgraphviz')\n");
			c.eval("inchi = '"+inchi+"'\n");
			
			List<List<String>> peakIDList = new ArrayList<List<String>>();
			List<List<String>> MSnParentPeakIDList = new ArrayList<List<String>>();
			List<List<String>> msLevelList = new ArrayList<List<String>>();
			List<List<String>> syList = new ArrayList<List<String>>();
			List<List<String>> massSList = new ArrayList<List<String>>(); 
			List<List<String>> massSLoList = new ArrayList<List<String>>(); 
			List<List<String>> losList = new ArrayList<List<String>>();
			List<List<String>> colorList = new ArrayList<List<String>>();
			List<List<String>> colorLoList = new ArrayList<List<String>>();
			
			List<String> peakIDL = new ArrayList<String>();
			List<String> MSnParentPeakIDL = new ArrayList<String>();
			List<String> msLevelL = new ArrayList<String>();
			List<String> syL = new ArrayList<String>();
			List<String> massSL = new ArrayList<String>(); 
			List<String> massSLoL = new ArrayList<String>(); 
			List<String> losL = new ArrayList<String>();
			List<String> colorL = new ArrayList<String>();
			List<String> colorLoL = new ArrayList<String>();
			
			boolean conEC = false;
			String group = "0";
			while ((str = reader.readLine()) != null) {
				String[] ll = str.split(",");
				if(!group.equals(ll[0])){
					// store
					if(!group.equals("0")){
						peakIDList.add(peakIDL);
						MSnParentPeakIDList.add(MSnParentPeakIDL);
						msLevelList.add(msLevelL);
						syList.add(syL);
						massSList.add(massSL); 
						massSLoList.add(massSLoL); 
						losList.add(losL); 
						colorList.add(colorL);
						colorLoList.add(colorLoL);
					}
					
					// initiate again
					peakIDL = new ArrayList<String>();
					MSnParentPeakIDL = new ArrayList<String>();
					msLevelL = new ArrayList<String>();
					syL = new ArrayList<String>();
					massSL = new ArrayList<String>(); 
					massSLoL = new ArrayList<String>(); 
					losL = new ArrayList<String>();
					colorL = new ArrayList<String>();
					colorLoL = new ArrayList<String>();
					group = ll[0];
				}
				peakIDL.add(ll[1].replace("peak", ""));
				MSnParentPeakIDL.add(ll[2].replace("peak", ""));
				msLevelL.add(ll[8]);
				if(ll.length > 9){
					conEC = true;
					syL.add(ll[9]);
				}else
					syL.add("");
		        massSL.add(currency.format(Double.valueOf(ll[3])));
		        if(ll[14].equals("?"))
		        	massSLoL.add("");
		        else
		        	massSLoL.add(currency.format(Double.valueOf(ll[14])));
		        
//		        if(ll.length > 9)
					losL.add(ll[11]);
//		        else
//		        	losL.add("");
		        
//		        if(ll.length > 15){
		        if(!ll[15].equals("0"))
		        	colors = true;
		        colorL.add(ll[15]);
//		        }else
//		        	colorL.add("0");
		        
//		        if(ll.length > 16){
			        colorLoL.add(ll[16]);
//		        }else
//		        	colorLoL.add("0");
		        	
		       
			}
			// adding the last group
			{
				peakIDList.add(peakIDL);
				MSnParentPeakIDList.add(MSnParentPeakIDL);
				msLevelList.add(msLevelL);
				syList.add(syL);
				massSList.add(massSL); 
				massSLoList.add(massSLoL); 
				losList.add(losL);
				colorList.add(colorL);
				colorLoList.add(colorLoL);
			}
			reader.close();
			
			c.eval("pdf('"+output+"', height = 10, width = 10); \n");
			for(int j = 0 ; j < losList.size(); j++){
				c.eval("graphEx <- matrix(ncol=9,nrow="+losList.get(j).size()+")");
				for(int i = 0 ; i < losList.get(j).size(); i++){
//					System.out.println(i+" c("+peakIDList.get(j).get(i)+","+MSnParentPeakIDList.get(j).get(i)+","+msLevelList.get(j).get(i)+",'"+syList.get(j).get(i)+"',"+massSList.get(j).get(i)+",'"+losList.get(j).get(i)+"',"+massSLoList.get(j).get(i)+","+colorList.get(j).get(i)+","+colorLoList.get(j).get(i)+")");
					c.eval("graphEx["+(i+1)+",] <- c("+peakIDList.get(j).get(i)+","+MSnParentPeakIDList.get(j).get(i)+","+msLevelList.get(j).get(i)+",'"+syList.get(j).get(i)+"',"+massSList.get(j).get(i)+",'"+losList.get(j).get(i)+"',"+massSLoList.get(j).get(i)+","+colorList.get(j).get(i)+","+colorLoList.get(j).get(i)+")");
				}
				String script = "colnames(graphEx) <- c('peakID','MSnParentPeakID','msLevel','sy','mass','los','masslos','color','colorLo') \n"+
					"type = '"+type+"'\n";
				if(colors)
					script += extractRScript3();
				else if(conEC)
					script += extractRScript();
				else
					script += extractRScript2();
				c.eval(script);
			}
			c.eval("title(main = '"+title+"',sub = '"+subtitle+"')");
			c.eval("dev.off()");
			c.close();
			
		} catch (RserveException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	/**
	 * Set the digits of the mass to be show in the pdf plot.
	 * e.g. "#" for nominal mass or "#.##" two digits after coma.
	 * 
	 * @param digits A string setting the digits
	 */
	public void setDigitsMass(String digits){
		this.digits = digits;
	}
    /**
     * Flushes the output and closes this object
     */
	public void close() throws IOException {
	}
	
	/**
	 * Script to crate the pdf plot
	 * 
	 * @return A string with the script
	 */
	private String extractRScript(){
		String script = 
		"gm<-NULL \n"+
		"miniPlotTree <- function(graphEx,level,peak,gm){ \n"+
		"	if (graphEx[peak,'msLevel']>1) { \n"+
		"			gm <-rbind(gm,c(graphEx[peak,'MSnParentPeakID'], \n"+
		"							graphEx[peak,'peakID'])) \n"+
		"	} \n"+
		"	for (b in which((graphEx[,'msLevel']>level) \n"+
		"					& (graphEx[,'MSnParentPeakID'] == graphEx[peak,'peakID']) )) { \n"+
		"		gm <- miniPlotTree(graphEx,graphEx[b,'msLevel'],b,gm) \n"+
		"	} \n"+
		"	return(gm) \n"+
		"} \n"+
		"nodeNames=NULL \n"+
		"nodeLabels=NULL \n"+
		"nodeColors=NULL \n"+
		"for (a in 1:length(graphEx[,'peakID'])) { \n"+
		"    nodeNames <-c(nodeNames, paste(graphEx[a,'peakID'],sep='')) \n"+
		"    if(type=='nmc:losses'){\n" +
		"		nodeLabels <- c(nodeLabels, paste(a,'\n',graphEx[a,'los'],'\n',graphEx[a,'masslos'],sep='')) \n"+
		"    }else{\n" +
		"		nodeLabels <- c(nodeLabels, paste(a,'\n',graphEx[a,'sy'],'\n',graphEx[a,'mass'],sep='')) \n"+
		"    }\n"+
		"	 if(a == 1 && type=='nmc:losses'){\n" +
		"		nodeColors <- c(nodeColors, 'black') \n" +
		"	 }else{\n" +
		"       nodeColors <- c(nodeColors, 'white')\n" +
		"	 }\n"+
		"} \n"+
		"#print(as.data.frame(graphEx));\n" +
		"nodes = list(label = nodeNames) \n"+
		"nodes$label[nodeNames]  = nodeLabels \n"+
		"nodes$fixedsize[nodeNames] <- TRUE \n"+
		"nodes$fill[nodeNames] = nodeColors \n"+
		"nodes$height[nodeNames] <- nodes$width[nodeNames] <- 1 \n "+
		"gm <- miniPlotTree(graphEx,1,1,gm) \n"+
		"colnames(gm) = c('fromID','toID') \n"+
		"V <- as.character(unique(as.vector(gm))) \n"+
		"OG <- new('graphNEL',nodes=V) \n"+
		"edgemode(OG) <- 'directed' \n"+
		"OG <- addEdge(from=as.character(gm[,'fromID']), \n"+
		"              to=as.character(gm[,'toID']), \n"+
		"              graph=OG,weights=1) \n"+
		"edgeNames=NULL \n"+
		"edgeLabels=NULL \n"+
		"for (a in 2:length(graphEx[,'peakID'])) { \n"+
		"    edgeNames <-c(edgeNames, edgeNames(OG)[a-1]) \n"+
		"    edgeLabels <- c(edgeLabels,graphEx[a,'los'] ) \n"+
		"} \n"+
		"edges = list(label = edgeNames) \n"+
		"#edges$label[edgeNames]  = edgeLabels \n"+
		"#nodeRenderInfo(OG) <- list(shape = 'ellipse') \n"+
		"#OG <- layoutGraph(OG, edgeAttrs=edges,nodeAttrs=nodes) \n"+
		"nodeRenderInfo(OG) <- list(fill=nodes$fill[nodeNames])\n"+
		"attrs <- list(node=list(fixedsize=FALSE))\n"+
		"OG <- layoutGraph(OG, nodeAttrs=nodes, attrs=attrs) \n"+
		"graph.par(list(graph = list(main = '', sub = inchi, cex.main = 1.8, cex.sub = 1.4))) \n" +
		"renderGraph(OG) \n";
		return script;
	}

	/**
	 * Script to crate the pdf plot
	 * 
	 * @return A string with the script
	 */
	private String extractRScript2(){
		String script = 
		"gm<-NULL \n"+
		"miniPlotTree <- function(graphEx,level,peak,gm){ \n"+
		"	if (graphEx[peak,'msLevel']>1) { \n"+
		"			gm <-rbind(gm,c(graphEx[peak,'MSnParentPeakID'], \n"+
		"							graphEx[peak,'peakID'])) \n"+
		"	} \n"+
		"	for (b in which((graphEx[,'msLevel']>level) \n"+
		"					& (graphEx[,'MSnParentPeakID'] == graphEx[peak,'peakID']) )) { \n"+
		"		gm <- miniPlotTree(graphEx,graphEx[b,'msLevel'],b,gm) \n"+
		"	} \n"+
		"	return(gm) \n"+
		"} \n"+
		"nodeNames=NULL \n"+
		"nodeLabels=NULL \n"+
		"nodeColors=NULL \n"+
		"for (a in 1:length(graphEx[,'peakID'])) { \n"+
		"    nodeNames <-c(nodeNames, paste(graphEx[a,'peakID'],sep='')) \n"+
		"    if(type=='nmc:losses'){" +
		"		nodeLabels <- c(nodeLabels, paste(a,'\n',graphEx[a,'los'],'\n',graphEx[a,'masslos'],sep='')) \n"+
		"    }else{\n" +
		"		nodeLabels <- c(nodeLabels, paste(a,'\n',graphEx[a,'sy'],'\n',graphEx[a,'mass'],sep='')) \n"+
		"    }\n"+
		"	 if(a == 1 && type=='nmc:losses')" +
		"		nodeColors <- c(nodeColors, 'black') \n" +
		"	 else" +
		"       nodeColors <- c(nodeColors, 'white')\n"+
		"} \n"+
		"nodes = list(label = nodeNames) \n"+
		"nodes$label[nodeNames]  = nodeLabels \n"+
		"nodes$fixedsize[nodeNames] <- TRUE \n"+
		"nodes$fill[nodeNames] = nodeColors \n"+
		"nodes$height[nodeNames] <- nodes$width[nodeNames] <- 1 \n "+
		"gm <- miniPlotTree(graphEx,1,1,gm) \n"+
		"colnames(gm) = c('fromID','toID') \n"+
		"V <- as.character(unique(as.vector(gm))) \n"+
		"OG <- new('graphNEL',nodes=V) \n"+
		"edgemode(OG) <- 'directed' \n"+
		"OG <- addEdge(from=as.character(gm[,'fromID']), \n"+
		"              to=as.character(gm[,'toID']), \n"+
		"              graph=OG,weights=1) \n"+
		"edgeNames=NULL \n"+
		"edgeLabels=NULL \n"+
		"for (a in 2:length(graphEx[,'peakID'])) { \n"+
		"    edgeNames <-c(edgeNames, edgeNames(OG)[a-1]) \n"+
		"    edgeLabels <- c(edgeLabels,graphEx[a,'los'] ) \n"+
		"} \n"+
		"edges = list(label = edgeNames) \n"+
//		"#edges$label[edgeNames]  = edgeLabels \n"+
//		"#nodeRenderInfo(OG) <- list(shape = 'ellipse') \n"+
//		"#OG <- layoutGraph(OG, edgeAttrs=edges,nodeAttrs=nodes) \n"+
		"nodeRenderInfo(OG) <- list(fill=nodes$fill[nodeNames])\n"+
		"attrs <- list(node=list(fixedsize=FALSE))\n"+
		"OG <- layoutGraph(OG, nodeAttrs=nodes, attrs=attrs) \n"+
		"graph.par(list(graph = list(main = '', sub = inchi, cex.main = 1.8, cex.sub = 1.4))) \n" +
		"renderGraph(OG) \n";
		return script;
	}

	/**
	 * Script to crate the pdf plot with colors
	 * 
	 * @return A string with the script
	 */
	private String extractRScript3(){
		String script = 
			"gm<-NULL \n"+
			"miniPlotTree <- function(graphEx,level,peak,gm){ \n"+
			"	if (graphEx[peak,'msLevel']>1) { \n"+
			"			gm <-rbind(gm,c(graphEx[peak,'MSnParentPeakID'], \n"+
			"							graphEx[peak,'peakID'])) \n"+
			"	} \n"+
			"	for (b in which((graphEx[,'msLevel']>level) \n"+
			"					& (graphEx[,'MSnParentPeakID'] == graphEx[peak,'peakID']) )) { \n"+
			"		gm <- miniPlotTree(graphEx,graphEx[b,'msLevel'],b,gm) \n"+
			"	} \n"+
			"	return(gm) \n"+
			"} \n"+
			"nodeNames=NULL \n"+
			"nodeLabels=NULL \n"+
			"nodeColors=NULL \n"+
			"for (a in 1:length(graphEx[,'peakID'])) { \n"+
			"    nodeNames <-c(nodeNames, paste(graphEx[a,'peakID'],sep='')) \n"+
			"    if(type=='nmc:losses'){" +
			"		nodeLabels <- c(nodeLabels, paste(a,'\n',graphEx[a,'los'],'\n',graphEx[a,'masslos'],sep='')) \n"+
			"    }else{\n" +
			"		nodeLabels <- c(nodeLabels, paste(a,'\n',graphEx[a,'sy'],'\n',graphEx[a,'mass'],sep='')) \n"+
			"    }\n"+
			"    if(type=='nmc:losses'){" +
			"		if(a == 1){" +
			"			nodeColors <- c(nodeColors, 'black') \n"+
			"		}else if(graphEx[a,'colorLo'] == 0){ \n" +
			"			nodeColors <-c(nodeColors, 'white') \n" +
			"	 	}else if(graphEx[a,'colorLo'] == 1){ \n" +
			"			nodeColors <-c(nodeColors, '#FFFF66') \n" +
			"	 	}else if(graphEx[a,'colorLo'] == 2){ \n" +
			"			nodeColors <-c(nodeColors, '#FFCC66') \n" +
			"	 	}else if(graphEx[a,'colorLo'] == 3){ \n" +
			"			nodeColors <-c(nodeColors, '#FF6600') \n" +
			"	 	}else if(graphEx[a,'colorLo'] == 4){ \n" +
			"			nodeColors <-c(nodeColors, '#FF0000') \n" +
			"	 	}else if(graphEx[a,'colorLo'] == 5){ \n" +
			"			nodeColors <-c(nodeColors, '#660000') \n" +
			"	 	}" +
			"    }else{" +
			"	 	if(graphEx[a,'color'] == 0){ \n" +
			"			nodeColors <-c(nodeColors, 'white') \n" +
			"	 	}else if(graphEx[a,'color'] == 1){ \n" +
			"			nodeColors <-c(nodeColors, '#FFFF66') \n" +
			"	 	}else if(graphEx[a,'color'] == 2){ \n" +
			"			nodeColors <-c(nodeColors, '#FFCC66') \n" +
			"	 	}else if(graphEx[a,'color'] == 3){ \n" +
			"			nodeColors <-c(nodeColors, '#FF6600') \n" +
			"	 	}else if(graphEx[a,'color'] == 4){ \n" +
			"			nodeColors <-c(nodeColors, '#FF0000') \n" +
			"	 	}else if(graphEx[a,'color'] == 5){ \n" +
			"			nodeColors <-c(nodeColors, '#660000') \n" +
			"	 	}" +
			"	 } \n"+
			"} \n"+
			"nodes = list(label = nodeNames) \n"+
			"nodes$label[nodeNames]  = nodeLabels \n"+
			"nodes$fixedsize[nodeNames] <- TRUE \n"+
			"nodes$fill[nodeNames] = nodeColors \n"+
			"nodes$height[nodeNames] <- nodes$width[nodeNames] <- 1 \n "+
			"gm <- miniPlotTree(graphEx,1,1,gm) \n"+
			"colnames(gm) = c('fromID','toID') \n"+
			"V <- as.character(unique(as.vector(gm))) \n"+
			"OG <- new('graphNEL',nodes=V) \n"+
			"edgemode(OG) <- 'directed' \n"+
			"OG <- addEdge(from=as.character(gm[,'fromID']), \n"+
			"              to=as.character(gm[,'toID']), \n"+
			"              graph=OG,weights=1) \n"+
			"edgeNames=NULL \n"+
			"edgeLabels=NULL \n"+
			"for (a in 2:length(graphEx[,'peakID'])) { \n"+
			"    edgeNames <-c(edgeNames, edgeNames(OG)[a-1]) \n"+
			"    edgeLabels <- c(edgeLabels,graphEx[a,'los'] ) \n"+
			"} \n"+
			"edges = list(label = edgeNames) \n"+
			"#edges$label[edgeNames]  = edgeLabels \n"+
			"#nodeRenderInfo(OG) <- list(shape = 'ellipse') \n"+
			"#OG <- layoutGraph(OG, edgeAttrs=edges,nodeAttrs=nodes) \n"+
			"nodeRenderInfo(OG) <- list(fill=nodes$fill[nodeNames])\n"+
			"attrs <- list(node=list(fixedsize=FALSE))\n"+
			"OG <- layoutGraph(OG, nodeAttrs=nodes, attrs=attrs) \n"+
			"graph.par(list(graph = list(main = '', sub = inchi, cex.main = 1.8, cex.sub = 1.4)))\n" +
			"renderGraph(OG) \n";
			return script;
	}

	/**
	 * Set a title in the pdf.
	 * 
	 * @param title    A string with the title to add
	 * @param subtitle A string with the subtitle to add
	 */
	public void setTitle(String title, String subtitle) {
		this.title = title;
		this.subtitle = subtitle;
	}

	/**
	 * Set the type of fragmentation tree to be displayed
	 * 
	 * @param fragments
	 */
	public void setType(String type) {
		this.type  = type;
	}
}
