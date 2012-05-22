package org.sams.main;

import java.io.BufferedInputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;

import ml.options.OptionData;
import ml.options.OptionSet;
import ml.options.Options;
import ml.options.Options.Multiplicity;
import ml.options.Options.Separator;

import org.openscience.cdk.formula.MolecularFormulaRange;
import org.sams.MEFgenerator;
import org.sams.MZData;
import org.sams.MZDataConstants;
import org.sams.io.CMLReader;
import org.sams.io.CMLWriter;
import org.sams.io.MZXMLReader;
import org.sams.io.PDFWriter;
import org.sams.manipulator.MZDataManipulator;

/**
 * Main class to process spectral trees
 * 
 * @author Miguel Rojas-Cherto
 *
 */
public class ProcessMZData {
	/**
	 * E.g. -sn=1 -mzgap=0.5 -rint=0 -acc=5 -e=[C1..15,H1..9,O0..4,N0..2] -rules=[RDBER,nitrogenR] -imzXML  src/test/data/mzXML/F002169_C15H9O4.mzXML -ocml  test.cml process
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		PackageVersion version = new PackageVersion();
		if(args.length == 1){
			if(args[0].equals("-H")){
				System.out.println("MEF processes and enriches spectral tree data");
				System.out.println("");
				System.out.println("Usage: java -jar file.jar [Parameters] <input spec> <ouput spec> [action]");
				System.out.println("");
				System.out.println("Each spec can be a file whose extension decides the format. Optionally,");
				System.out.println("the format can be specified by preceding the file by -i<format-type> ");
				System.out.println("e.g. -imzXML, for input and -o<format-type> for output.");
				System.out.println("");
				System.out.println("See below for available format-types, which are the same as the ");
				System.out.println("file extensions and are case independent.");
				System.out.println("If no output file is defined stdout is used instead.");
				System.out.println("");
				System.out.println("Parameters");
				System.out.println("-acc=    Maximal accuracy range. Set in ppm E.g. -acc=5");
				System.out.println("-e=      Elements to be include, together with the");
				System.out.println("         the upper-/lower-limit of the number of atoms");
				System.out.println("         e.g. -e=[C1..15,H1..9,O0..4,N0..2]");
				System.out.println("-mzGap=  Minimal distance between adjacent peaks. E.g. -mzgap=0.5");
				System.out.println("-rules=  Constrain rules applied to the formula");
				System.out.println("         e.g. -rules=[RDBER,nitrogenR]");
				System.out.println("-rt=     Relative intensity threshold. E.g. -rint=0");
				System.out.println("-sn=     Signal to noise threshold. E.g. -sn=1");
				System.out.println("-F       List of recognized file formats use");
				System.out.println("-R       List of recognized constraint formula rules");
				System.out.println("-v       Display version information");
				System.out.println("-H       Outputs this help text");
				System.out.println("If no parameters are defined default settings are used.");
				System.out.println("");
				System.out.println("Actions");
				System.out.println("convert -- Convert into other format mass spectral tree data");
				System.out.println("process -- Processes and enriches mass spectral tree data");
				System.out.println("compare -- Compare two mass spectral tree data");
				System.out.println("");
				
				System.exit(1);
			}else if(args[0].equals("-v")){
	            version.printVersion(); 
				System.exit(1);
				
			}else if(args[0].equals("-R")){
				System.out.println("RDBER -- Rings plus Double Bonds Equivalent rule");
				System.out.println("nitrogenR -- Nitrogen rule");
				
				System.exit(1);
			}else if(args[0].equals("-F")){
				System.out.println("txt -- Plain Text");
				System.out.println("cml -- Chemical Markup Language");
				System.out.println("mzXML -- mzXML format [Read-only]");
				System.out.println("pdf -- pdf format [Write-only]");
				
				System.exit(1);
			}
		}
		
		Options opt1 = new Options(args,1);
		opt1.addSet("cset")
			.addOption("sn", false, Separator.EQUALS, Multiplicity.ZERO_OR_ONE)
			.addOption("mzgap", false, Separator.EQUALS, Multiplicity.ZERO_OR_ONE)
			.addOption("rint", false, Separator.EQUALS, Multiplicity.ZERO_OR_ONE)
			.addOption("acc", false, Separator.EQUALS, Multiplicity.ZERO_OR_ONE)
			.addOption("e", false, Separator.EQUALS, Multiplicity.ZERO_OR_ONE)
			.addOption("rules", false, Separator.EQUALS, Multiplicity.ZERO_OR_ONE)
			.addOption("occurr", false, Separator.EQUALS, Multiplicity.ZERO_OR_ONE)
			.addOption("i", true, Separator.BLANK,Multiplicity.ONCE)
			.addOption("o", true, Separator.BLANK,Multiplicity.ZERO_OR_ONE);
		opt1.addSet("pset")
		.addOption("i1", true, Separator.BLANK,Multiplicity.ONCE)
		.addOption("i2", true, Separator.BLANK,Multiplicity.ONCE);

		// [Action] is the data
		
		OptionSet set = opt1.getMatchingSet();
		if (set == null) {
			System.out.println("No output file or format correctly spec!");
			System.out.print("SAMS ");
			version.printVersion(); 
			System.out.println("Usage: java -jar file.jar [ Parameters ... ] [-i<input-type>] <inputfile> [-o<ouput-type>] <outputfile> [Action]");
			System.out.println("Try  -H option for more information."); 
			System.exit(1);
		}else if(set.getSetName().equals("pset")) {
			String format_input1 = null;
			String file_input1 = null;
			String format_input2 = null;
			String file_input2 = null;
			String currentDir = new File("").getAbsolutePath();
			OptionData d1 = set.getOption("i1");
			for (int i = 0; i < d1.getResultCount(); i++) {
				format_input1 = d1.getResultDetail(i);
				file_input1 = d1.getResultValue(i);
				if(!(new File(file_input1)).exists()){
					file_input1 = currentDir+"/"+file_input1;
				}
				if(!(new File(file_input1)).exists()){
					System.out.println("The current file doesn't exist: "+file_input2); 
					System.exit(1);
				}
			}
			System.out.print(".");
			OptionData d2 = set.getOption("i2");
			for (int i = 0; i < d2.getResultCount(); i++) {
				format_input2 = d2.getResultDetail(i);
				file_input2 = d2.getResultValue(i);
				if(!(new File(file_input2)).exists()){
					file_input2 = currentDir+"/"+file_input2;
				}
				if(!(new File(file_input2)).exists()){
					System.out.println("The current file doesn't exist: "+file_input2); 
					System.exit(1);
				}
					
			}
			System.out.print(".");
			InputStream input1 = new BufferedInputStream(new FileInputStream(file_input1));
			CMLReader mzReader1 = new CMLReader(input1);
			MZData mzData1 = new MZData();
			mzData1 = mzReader1.read(mzData1);
			
			System.out.print(".");
			InputStream input2 = new BufferedInputStream(new FileInputStream(file_input2));
			CMLReader mzReader2 = new CMLReader(input2);
			MZData mzData2 = new MZData();
			mzData2 = mzReader2.read(mzData2);
			
			String[] nam1 = file_input1.split("/");
			String fileName1 = nam1[nam1.length-1].replace(".cml", "");
			String[] nam2 = file_input2.split("/");
			String fileName2 = nam2[nam2.length-1].replace(".cml", "");;
			String output12 = currentDir+"/"+fileName1+"_to_"+fileName2+".pdf";
			String output21 = currentDir+"/"+fileName2+"_to_"+fileName1+".pdf";
			String output12L = currentDir+"/"+fileName1+"_to_"+fileName2+"-los.pdf";
			String output21L = currentDir+"/"+fileName2+"_to_"+fileName1+"-los.pdf";
			
			System.out.print(".");
			MZData mzDataCol12 = MZDataManipulator.getColored(mzData1,mzData2);
			PDFWriter pdfWriter1 = new PDFWriter(output12);
			pdfWriter1.write(mzDataCol12);
			pdfWriter1.close();
			System.out.println("");
			System.out.println("Process finished, the "+output12+" file for fragments was created");
			
			MZData mzDataCol21 = MZDataManipulator.getColored(mzData2,mzData1);
			PDFWriter pdfWriter2 = new PDFWriter(output21);
			pdfWriter2.write(mzDataCol21);
			pdfWriter2.close();
			System.out.println("");
			System.out.println("Process finished, the "+output21+" file for fragments was created");
			
			PDFWriter pdfWriter1L = new PDFWriter(output12L);
			pdfWriter1L.setType(MZDataConstants.LOSSES);
			pdfWriter1L.write(mzDataCol12);
			pdfWriter1L.close();
			System.out.println("");
			System.out.println("Process finished, the "+output12L+" file for losses was created");
			
			PDFWriter pdfWriter2L = new PDFWriter(output21L);
			pdfWriter2L.setType(MZDataConstants.LOSSES);
			pdfWriter2L.write(mzDataCol21);
			pdfWriter2L.close();
			System.out.println("");
			System.out.println("Process finished, the "+output21L+" file for losses was created");
			
			
		}else if(set.getSetName().equals("cset")) {
			String format_input = "mzXML";
			String format_ouput = "cml";
			String file_input = "";
			String file_output = null;
			String sn = "1";
			String mzgap = "0.5";
			String rint = "0.0";
			String acc = "5";
			String occurr = null;
			String rules = "[nitrogenR,RDBER]";
			String e = "[C1..5,H1..15,O0..2,N0..2]";

			OptionData d1 = set.getOption("i");
			for (int i = 0; i < d1.getResultCount(); i++) {
				format_input = d1.getResultDetail(i);
				file_input = d1.getResultValue(i);
			}
			OptionData d2 = set.getOption("o");
			for (int i = 0; i < d2.getResultCount(); i++) {
				format_ouput = d2.getResultDetail(i);
				file_output = d2.getResultValue(i);
			}
			if(set.isSet("sn")){
				OptionData d3 = set.getOption("sn");
				for (int i = 0; i < d3.getResultCount(); i++) {
					sn = d3.getResultValue(i);
				}
			}
			if(set.isSet("mzgap")){
				OptionData d3 = set.getOption("mzgap");
				for (int i = 0; i < d3.getResultCount(); i++) {
					mzgap = d3.getResultValue(i);
				}
			}
			if(set.isSet("rint")){
				OptionData d3 = set.getOption("rint");
				for (int i = 0; i < d3.getResultCount(); i++) {
					rint= d3.getResultValue(i);
				}
			}
			if(set.isSet("acc")){
				OptionData d3 = set.getOption("acc");
				for (int i = 0; i < d3.getResultCount(); i++) {
					acc= d3.getResultValue(i);
				}
			}
			if(set.isSet("rules")){
				OptionData d3 = set.getOption("rules");
				for (int i = 0; i < d3.getResultCount(); i++) {
					rules= d3.getResultValue(i);
				}
			}
			if(set.isSet("occurr")){
				OptionData d3 = set.getOption("occurr");
				for (int i = 0; i < d3.getResultCount(); i++) {
					occurr= d3.getResultValue(i);
				}
			}
			if(set.isSet("e")){
				OptionData d3 = set.getOption("e");
				for (int i = 0; i < d3.getResultCount(); i++) {
					e = d3.getResultValue(i);
				}
			}
			
			System.out.println("Processing the "+file_input+" started ");

			System.out.print(" .");
			if(!(new File(file_input)).exists()){
				String currentdir = System.getProperty("user.dir");
				file_input = currentdir+"/"+file_input;
				if(!(new File(file_input)).exists()){
					System.out.println("The process was interrupted. The "+file_input+" file doesn't exist. ");
					System.exit(1);
				}
			}
			
			if(file_output == null){
				File fi = new File(file_input);
				String name = fi.getName();
				int dotPos = name.lastIndexOf(".");
		        String extension = name.substring(dotPos);
				file_output = fi.getName().replace(extension, "."+format_ouput);
			}else{
				// looking if contains directory path
				if(!file_output.contains("/")){
//					if(file_output.contains("\\\")){
						System.out.println("file_output testing: "+file_output);
						String currentdir = System.getProperty("user.dir");
						file_output = currentdir+"/"+file_output;
				}
			}
			if(!file_output.contains("/")){ // search for current directory
				String currentDir = new File("").getAbsolutePath();
				file_output = currentDir+"/"+file_output;
			}
			
			System.out.print(".");
			
			MZData mzData = new MZData();
			if(format_input.equals("mzXML")){
				InputStream input = new BufferedInputStream(new FileInputStream(file_input));
				System.out.print(".");
				MZXMLReader mzReader = new MZXMLReader(input);
				System.out.print(".");
				mzReader.setMZGap(new Double(mzgap));
				mzReader.setSNThresh(new Double(sn));
				mzReader.setRInt(new Double(rint));
				mzData = mzReader.read(mzData);
				input.close();
				System.out.print(".");
			}else if(format_input.equals("cml")){
				InputStream input = new BufferedInputStream(new FileInputStream(file_input));
				System.out.print(".");
				CMLReader mzReader = new CMLReader(input);
				mzData = mzReader.read(mzData);
				System.out.print(".");
			}else if(format_input.equals("txt")){
				double[][] matrix = MZDataManipulator.getMZDataMatrixFromFile(file_input);
				System.out.print(".");
				mzData = MZDataManipulator.getMZData(matrix, mzData);
				System.out.print(".");
			}
			
			if(set.getData().get(0).equals("convert")){
				if(occurr != null){
					mzData = MZDataManipulator.group(mzData, 0.4);
				}
				if(format_ouput.equals("cml"))
					printCML(mzData, file_output);
				else if (format_ouput.equals("pdf"))
					printPDF(mzData, file_output);
				else if (format_ouput.equals("txt"))
					printTXT(mzData, file_output);
						
			}else if(set.getData().get(0).equals("process")){
				MEFgenerator gen = new MEFgenerator(mzData);
				
				List<Integer> listMA = new ArrayList<Integer>();
				listMA.add(new Integer(acc));
				listMA.add(new Integer(acc));
				listMA.add(new Integer(acc));
				listMA.add(new Integer(acc));
				listMA.add(new Integer(acc));
				gen.setMassAccuracy(listMA);
				
				String[] rulesA = rules.replace("[", "").replace("]", "").split(",");
			 	gen.setRules(rulesA);

				System.out.print(".");
			 	String[] atoms = e.replace("[", "").replace("]", "").split(",");
			 	String[] stringElements = new String[atoms.length];
			 	int[] elementsMax = new int[atoms.length];
			 	int[] elementsMin = new int[atoms.length];
			 	for(int i = 0 ; i < atoms.length ; i++){
			 		String atomT = atoms[i];
			 		String[] maxmin = atomT.split("\\.\\.");
			 		String atom = "";
			 		String min = "";
			 		for (int x = 0; x < maxmin[0].length(); x++) {
			 		    final char c = maxmin[0].charAt(x);
			 		    if ((c >= '0') && (c <= '9')) {
			 		    	// 0 - 9
			 		    	min += c;
			 		    }else{
			 		    	// string
			 		    	atom += c;
			 		    }
			 		  }
			 		stringElements[i] = atom;
			 		elementsMax[i] = new Integer(maxmin[1]);
			 		elementsMin[i] = new Integer(min);
			 	}
			 	MolecularFormulaRange range = MEFgenerator.creatingRange(stringElements, elementsMax, elementsMin);
				gen.setRangeMolecularFormula(range);

				System.out.print(".");
				MZData mzData2 = gen.init();
				System.out.print(".");
				if(occurr != null){
					mzData2 = MZDataManipulator.group(mzData2, Double.parseDouble(occurr));
				}
				
				if(format_ouput.equals("cml"))
					printCML(mzData2, file_output);
				else if (format_ouput.equals("pdf"))
					printPDF(mzData2, file_output);
				else if (format_ouput.equals("txt"))
					printTXT(mzData, file_output);
						
			}
		}
	}

	/**
	 * Copy the mzData into a PDF file
	 * 
	 * @param mzData2
	 * @param fileOutput
	 */
	private static void printPDF(MZData mzData, String fileOutput) {
		System.out.print(".");
		try {
			PDFWriter pdfWriter = new PDFWriter(fileOutput);
			System.out.print(".");
			pdfWriter.write(mzData);
			pdfWriter.close();
		} catch (IOException e){
			e.printStackTrace();
		}
		System.out.println("");
		System.out.println("Process finished, the "+fileOutput+" file was created");
	}

	/**
	 * Copy the mzData into a CML file
	 * 
	 * @param mzData       The MZData object
	 * @param fileOutput   String name of the cml file
	 * @throws IOException
	 */
	private static void printCML(MZData mzData, String fileOutput) throws IOException {
		StringWriter output3 = new StringWriter();
        CMLWriter cmlWriter3 = new CMLWriter(output3);
        cmlWriter3.write(mzData);
        cmlWriter3.close();

		System.out.print(".");
		FileWriter fstream3 = new FileWriter(fileOutput);
		BufferedWriter out3 = new BufferedWriter(fstream3);
		out3.write(output3.toString());
		out3.close();
		fstream3.close();
		System.out.println("");
		System.out.println("Process finished, the "+fileOutput+" file was created");

	}

	/**
	 * Copy the mzData into a txt file
	 * 
	 * @param mzData       The MZData object
	 * @param fileOutput   String name of the txt file
	 * @throws IOException
	 */
	private static void printTXT(MZData mzData, String fileOutput) throws IOException {
		
		System.out.print(".");
		FileWriter fstream3 = new FileWriter(fileOutput);
		BufferedWriter out3 = new BufferedWriter(fstream3);
		out3.write(MZDataManipulator.toString(mzData,false));
		out3.close();
		fstream3.close();
		System.out.println("");
		System.out.println("Process finished, the "+fileOutput+" file was created");

	}

}
