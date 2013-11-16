package org.sams.manipulator;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.Set;

import org.openscience.cdk.formula.IsotopePattern;
import org.openscience.cdk.formula.IsotopePatternGenerator;
import org.openscience.cdk.formula.IsotopePatternSimilarity;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.interfaces.IMolecularFormulaSet;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.sams.MEFgenerator;
import org.sams.MEFgenerator.Polarity;
import org.sams.spect.ParentIon;

public class PrintTools {
	private FileOutputStream fout;
	private static PrintStream print2;
	private PrintStream print;
	private int toleranceISO;
	private double minAbISO;
	private static Polarity polarity;
	public PrintTools(String pathOutput){
		try {
			fout = new FileOutputStream (pathOutput);
			print = new PrintStream(fout);
		} catch (IOException e) {
			e.printStackTrace();
		}	
	}
	public void close(){
	    try {
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}	
	}
	public static void printOutPutCompress(ParentIon ionFinal,PrintStream print) {
		print2 = print;
	}
	/**
	 * Print the output in a tree way. 
	 * 
	 * @param ionFinal The ParentIon
	 */
	public void printOutPutCompress(ParentIon ionFinal, Polarity polari) {
		polarity = polari;
		printingLine(ionFinal);
	    
	    for(int f2 = 0 ; f2 < ionFinal.getFragments().size(); f2++){
        	ParentIon ion2 = ionFinal.getFragments().get(f2);
        	printingLine(ion2);
        	for(int f3 = 0 ; f3 < ion2.getFragments().size(); f3++){
            	ParentIon ion3 = ion2.getFragments().get(f3);
            	printingLine(ion3);
            	for(int f4 = 0 ; f4 < ion3.getFragments().size(); f4++){
                	ParentIon ion4 = ion3.getFragments().get(f4);
                	printingLine(ion4);
                	for(int f5 = 0 ; f5 < ion4.getFragments().size(); f5++){
                    	ParentIon ion5 = ion4.getFragments().get(f5);
                    	printingLine(ion5);
                    	for(int f6 = 0 ; f6 < ion5.getFragments().size(); f6++){
                        	ParentIon ion6 = ion5.getFragments().get(f6);
                        	printingLine(ion6);
                        	for(int f7 = 0 ; f7 < ion6.getFragments().size(); f7++){
                            	ParentIon ion7 = ion6.getFragments().get(f7);
                            	printingLine(ion7);
                        	}
                    	}
                	}
            	}
        	}
	    }
	   
		
	}
	/**
	 * Print the output in a tree way. 
	 * 
	 * @param ionFinal The ParentIon
	 * @param isoP 
	 * @param scoreISO 
	 * @param toleranceISO 
	 * @param polarity 
	 */
	public String printOutPutCompress2(ParentIon ionFinal, IsotopePattern isoP, int toleranceISO, double minAbISO, Polarity polarity) {
		this.toleranceISO = toleranceISO;
		this.minAbISO = minAbISO;
		this.polarity = polarity;
		String output = printingLine(ionFinal,isoP);
	    
	    for(int f2 = 0 ; f2 < ionFinal.getFragments().size(); f2++){
        	ParentIon ion2 = ionFinal.getFragments().get(f2);
        	output += printingLine(ion2);
        	for(int f3 = 0 ; f3 < ion2.getFragments().size(); f3++){
            	ParentIon ion3 = ion2.getFragments().get(f3);
            	output += printingLine(ion3);
            	for(int f4 = 0 ; f4 < ion3.getFragments().size(); f4++){
                	ParentIon ion4 = ion3.getFragments().get(f4);
                	output += printingLine(ion4);
                	for(int f5 = 0 ; f5 < ion4.getFragments().size(); f5++){
                    	ParentIon ion5 = ion4.getFragments().get(f5);
                    	output += printingLine(ion5);
                    	for(int f6 = 0 ; f6 < ion5.getFragments().size(); f6++){
                        	ParentIon ion6 = ion5.getFragments().get(f6);
                        	output += printingLine(ion6);
                        	for(int f7 = 0 ; f7 < ion6.getFragments().size(); f7++){
                            	ParentIon ion7 = ion6.getFragments().get(f7);
                            	output += printingLine(ion7);
                        	}
                    	}
                	}
            	}
        	}
	    }
	   return output;
	}
	public void printThis(String thisPrint){
		print.println(thisPrint);
	}
	private String printingLine(ParentIon ionFinal) {
		return printingLine(ionFinal, null);
	}
	/**
	 * print Line 
	 * @param ionFinal The ParentIon
	 * @param isoP 
	 */
	private String printingLine(ParentIon ionFinal, IsotopePattern isoP) {
		String output = "";
		String form = "";
		String spaceIni = "";
    	double massExp = MEFgenerator.correctingMass(ionFinal.getMass(), polarity);

    	int charge = 0;
		if(polarity.equals(Polarity.positive))
			charge = 1;
		else if(polarity.equals(Polarity.negative))
			charge = -1;
		
		if(ionFinal.getLevel() != 1)
			for(int i = 0; i < ionFinal.getLevel()-1; i++)
				spaceIni += "   ";
		if(ionFinal.getFormulaSet() != null)
		    for(int i = 0; i < ionFinal.getFormulaSet().size(); i++){
		    	IMolecularFormula formula = ionFinal.getFormulaSet().getMolecularFormula(i);
		    	formula.setCharge(charge);
		    	double diff = (massExp-MolecularFormulaManipulator.getTotalExactMass(formula));
		    	
		    	DecimalFormat df1 = new DecimalFormat("####.0");
		    	String result = df1.format(diff*1000000/massExp);
		    	form += MolecularFormulaManipulator.getString(formula)+"<ACC("+result+")";

		    	if(isoP != null){
					IsotopePatternSimilarity is = new IsotopePatternSimilarity();
		    		IsotopePatternGenerator isotopeGe = new IsotopePatternGenerator(minAbISO);
	    			IsotopePattern patternIsoPredicted = isotopeGe.getIsotopes(formula);
		    		is.seTolerance(toleranceISO);
		    		double tempScore = is.compare(isoP, patternIsoPredicted);
		    		form += "/ISO("+tempScore+")";
		    	}
//		    	double iso = ruleISO(ionFinal,formula,ionFinal.getAccuracy());
//		    	df1 = new DecimalFormat("####.00");
//		    	if(iso != -1)
//		    		form += "/ISO("+df1.format(iso)+")";
		    	
		    	if(i != ionFinal.getFormulaSet().size()-1)
		    		form += ", ";
		    }
		else form += " ? ";
	    print.print(spaceIni+ionFinal.getID()+":"+massExp+"("+ionFinal.getIntensity_Abs()+")"+" ["+form+"] ");
	    output += spaceIni+ionFinal.getID()+":"+massExp+"("+ionFinal.getIntensity_Abs()+")"+" ["+form+"] ";
	    IMolecularFormulaSet lossSet = ionFinal.getFormulaLossSet();
//	    System.out.println(lossSet.size());
	    if(ionFinal.getID().equals("1_")){
	    	print.println("");
	    	output += "\n";
	    }else{
	    	print.print(" > ");
	    	output += " > ";
	    	if(lossSet == null){
	    		print.println("");
	    		output += "\n";
	    	}else if(lossSet.size() == 1){
	    		print.println(MolecularFormulaManipulator.getString(lossSet.getMolecularFormula(0)));
	    		output += MolecularFormulaManipulator.getString(lossSet.getMolecularFormula(0))+"\n";
	    	}else{
			    	for(IMolecularFormula formulaLoss : lossSet.molecularFormulas()){
			    		print.print(MolecularFormulaManipulator.getString(formulaLoss) +", ");
			    		output += MolecularFormulaManipulator.getString(formulaLoss) +", ";
			    	}
			    	print.println("");
			    	output += "\n";
			    }
	    }
		return output;
	}
	public void table(HashMap<String, ParentIon> hashMap) {
		String printOut = tableIn(hashMap);
		print.println(printOut);
		
	}
	public String tableIn(HashMap<String, ParentIon> hashMap) {
		String printOut = "\n";
		printOut += "table: \n";
		Set<Entry<String, ParentIon>> set = hashMap.entrySet();
		for(Iterator<Entry<String, ParentIon>> it = set.iterator(); it.hasNext();){
			Entry<String, ParentIon> key = it.next();
			ParentIon ion = key.getValue();
			String idParent = "";
			if(ion.getParent() != null)
				idParent = ion.getParent().getID();
			printOut += ion.getID_XCMS()+" - "+ion.getID()+"<-"+idParent+", M:"+ion.getMass()+", L:"+ion.getLevel()+"\n";
		}
		printOut += "\n";
		return printOut;
	}
}