/*  $RCSfile$
 *  $Author: egonw $
 *  $Date: 2007-09-03 12:53:05 +0200 (Mon, 03 Sep 2007) $
 *  $Revision: 8848 $
 *
 *  Copyright (C) 2005-2007  Miguel Rojas <miguelrojasch@users.sf.net>
 *
 *  Contact: cdk-devel@lists.sourceforge.net
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public License
 *  as published by the Free Software Foundation; either version 2.1
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.sams;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.config.IsotopeFactory;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.formula.IsotopePattern;
import org.openscience.cdk.formula.IsotopePatternGenerator;
import org.openscience.cdk.formula.IsotopePatternSimilarity;
import org.openscience.cdk.formula.MassToFormulaTool;
import org.openscience.cdk.formula.MolecularFormula;
import org.openscience.cdk.formula.MolecularFormulaRange;
import org.openscience.cdk.formula.MolecularFormulaSet;
import org.openscience.cdk.formula.rules.ElementRule;
import org.openscience.cdk.formula.rules.IRule;
import org.openscience.cdk.formula.rules.NitrogenRule;
import org.openscience.cdk.formula.rules.RDBERule;
import org.openscience.cdk.formula.rules.ToleranceRangeRule;
import org.openscience.cdk.interfaces.IIsotope;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.interfaces.IMolecularFormulaSet;
import org.openscience.cdk.nonotify.NNMolecularFormulaSet;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaRangeManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaSetManipulator;
import org.sams.MEFgenerator.Polarity;
import org.sams.spect.ParentIon;

 /**
 *  Class that dereplicate elemental composition of the fragments from a spectral tree
 *  
 * @author Miguel Rojas-Cherto
 */
public class Generator2TreeMF {
		DefaultChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		public boolean flagPrint = false;
		/** Detect if the number of elements is changing*/
		HashMap<String, Integer> mapOccur;
		/** Count the times which produce the same result. 3 times it will be stoped*/
		int countInt = 0;
		/** Maximal level to analyze */
		private int levels;
		/** PrintStream to print the process took account*/
//		private PrintStream out;
		/** Polarity of the experiments. Important for the constraints*/
		private Polarity polarity;
		private IsotopePattern isotoPattern;
		private boolean nitrogenR = false;
		private boolean isotopePatternR = false;
		private boolean RDBER = false;
		private String[] rules;
		private int cyclesMax;
//		private Boolean fixParent;
		private String outPrint = "";
		private double minAbISO;
		private double toleranceISO;
		private double scoreISO;
		
        /**
         * Main class to initiate the process
         * @param scoreISO 
         * @param toleranceMass 
         * @param minAb 
         * @param args
         * @throws IOException 
         * @throws CDKException 
         */
        public Generator2TreeMF(ParentIon compound1, 
        		MolecularFormulaRange mfRangeInitial,
        		double toleranceG, 
        		int levelsMax, 
        		int cyclesMax,
        		Polarity polar, 
        		IsotopePattern spExp,
        		String[] rules, double minAbISO, double toleranceISO, double scoreISO
//        		,Boolean fixParent
        		) {

			try {
				
	        	this.levels = levelsMax;
	        	this.polarity = polar;
	        	this.isotoPattern = spExp;
	        	this.rules = rules;
	        	this.cyclesMax = cyclesMax;
	//        	this.fixParent = fixParent;
	        	this.minAbISO = minAbISO;
	        	this.toleranceISO = toleranceISO;
	        	this.scoreISO = scoreISO;
	        	
	        	boolean flagRedAccur = true;
	        	int countRedAccur = 0;
	
	        	printingINITIAL_COND(mfRangeInitial);
	      		 
	        	while(flagRedAccur){
	        		
	        		////////////////// PRINTING ///////////////////////////////////////////
	        		printing1Step(countRedAccur);
		            printingResultTitle(compound1);
	          		////////////////// PRINTING ///////////////////////////////////////////
	             	  
					countInt = 0;
					mapOccur = new HashMap<String, Integer>();
	
					IMolecularFormulaSet resultsMF1 = null;
					if(countRedAccur == 0){
	//					if(fixParent){ // when is fixed parent the range is the EF
	//						MolecularFormulaRange mfRangeInitialParent = new MolecularFormulaRange();
	//			        	for(IIsotope isotop:mfRangeInitial.isotopes()){
	//			        		mfRangeInitialParent.addIsotope(isotop, mfRangeInitial.getIsotopeCountMax(isotop), mfRangeInitial.getIsotopeCountMax(isotop));
	//			        	}
	//						resultsMF1 = generateMF(compound1,mfRangeInitialParent,compound1.getAccuracy());
	//					}else
								resultsMF1 = generateMF(compound1,mfRangeInitial,compound1.getAccuracy());
					}else
						resultsMF1 = compound1.getFormulaSet();
					
					// Constraints to apply
	            	if(resultsMF1 != null)
	            		resultsMF1 = ruleNIon(resultsMF1); 
	            	if(resultsMF1 != null)
	                	resultsMF1 = ruleRDBE(resultsMF1);
	            	 
	            	compound1.addFormulaSet(resultsMF1);
	            	if(compound1.getFragments().size() == 0){
		            	break;
		            }
	            	
		            IMolecularFormula maximalEle2T = null;
		    		if(resultsMF1 != null){
	
			            MolecularFormulaRange mfRange1 = MolecularFormulaRangeManipulator.getRange(resultsMF1);
				        maximalEle2T = MolecularFormulaRangeManipulator.getMaximalFormula(mfRange1,builder);
				        // TODO: why add to max and min when we have addRange
				        compound1.addMaximalFormula(maximalEle2T);
				        IMolecularFormula minimalEle11 = MolecularFormulaRangeManipulator.getMinimalFormula(mfRange1,builder);
				        compound1.addMinimalFormula(minimalEle11);
	
		                outPrint += "Final Maximal Range = "+MolecularFormulaManipulator.getString(maximalEle2T)+"\n";
			            outPrint += "Final Minimal Range = "+MolecularFormulaManipulator.getString(minimalEle11)+"\n";
					}
		            /*2 level *******************/
		            if(levels >= 2){
		                List<ParentIon> listCompounds2 = compound1.getFragments();
		                for(int i = 0; i < listCompounds2.size(); i++){
		                	ParentIon compound22 = listCompounds2.get(i);
		                	////////////////// PRINTING ///////////////////////////////////////////
	    	            	printingResultTitle(compound22);
	                 		////////////////// PRINTING ///////////////////////////////////////////
		                	IMolecularFormulaSet resultsMF2 = null;
							
		                	if(maximalEle2T != null){
			                	compound22.addMaximalFormula(maximalEle2T);
		                        
		                        //removing those MF are not theirs elements into the range
		                        if(countRedAccur == 0){
									MolecularFormulaRange newRange2 = MolecularFormulaRangeManipulator.getRange(compound1.getFormulaSet());
			                        IMolecularFormula newMax2 = MolecularFormulaRangeManipulator.getMaximalFormula(newRange2,builder);
			                        //TODO: why should it start again from the minimal
			                        IMolecularFormula newMin2 = creatingMinimal0(newMax2);
			                        
			                        IMolecularFormulaSet newFSet2 = new MolecularFormulaSet();
			                        newFSet2.addMolecularFormula(newMax2);
			                        newFSet2.addMolecularFormula(newMin2);
			                        MolecularFormulaRange mfRange2 = MolecularFormulaRangeManipulator.getRange(newFSet2);
		                         
			                        resultsMF2 = generateMF(compound22,mfRange2,compound22.getAccuracy());
		
		        	            	// Constraints to apply
		        	            	if(resultsMF2 != null)
		        	            		resultsMF2 = ruleNIon(resultsMF2);
		        	            	
		        	            }else {
		        	            	if(compound22.getFormulaSet() == null)
		                        		continue;
		                        	
		        	            	IMolecularFormula minimalElea2 = MolecularFormulaRangeManipulator.getMinimalFormula(MolecularFormulaRangeManipulator.getRange(compound22.getFormulaSet()),builder);
		                        	IMolecularFormula maximalElea1 = MolecularFormulaRangeManipulator.getMaximalFormula(MolecularFormulaRangeManipulator.getRange(compound1.getFormulaSet()),builder);
		                        	
		                            IMolecularFormulaSet set = new MolecularFormulaSet();
		                            set.addMolecularFormula(minimalElea2);
		                            set.addMolecularFormula(maximalElea1);
		                            
		                            resultsMF2 = MolecularFormulaSetManipulator.remove(compound22.getFormulaSet(), MolecularFormulaRangeManipulator.getRange(set));;
		                            
		                        }
		                	}
				        	compound22.addFormulaSet(resultsMF2);
				            
				        	IMolecularFormula minimalEle1 = null;
							IMolecularFormula maximalEle3T = null;
							
							if(resultsMF2 != null){
								MolecularFormulaRange mfRange2 = MolecularFormulaRangeManipulator.getRange(resultsMF2);
			
		                        minimalEle1 = MolecularFormulaRangeManipulator.getMinimalFormula(mfRange2,builder);
		                        maximalEle3T = MolecularFormulaRangeManipulator.getMaximalFormula(mfRange2,builder);
		                        
		                        outPrint += "Final Maximal Range = "+MolecularFormulaManipulator.getString(maximalEle3T)+"\n";
					            outPrint += "Final Minimal Range = "+MolecularFormulaManipulator.getString(minimalEle1)+"\n";
					            
		                        compound22.addMaximalFormula(maximalEle3T);
		                        compound22.addMinimalFormula(minimalEle1);
							}                
	                        
							/*3 level *******************/
	    	                if(levels >= 3){
		                        List<ParentIon> listCompounds3 = compound22.getFragments();
		                        for(int j = 0; j < listCompounds3.size(); j++){
		                                ParentIon compound33 = listCompounds3.get(j);
		        	                	////////////////// PRINTING ///////////////////////////////////////////
		            	            	printingResultTitle(compound33);
		                         		////////////////// PRINTING ///////////////////////////////////////////
		                                IMolecularFormulaSet resultsMF3 = null;
		        	    	                        
		                                if(maximalEle3T != null){
			                                compound33.addMaximalFormula(maximalEle3T);
			                                
			                                //removing those MF are not theirs elements into the range
			    	                        if(countRedAccur == 0){
				                                
				                                MolecularFormulaRange newRange3 = MolecularFormulaRangeManipulator.getRange(compound22.getFormulaSet());
				                                IMolecularFormula newMax3 = MolecularFormulaRangeManipulator.getMaximalFormula(newRange3,builder);
						                        IMolecularFormula newMin3 = creatingMinimal0(newMax3);
				                                IMolecularFormulaSet newFSet3 = new MolecularFormulaSet();
				                                newFSet3.addMolecularFormula(newMax3);
				                                newFSet3.addMolecularFormula(newMin3);
				                                
				                                MolecularFormulaRange mfRange3 = MolecularFormulaRangeManipulator.getRange(newFSet3);
				                                
						                        resultsMF3 = generateMF(compound33,mfRange3,compound33.getAccuracy());
		
						                    	// Constraints to apply
						        	            if(resultsMF3 != null)
						        	            	resultsMF3 = ruleNIon(resultsMF3);
				                                
			    	                        }else{
			    	                        	if(compound33.getFormulaSet() == null)
	            	                        		continue;
	            	                        	IMolecularFormula minimalElea3 = MolecularFormulaRangeManipulator.getMinimalFormula(MolecularFormulaRangeManipulator.getRange(compound33.getFormulaSet()),builder);
			    	                        	IMolecularFormula maximalElea2 = MolecularFormulaRangeManipulator.getMaximalFormula(MolecularFormulaRangeManipulator.getRange(compound22.getFormulaSet()),builder);
			    	                        	
				                                IMolecularFormulaSet set = new MolecularFormulaSet();
		                                        set.addMolecularFormula(minimalElea3);
		                                        set.addMolecularFormula(maximalElea2);
		                                        
				                                resultsMF3 = MolecularFormulaSetManipulator.remove(compound33.getFormulaSet(), MolecularFormulaRangeManipulator.getRange(set));;
				                                
			    	                        }
		                                }
			    	                    compound33.addFormulaSet(resultsMF3);
				        	            
			    	                    IMolecularFormula minimalEle2 = null;
		                                IMolecularFormula maximalEle3 = null;
			    	                    if(resultsMF3 != null){
			                                MolecularFormulaRange mfRange3 = MolecularFormulaRangeManipulator.getRange(resultsMF3);
			                                
			                                minimalEle2 = MolecularFormulaRangeManipulator.getMinimalFormula(mfRange3,builder);
			                                maximalEle3 = MolecularFormulaRangeManipulator.getMaximalFormula(mfRange3,builder);
	
			    	                        outPrint += "Final Maximal Range = "+MolecularFormulaManipulator.getString(maximalEle3)+"\n";
			    				            outPrint += "Final Minimal Range = "+MolecularFormulaManipulator.getString(minimalEle2)+"\n";
			    				                       
			                                compound33.addMaximalFormula(maximalEle3);
			                                compound33.addMinimalFormula(minimalEle2);
			    	                    }		                                
		        		                
		                                /*4 level *******************/
		            	                if(levels >= 4){
			                                List<ParentIon> listCompounds4 = compound33.getFragments();
			                                for(int z = 0; z < listCompounds4.size(); z++){
			                                        ParentIon compound44 = listCompounds4.get(z);
			                	                	////////////////// PRINTING ///////////////////////////////////////////
			                    	            	printingResultTitle(compound44);
			                                 		////////////////// PRINTING ///////////////////////////////////////////
			                                        IMolecularFormulaSet resultsMF4 = null;
			                                        
			            	                        if(maximalEle3 != null){
				                                        compound44.addMaximalFormula(maximalEle3);
						                                
				                                         //removing those MF are not theirs elements into the range
				            	                        if(countRedAccur == 0){
					                                        
					                                        
					                                        MolecularFormulaRange newRange4 = MolecularFormulaRangeManipulator.getRange(compound33.getFormulaSet());
					                                        IMolecularFormula newMax4 = MolecularFormulaRangeManipulator.getMaximalFormula(newRange4,builder);
					                                        IMolecularFormula newMin4 = creatingMinimal0(newMax4);
							                                IMolecularFormulaSet newFSet4 = new MolecularFormulaSet();
					                                        newFSet4.addMolecularFormula(newMax4);
					                                        newFSet4.addMolecularFormula(newMin4);
					                                        
					                                        MolecularFormulaRange mfRange4 = MolecularFormulaRangeManipulator.getRange(newFSet4);
							                                
									                        resultsMF4 = generateMF(compound44,mfRange4,compound44.getAccuracy());
					
									                    	// Constraints to apply
									                        if(resultsMF4 != null)
									        	            	resultsMF4 = ruleNIon(resultsMF4);
									                        
				            	                        }else  {
				            	                        	if(compound44.getFormulaSet() == null)
				            	                        		continue;
				            	                        	IMolecularFormula minimalElea4 = MolecularFormulaRangeManipulator.getMinimalFormula(MolecularFormulaRangeManipulator.getRange(compound44.getFormulaSet()),builder);
						    	                        	IMolecularFormula maximalElea3 = MolecularFormulaRangeManipulator.getMaximalFormula(MolecularFormulaRangeManipulator.getRange(compound33.getFormulaSet()),builder);
						    	                        	
							                                IMolecularFormulaSet set = new MolecularFormulaSet();
					                                        set.addMolecularFormula(minimalElea4);
					                                        set.addMolecularFormula(maximalElea3);
					                                        
							                                resultsMF4 = MolecularFormulaSetManipulator.remove(compound44.getFormulaSet(), MolecularFormulaRangeManipulator.getRange(set));;
							                                
						    	                        }
				            	                    }		                                        
			                                        compound44.addFormulaSet(resultsMF4);
			                                        
			                                        IMolecularFormula minimalEle3 = null;
			                                        IMolecularFormula maximalEle4 = null;
			                                        if(resultsMF4 != null){
				                                        MolecularFormulaRange mfRange4 = MolecularFormulaRangeManipulator.getRange(resultsMF4);
				                                        
				                                        minimalEle3 = MolecularFormulaRangeManipulator.getMinimalFormula(mfRange4,builder);
				                                        maximalEle4 = MolecularFormulaRangeManipulator.getMaximalFormula(mfRange4,builder);
	
						    	                        outPrint += "Final Maximal Range = "+MolecularFormulaManipulator.getString(maximalEle4)+"\n";
						    				            outPrint += "Final Minimal Range = "+MolecularFormulaManipulator.getString(minimalEle3)+"\n";
						    				            
				                                        compound44.addMaximalFormula(maximalEle4);
				                                        compound44.addMinimalFormula(minimalEle3);
			                                      
			                                        }
			                                        /*5 level *******************/
			                                    	if(levels >= 5){
						                                List<ParentIon> listCompounds5 = compound44.getFragments();
						                                for(int z5 = 0; z5 < listCompounds5.size(); z5++){
						                                        ParentIon compound55 = listCompounds5.get(z5);
						                	                	////////////////// PRINTING ///////////////////////////////////////////
						                    	            	printingResultTitle(compound55);
						                                 		////////////////// PRINTING ///////////////////////////////////////////
						                                        IMolecularFormulaSet resultsMF5 = null;
						            	                        
						                                        if(maximalEle4 != null){
							                                        compound55.addMaximalFormula(maximalEle4);
		//							                                
		//							                                 //removing those MF are not theirs elements into the range
							            	                        if(countRedAccur == 0){
								                                        
								                                        MolecularFormulaRange newRange5 = MolecularFormulaRangeManipulator.getRange(compound44.getFormulaSet());
								                                        IMolecularFormula newMax5 = MolecularFormulaRangeManipulator.getMaximalFormula(newRange5,builder);
								                                        IMolecularFormula newMin5 = creatingMinimal0(newMax5);
										                                IMolecularFormulaSet newFSet5 = new MolecularFormulaSet();
								                                        newFSet5.addMolecularFormula(newMax5);
								                                        newFSet5.addMolecularFormula(newMin5);
		
								                                        MolecularFormulaRange mfRange5 = MolecularFormulaRangeManipulator.getRange(newFSet5);
										                                
												                        resultsMF5 = generateMF(compound55,mfRange5,compound55.getAccuracy());
								
												                    	// Constraints to apply
												                        if(resultsMF5 != null)
													                        resultsMF5 = ruleNIon(resultsMF5);
												                        
							            	                        }else {
							            	                        	if(compound55.getFormulaSet() == null)
							            	                        		continue;
									    	                        	IMolecularFormula minimalElea5 = MolecularFormulaRangeManipulator.getMinimalFormula(MolecularFormulaRangeManipulator.getRange(compound55.getFormulaSet()),builder);
									    	                        	IMolecularFormula maximalElea4 = MolecularFormulaRangeManipulator.getMaximalFormula(MolecularFormulaRangeManipulator.getRange(compound44.getFormulaSet()),builder);
									    	                        	
										                                IMolecularFormulaSet set = new MolecularFormulaSet();
								                                        set.addMolecularFormula(minimalElea5);
								                                        set.addMolecularFormula(maximalElea4);
								                                        
										                                resultsMF5 = MolecularFormulaSetManipulator.remove(compound55.getFormulaSet(), MolecularFormulaRangeManipulator.getRange(set));;
										                                
									    	                        }
						                                        }
						            	                        compound55.addFormulaSet(resultsMF5);
						            	                        
						            	                        IMolecularFormula minimalEle4 = null;
						                                        IMolecularFormula maximalEle5 = null;
						            	                        if(resultsMF5 != null){
							            	                        MolecularFormulaRange mfRange5 = MolecularFormulaRangeManipulator.getRange(resultsMF5);
							                                        
							                                        minimalEle4 = MolecularFormulaRangeManipulator.getMinimalFormula(mfRange5,builder);
							                                        maximalEle5 = MolecularFormulaRangeManipulator.getMaximalFormula(mfRange5,builder);
	
									    	                        outPrint += "Final Maximal Range = "+MolecularFormulaManipulator.getString(maximalEle4)+"\n";
									    				            outPrint += "Final Minimal Range = "+MolecularFormulaManipulator.getString(minimalEle3)+"\n";
									    				            
							                                        compound55.addMaximalFormula(maximalEle5);
							                                        compound55.addMinimalFormula(minimalEle4);
							                                        
						            	                        }							                                      
						                                        /*6 level *******************/
						                                    	if(levels >= 6){
									                                List<ParentIon> listCompounds6 = compound55.getFragments();
									                                for(int z6 = 0; z6 < listCompounds6.size(); z6++){
									                                        ParentIon compound66 = listCompounds6.get(z6);
									                	                	////////////////// PRINTING ///////////////////////////////////////////
									                    	            	printingResultTitle(compound66);
									                                 		////////////////// PRINTING ///////////////////////////////////////////
									                    	            	IMolecularFormulaSet resultsMF6 = null;
									                    	            	
									                    	            	if(maximalEle5 != null){
										            	                        compound66.addMaximalFormula(maximalEle5);
										            	                        
										                                         //removing those MF are not theirs elements into the range
										            	                        if(countRedAccur == 0){
											                                        
											                                        
											                                        MolecularFormulaRange newRange6 = MolecularFormulaRangeManipulator.getRange(compound55.getFormulaSet());
											                                        IMolecularFormula newMax6 = MolecularFormulaRangeManipulator.getMaximalFormula(newRange6,builder);
											                                        IMolecularFormula newMin6 = creatingMinimal0(newMax6);
													                                IMolecularFormulaSet newFSet6 = new MolecularFormulaSet();
											                                        newFSet6.addMolecularFormula(newMax6);
											                                        newFSet6.addMolecularFormula(newMin6);
											                                        MolecularFormulaRange mfRange6 = MolecularFormulaRangeManipulator.getRange(newFSet6);
													                                
															                        resultsMF6 = generateMF(compound66,mfRange6,compound66.getAccuracy());
											
															                    	// Constraints to apply
															                        if(resultsMF6 != null)
															        	            	resultsMF6 = ruleNIon(resultsMF6);
											                                        									                                        
										            	                        }else   {
	
										            	                        	if(compound66.getFormulaSet() == null)
										            	                        		continue;
										            	                        	
												    	                        	IMolecularFormula minimalElea6 = MolecularFormulaRangeManipulator.getMinimalFormula(MolecularFormulaRangeManipulator.getRange(compound66.getFormulaSet()),builder);
												    	                        	IMolecularFormula maximalElea5 = MolecularFormulaRangeManipulator.getMaximalFormula(MolecularFormulaRangeManipulator.getRange(compound55.getFormulaSet()),builder);
												    	                        	
													                                IMolecularFormulaSet set = new MolecularFormulaSet();
											                                        set.addMolecularFormula(minimalElea6);
											                                        set.addMolecularFormula(maximalElea5);
											                                        
													                                resultsMF6 = MolecularFormulaSetManipulator.remove(compound66.getFormulaSet(), MolecularFormulaRangeManipulator.getRange(set));;
													                                
												    	                        }
									                    	            	}								                                        
									            	                        compound66.addFormulaSet(resultsMF6);
									                                        
									            	                        IMolecularFormula minimalEle5 = null;
									                                        IMolecularFormula maximalEle6 = null;
									            	                        if(resultsMF6 != null){
										            	                        MolecularFormulaRange mfRange6 = MolecularFormulaRangeManipulator.getRange(resultsMF6);
										                                        
										                                        minimalEle5 = MolecularFormulaRangeManipulator.getMinimalFormula(mfRange6,builder);
										                                        maximalEle6 = MolecularFormulaRangeManipulator.getMaximalFormula(mfRange6,builder);
	
												    	                        outPrint += "Final Maximal Range = "+MolecularFormulaManipulator.getString(maximalEle6)+"\n";
												    				            outPrint += "Final Minimal Range = "+MolecularFormulaManipulator.getString(minimalEle5)+"\n";
												    				            
										                                        compound66.addMaximalFormula(maximalEle6);
										                                        compound66.addMinimalFormula(minimalEle5);
									            	                        }										                		                
									                		                /*7 level *******************/
									                                    	if(levels >= 7){
												                                List<ParentIon> listCompounds7 = compound66.getFragments();
												                                for(int z7 = 0; z7 < listCompounds7.size(); z7++){
												                                        ParentIon compound77 = listCompounds7.get(z7);
												                	                	////////////////// PRINTING ///////////////////////////////////////////
												                    	            	printingResultTitle(compound77);
												                                 		////////////////// PRINTING ///////////////////////////////////////////
												                    	            	IMolecularFormulaSet resultsMF7 = null;
												            	                        if(maximalEle6 != null){ 
													                    	            	compound77.addMaximalFormula(maximalEle6);
		//													                                    
													                                        //removing those MF are not theirs elements into the range
													            	                        if(countRedAccur == 0){
													                                       
													                                        
													                                        MolecularFormulaRange newRange7 = MolecularFormulaRangeManipulator.getRange(compound66.getFormulaSet());
													                                        IMolecularFormula newMax7 = MolecularFormulaRangeManipulator.getMaximalFormula(newRange7,builder);
													                                        IMolecularFormula newMin7 = creatingMinimal0(newMax7);
															                                IMolecularFormulaSet newFSet7 = new MolecularFormulaSet();
													                                        newFSet7.addMolecularFormula(newMax7);
													                                        newFSet7.addMolecularFormula(newMin7);
		
													                                        MolecularFormulaRange mfRange7 = MolecularFormulaRangeManipulator.getRange(newFSet7);
															                                															                        
																	                        resultsMF7 = generateMF(compound77,mfRange7,compound77.getAccuracy());
													
																	                    	// Constraints to apply
																	                        if(resultsMF7 != null)
																	        	            	resultsMF7 = ruleNIon(resultsMF7);
													                                        
													            	                        }else{
													            	                        	if(compound77.getFormulaSet() == null)
													            	                        		continue;
															    	                        	IMolecularFormula minimalElea7 = MolecularFormulaRangeManipulator.getMinimalFormula(MolecularFormulaRangeManipulator.getRange(compound77.getFormulaSet()),builder);
															    	                        	IMolecularFormula maximalElea6 = MolecularFormulaRangeManipulator.getMaximalFormula(MolecularFormulaRangeManipulator.getRange(compound66.getFormulaSet()),builder);
															    	                        	
																                                IMolecularFormulaSet set = new MolecularFormulaSet();
														                                        set.addMolecularFormula(minimalElea7);
														                                        set.addMolecularFormula(maximalElea6);
														                                        
																                                resultsMF7 = MolecularFormulaSetManipulator.remove(compound77.getFormulaSet(), MolecularFormulaRangeManipulator.getRange(set));;
																                                
															    	                        }
												            	                        }											                                        
												            	                        compound77.addFormulaSet(resultsMF7);
												                                        
												            	                        if(resultsMF7 != null){
													            	                        MolecularFormulaRange mfRange7 = MolecularFormulaRangeManipulator.getRange(resultsMF7);
													                                        
													                                        IMolecularFormula minimalEle6 = MolecularFormulaRangeManipulator.getMinimalFormula(mfRange7,builder);
													                                        IMolecularFormula maximalEle7 = MolecularFormulaRangeManipulator.getMaximalFormula(mfRange7,builder);
	
															    	                        outPrint += "Final Maximal Range = "+MolecularFormulaManipulator.getString(maximalEle7)+"\n";
															    				            outPrint += "Final Minimal Range = "+MolecularFormulaManipulator.getString(minimalEle6)+"\n";
															    				            
													                                        compound77.addMaximalFormula(maximalEle7);
													                                        compound77.addMinimalFormula(minimalEle6);
												            	                        }												                		                
												                                }//listCompounds77
											            	                }//levels[7]	
									                		                
									                                }//listCompounds66
								            	                }//levels[6]	
						                                        
						                                }//listCompounds55
					            	                }//levels[5]		
			                                        
			                                }//listCompounds44
		            	                }//levels[4]		                                
		                        }//listCompounds33
	    	                }// levels[3]
	                }//listCompounds22
	                
	            }//levels[2]
	                
	            outPrint += "\n";
	            outPrint += "#########################################\n";
	            outPrint += "#\n";
	            outPrint += "#  2er. STEP: Generate EF for all losses\n";
	            outPrint += "#\n";
	            outPrint += "#########################################\n";
	            outPrint += "\n";
	            
	            if(levels >= 2){ 
	                for(int i = 0; i < compound1.getFragments().size(); i++){
	                    ParentIon compound22 = compound1.getFragments().get(i);
	            		//////////////////PRINTING ///////////////////////////////////////////
	                	printingResultTitleLoss(compound1,compound22);
	             		////////////////// PRINTING ///////////////////////////////////////////
	                    
	                	if(compound22.getFormulaSet() == null){
	                		outPrint += "Non containg MF\n";
	                		continue;
	                	}if(compound1.getFormulaSet() == null){
	                		compound22.addFormulaSet(null);
	                	}else
	                		applyingLossConstraints(compound1,compound22);
	                    
	                    // 3 ---------------
		                if(levels >= 3){
	                        List<ParentIon> ListCompounds3 = compound22.getFragments();
	                        for(int j = 0; j < ListCompounds3.size(); j++){
	                            ParentIon compound33 = ListCompounds3.get(j);
	                            //////////////////PRINTING ///////////////////////////////////////////
	                            printingResultTitleLoss(compound22,compound33);
	                     		////////////////// PRINTING ///////////////////////////////////////////
	                        	if(compound33.getFormulaSet() == null){
	                        		outPrint += "Non containg MF\n";
	                        		continue;
	                        	}else if(compound22.getFormulaSet() == null){
	                        		compound33.addFormulaSet(null);
	                        	}else
	                        		applyingLossConstraints(compound22,compound33);
	                            
	                       	   // 4 ---------------
		    	                if(levels >= 4){
		                                List<ParentIon> ListCompounds4 = compound33.getFragments();
		                                for(int j4 = 0; j4 < ListCompounds4.size(); j4++){
		                                        
		                                    ParentIon compound44 = ListCompounds4.get(j4);
		                                    //////////////////PRINTING ///////////////////////////////////////////
		                                    printingResultTitleLoss(compound33,compound44);
		                             		////////////////// PRINTING ///////////////////////////////////////////
		                                    
		                                	if(compound44.getFormulaSet() == null){
		                                		outPrint += "Non containg MF\n";
		                                		continue;
		                                	}else if(compound33.getFormulaSet() == null){
		                                		compound44.addFormulaSet(null);
		                                	}else
		                                		applyingLossConstraints(compound33,compound44);
		                                    
		                               	   // 5 ---------------
		                                	if(levels >= 5){
					                                List<ParentIon> ListCompounds5 = compound44.getFragments();
					                                for(int j5 = 0; j5 < ListCompounds5.size(); j5++){
					                                        
					                                    ParentIon compound55 = ListCompounds5.get(j5);
					                                    //////////////////PRINTING ///////////////////////////////////////////
					                                    printingResultTitleLoss(compound44,compound55);
					                             		////////////////// PRINTING ///////////////////////////////////////////
					                                    
					                                	if(compound55.getFormulaSet() == null){
					                                		outPrint += "Non containg MF\n";
					                                		continue;
					                                	}else if(compound44.getFormulaSet() == null){
					                                		compound55.addFormulaSet(null);
					                                	}else
					                                		applyingLossConstraints(compound44,compound55);
					                                    
						                                 // 6 ---------------
						                                if(levels >= 6){
								                                List<ParentIon> ListCompounds6 = compound55.getFragments();
								                                for(int j6 = 0; j6 < ListCompounds6.size(); j6++){
								                                        
								                                    ParentIon compound66 = ListCompounds6.get(j6);
								                                    //////////////////PRINTING ///////////////////////////////////////////
								                                    printingResultTitleLoss(compound55,compound66);
								                             		////////////////// PRINTING ///////////////////////////////////////////
								                                    
								                                	if(compound66.getFormulaSet() == null){
								                                		outPrint += "Non containg MF\n";
								                                		continue;
								                                	}else if(compound55.getFormulaSet() == null){
								                                		compound66.addFormulaSet(null);
								                                	}else
								                                		applyingLossConstraints(compound55,compound66);
								                                    			                                    
								                                    // 7 ---------------
								                                	if(levels >= 7){
											                                List<ParentIon> ListCompounds7 = compound66.getFragments();
											                                for(int j7 = 0; j7 < ListCompounds7.size(); j7++){
											                                        
											                                    ParentIon compound77 = ListCompounds7.get(j7);
											                                    //////////////////PRINTING ///////////////////////////////////////////
											                                    printingResultTitleLoss(compound66,compound77);
											                             		////////////////// PRINTING ///////////////////////////////////////////
											                                    
											                                	if(compound77.getFormulaSet() == null){
											                                		outPrint += "Non containg MF\n";
											                                		continue;
											                                	}else if(compound66.getFormulaSet() == null){
											                                		compound77.addFormulaSet(null);
											                                	}else
											                                		applyingLossConstraints(compound66,compound77);
											                                                                   	   
											                                }// compound7.size
								                                	}//level[7]
								                                	
									                         }// compound6.size
									    	            }//level[6]
						                                
					                                }// compound5.size
					    	                }//level[5]
		                                	
		                                }// compound4.size
		    	                }//level[4]
	                        
	                        }// compound3.size
		                }//level[3]
	                    
	                }// compound2.size
	            }// levels[2]
	                
	
	            
	            boolean flag = true;
	            /** Count the number of cycles realized */
	            int count = 1;
	            while(flag){
	            	
	        		outPrint += "\n";
	        		outPrint += "////////////////////////////////////////////////////////////////////\n";
	        		outPrint += "// \n";
	        		outPrint += "// 3er. STEP: Reduction of the fragments\n";
	        		outPrint += "// \n";
	        		outPrint += "// cycle("+(count)+")\n";
	        		outPrint += "// \n";
	        		outPrint += "////////////////////////////////////////////////////////////////////\n";
	        		outPrint += "\n";
	        		count++;
	        		
	        		//////////////////PRINTING ///////////////////////////////////////////
	            	printingResultTitle(compound1);
	                outPrint += "---------------------------------------------\n";
	                outPrint += "\n";
	         		////////////////// PRINTING ///////////////////////////////////////////
	            	
	                // looking for the minimal range which should be the minimal of all fragments
	                IMolecularFormula formulaMin11 = null;
	                if(levels == 1)
	                	formulaMin11 = compound1.getMinimalMF();
	                else{
	                	formulaMin11 = getMinimalMF(compound1);
	                }
	                
	                IMolecularFormula formulaMax1 = compound1.getMaximalMF();
	                if(formulaMin11 == null){
	                	outPrint += "Not containg Fragments\n";
	                	break;
	                }
	                IMolecularFormula formulaMin1 = formulaMin11;
	                MolecularFormulaRange mfRangeNew1 = new MolecularFormulaRange();
	                Iterator<IIsotope> ii = formulaMax1.isotopes().iterator();
	                while(ii.hasNext()){
	                        IIsotope isoto = ii.next();
	                        mfRangeNew1.addIsotope(isoto, formulaMin1.getIsotopeCount(isoto), formulaMax1.getIsotopeCount(isoto));
	                }
	                if(compound1.getFormulaSet() == null){
	                	outPrint += "No containing elemental formulas\n";
	                	break;
	                }
	                IMolecularFormulaSet newMFC1 = MolecularFormulaSetManipulator.remove(compound1.getFormulaSet(), mfRangeNew1);	
	                
	                ///////////////////////////////////////////////////
	               	int diff = compound1.getFormulaSet().size() - newMFC1.size();
	               	outPrint += "MF removed "+diff+", from "+compound1.getFormulaSet().size()+"\n";
	               	outPrint += "MF remaining "+newMFC1.size()+"\n";
	               	printingResultInitialCycle(newMFC1, mfRangeNew1);
	               	///////////////////////////////////////////////////
	               	                
	                compound1.addFormulaSet(newMFC1);
	                
	                MolecularFormulaRange rangeExtract = MolecularFormulaRangeManipulator.getRange(newMFC1);
	
	                IMolecularFormula formulaMin1E = MolecularFormulaRangeManipulator.getMinimalFormula(rangeExtract,builder);
	                IMolecularFormula formulaMax1E = MolecularFormulaRangeManipulator.getMaximalFormula(rangeExtract,builder);
	                compound1.addMaximalFormula(formulaMax1E);
	                compound1.addMinimalFormula(formulaMin1E);
	                
	                if(levels >= 2){
	                    List<ParentIon> listCompounds2 = compound1.getFragments();
	                    for(int i = 0; i < listCompounds2.size(); i++){
	                    	ParentIon compound22 = listCompounds2.get(i);
	                		//////////////////PRINTING ///////////////////////////////////////////
	                    	printingResultTitle(compound22);
	                 		////////////////// PRINTING ///////////////////////////////////////////
	                        
	                    	if(compound22.getFormulaSet() == null)
	                    		continue;
	                    	
	                        IMolecularFormula formulaMax2 = compound22.getMaximalMF();
	                        IMolecularFormula formulaMin2 = getMinimalMF(compound22);
	                        
	                        if(formulaMin2 == null){
	                        	formulaMin2 = compound22.getMinimalMF();
	                        }
	                        
	                        MolecularFormulaRange mfRangeNew2 = new MolecularFormulaRange();
	                        Iterator<IIsotope> ii2 = formulaMax2.isotopes().iterator();
	                        while(ii2.hasNext()){
	                                IIsotope isoto2 = ii2.next();
	                                mfRangeNew2.addIsotope(isoto2, formulaMin2.getIsotopeCount(isoto2), formulaMax2.getIsotopeCount(isoto2));
	                        }
	                                                
	                        IMolecularFormulaSet newMFC2 = MolecularFormulaSetManipulator.remove(compound22.getFormulaSet(), mfRangeNew2);
	
	                       	///////////////////////////////////////////////////
	                        diff = compound22.getFormulaSet().size() - newMFC2.size();
	                       	outPrint += "MF removed "+diff+", from "+compound22.getFormulaSet().size()+"\n";
	                       	outPrint += "MF remaining "+newMFC2.size()+"\n";
	                       	printingResultInitialCycle(newMFC2, mfRangeNew2);
	                       	///////////////////////////////////////////////////
	                       	
	                        compound22.addFormulaSet(newMFC2);
	                        
	                        MolecularFormulaRange rangeExtract2 = MolecularFormulaRangeManipulator.getRange(newMFC2);
	                        IMolecularFormula formulaMinE2 = MolecularFormulaRangeManipulator.getMinimalFormula(rangeExtract2,builder);
	                        IMolecularFormula formulaMaxE2 = MolecularFormulaRangeManipulator.getMaximalFormula(rangeExtract2,builder);
	                        compound22.addMinimalFormula(formulaMinE2);
	                        compound22.addMaximalFormula(formulaMaxE2);
	
	                        // constraining loss
	                        compound22.addFormulaLossSet(fittingRemovingLoss(compound1, compound22, compound22.getFormulaLossSet()));
	                        // constraining fragment
	                        compound1 = fittingRemovingParent1(compound1, compound22, compound22.getFormulaLossSet());
	                        
	    	                if(levels >= 3){
	                            List<ParentIon> ListCompounds3 = compound22.getFragments();
	                            for(int j = 0; j < ListCompounds3.size(); j++){
	                                    ParentIon compound33 = ListCompounds3.get(j);
	                                    //////////////////PRINTING ///////////////////////////////////////////
	                                	printingResultTitle(compound33);
	                             		////////////////// PRINTING ///////////////////////////////////////////
	                                    
	                                	if(compound33.getFormulaSet() == null)
	                                		continue;
	                                	
	                                	IMolecularFormula formulaMax3 =	compound33.getMaximalMF();
	                                    IMolecularFormula formulaMin3 = compound33.getMinimalMF();
	                                    
	                                    MolecularFormulaRange mfRangeNew3 = new MolecularFormulaRange();
	                                    Iterator<IIsotope> ii3 = formulaMax3.isotopes().iterator();
	                                    while(ii3.hasNext()){
	                                            IIsotope isoto3 = ii3.next();
	                                            mfRangeNew3.addIsotope(isoto3, formulaMin3.getIsotopeCount(isoto3), formulaMax3.getIsotopeCount(isoto3));
	                                    }
	                                    
	                                    IMolecularFormulaSet newMFC3 = MolecularFormulaSetManipulator.remove(compound33.getFormulaSet(), mfRangeNew3);
	
	                                   	///////////////////////////////////////////////////
	                                    diff = compound33.getFormulaSet().size() - newMFC3.size();
	                                   	outPrint += "MF removed "+diff+", from "+compound33.getFormulaSet().size()+"\n";
	                                   	outPrint += "MF remaining "+newMFC3.size()+"\n";
	                                   	printingResultInitialCycle(newMFC3, mfRangeNew3);
	                                   	///////////////////////////////////////////////////
	                                   	
	                                    compound33.addFormulaSet(newMFC3);
	                                    
	                                    MolecularFormulaRange rangeExtract3 = MolecularFormulaRangeManipulator.getRange(newMFC3);
	                                    IMolecularFormula formulaMinE3 = MolecularFormulaRangeManipulator.getMinimalFormula(rangeExtract3,builder);
	                                    IMolecularFormula formulaMaxE3 = MolecularFormulaRangeManipulator.getMaximalFormula(rangeExtract3,builder);
	                                    compound33.addMaximalFormula(formulaMaxE3);
	                                    compound33.addMinimalFormula(formulaMinE3);
	                                    
	                                    // constraining loss
	        	                        compound33.addFormulaLossSet(fittingRemovingLoss(compound22, compound33, compound33.getFormulaLossSet()));
	        	                        // constraining fragment
	        	                        compound22 = fittingRemovingParent1(compound22, compound33, compound33.getFormulaLossSet());
	        	                        
		            	                if(levels >= 4){
	                                        List<ParentIon> ListCompounds4 = compound33.getFragments();
	                                        for(int j4 = 0; j4 < ListCompounds4.size(); j4++){
	                                                ParentIon compound44 = ListCompounds4.get(j4);
	                                                //////////////////PRINTING ///////////////////////////////////////////
	                                            	printingResultTitle(compound44);
	                                         		////////////////// PRINTING ///////////////////////////////////////////
	                                                
	                                            	if(compound44.getFormulaSet() == null)
	                                            		continue;
	                                            	
	                                            	IMolecularFormula formulaMax4 = compound44.getMaximalMF();
	                                                IMolecularFormula formulaMin4 = compound44.getMinimalMF();
	                                                
	                                                MolecularFormulaRange mfRangeNew4 = new MolecularFormulaRange();
	                                                Iterator<IIsotope> ii4 = formulaMax4.isotopes().iterator();
	                                                while(ii4.hasNext()){
	                                                        IIsotope isoto4 = ii4.next();
	                                                        mfRangeNew4.addIsotope(isoto4, formulaMin4.getIsotopeCount(isoto4), formulaMax4.getIsotopeCount(isoto4));
	                                                }
	                                                IMolecularFormulaSet newMFC4 = MolecularFormulaSetManipulator.remove(compound44.getFormulaSet(), mfRangeNew4);
	
	                                               	///////////////////////////////////////////////////
	                                                diff = compound44.getFormulaSet().size() - newMFC4.size();
	                                               	outPrint += "MF removed "+diff+", from "+compound44.getFormulaSet().size()+"\n";
	                                               	outPrint += "MF remaining "+newMFC4.size()+"\n";
	                                               	printingResultInitialCycle(newMFC4, mfRangeNew4);
	                                               	///////////////////////////////////////////////////
	                                               	
	                                                compound44.addFormulaSet(newMFC4);
	                                                
	                                                MolecularFormulaRange rangeExtract4 = MolecularFormulaRangeManipulator.getRange(newMFC4);
	                                                IMolecularFormula formulaMinE4 = MolecularFormulaRangeManipulator.getMinimalFormula(rangeExtract4,builder);
	                                                IMolecularFormula formulaMaxE4 = MolecularFormulaRangeManipulator.getMaximalFormula(rangeExtract4,builder);
	                                                compound44.addMaximalFormula(formulaMaxE4);
	                                                compound44.addMinimalFormula(formulaMinE4);

	        	                                    // constraining loss
	        	        	                        compound44.addFormulaLossSet(fittingRemovingLoss(compound33, compound44, compound44.getFormulaLossSet()));
	        	        	                        // constraining fragment
	        	        	                        compound33 = fittingRemovingParent1(compound33, compound44, compound44.getFormulaLossSet());
	        	        	                        
	                                                if(levels >= 5){
				                                        List<ParentIon> ListCompounds5 = compound44.getFragments();
				                                        for(int j5 = 0; j5 < ListCompounds5.size(); j5++){
				                                                ParentIon compound55 = ListCompounds5.get(j5);
				                                                //////////////////PRINTING ///////////////////////////////////////////
				                                            	printingResultTitle(compound55);
				                                         		////////////////// PRINTING ///////////////////////////////////////////
				                                                
				                                            	if(compound55.getFormulaSet() == null)
				                                            		continue;
				                                            	
				                                            	IMolecularFormula formulaMax5 = compound55.getMaximalMF();
				                                                IMolecularFormula formulaMin5 = compound55.getMinimalMF();
				                                                
				                                                MolecularFormulaRange mfRangeNew5 = new MolecularFormulaRange();
				                                                Iterator<IIsotope> ii5 = formulaMax5.isotopes().iterator();
				                                                while(ii5.hasNext()){
				                                                        IIsotope isoto5 = ii5.next();
				                                                        mfRangeNew5.addIsotope(isoto5, formulaMin5.getIsotopeCount(isoto5), formulaMax5.getIsotopeCount(isoto5));
				                                                }
				                                                IMolecularFormulaSet newMFC5 = MolecularFormulaSetManipulator.remove(compound55.getFormulaSet(), mfRangeNew5);
				                                                compound55.addFormulaSet(newMFC5);
	
				                                               	///////////////////////////////////////////////////
				                                                diff = compound55.getFormulaSet().size() - newMFC5.size();
				                                               	outPrint += "MF removed "+diff+", from "+compound55.getFormulaSet().size()+"\n";
				                                               	outPrint += "MF remaining "+newMFC5.size()+"\n";
				                                               	printingResultInitialCycle(newMFC5, mfRangeNew5);
				                                               	///////////////////////////////////////////////////
				                                               	
				                                                MolecularFormulaRange rangeExtract5 = MolecularFormulaRangeManipulator.getRange(newMFC5);
				                                                IMolecularFormula formulaMinE5 = MolecularFormulaRangeManipulator.getMinimalFormula(rangeExtract5,builder);
				                                                IMolecularFormula formulaMaxE5 = MolecularFormulaRangeManipulator.getMaximalFormula(rangeExtract5,builder);
				                                                compound55.addMaximalFormula(formulaMaxE5);
				                                                compound55.addMinimalFormula(formulaMinE5);

				        	                                    // constraining loss
				        	        	                        compound55.addFormulaLossSet(fittingRemovingLoss(compound44, compound55, compound55.getFormulaLossSet()));
				        	        	                        // constraining fragment
				        	        	                        compound44 = fittingRemovingParent1(compound44, compound55, compound55.getFormulaLossSet());
				        	        	                        
				                                                if(levels >= 6){
							                                        List<ParentIon> ListCompounds6 = compound55.getFragments();
							                                        for(int j6 = 0; j6 < ListCompounds6.size(); j6++){
							                                                ParentIon compound66 = ListCompounds6.get(j6);
							                                                //////////////////PRINTING ///////////////////////////////////////////
							                                            	printingResultTitle(compound66);
							                                         		////////////////// PRINTING ///////////////////////////////////////////
							                                                
							                                            	if(compound66.getFormulaSet() == null)
							                                            		continue;
							                                            	
							                                            	IMolecularFormula formulaMax6 = compound66.getMaximalMF();
							                                                IMolecularFormula formulaMin6 = compound66.getMinimalMF();
							                                                
							                                                MolecularFormulaRange mfRangeNew6 = new MolecularFormulaRange();
							                                                Iterator<IIsotope> ii6 = formulaMax6.isotopes().iterator();
							                                                while(ii6.hasNext()){
							                                                        IIsotope isoto6 = ii6.next();
							                                                        mfRangeNew6.addIsotope(isoto6, formulaMin6.getIsotopeCount(isoto6), formulaMax6.getIsotopeCount(isoto6));
							                                                }
							                                                
							                                                IMolecularFormulaSet newMFC6 = MolecularFormulaSetManipulator.remove(compound66.getFormulaSet(), mfRangeNew6);
	
							                                               	///////////////////////////////////////////////////
							                                                diff = compound66.getFormulaSet().size() - newMFC6.size();
							                                               	outPrint += "MF removed "+diff+", from "+compound66.getFormulaSet().size()+"\n";
							                                               	outPrint += "MF remaining "+newMFC6.size()+"\n";
							                                               	printingResultInitialCycle(newMFC6, mfRangeNew6);
							                                               	///////////////////////////////////////////////////
							                                               	
							                                                compound66.addFormulaSet(newMFC6);
							                                                MolecularFormulaRange rangeExtract6 = MolecularFormulaRangeManipulator.getRange(newMFC6);
							
							                                                IMolecularFormula formulaMinE6 = MolecularFormulaRangeManipulator.getMinimalFormula(rangeExtract6,builder);
							                                                IMolecularFormula formulaMaxE6 = MolecularFormulaRangeManipulator.getMaximalFormula(rangeExtract6,builder);
							                                                compound66.addMaximalFormula(formulaMaxE6);
							                                                compound66.addMinimalFormula(formulaMinE6);

							        	                                    // constraining loss
							        	        	                        compound66.addFormulaLossSet(fittingRemovingLoss(compound55, compound66, compound66.getFormulaLossSet()));
							        	        	                        // constraining fragment
							        	        	                        compound55 = fittingRemovingParent1(compound55, compound66, compound66.getFormulaLossSet());
							        	        	                        
							                                                if(levels >= 7){
										                                        List<ParentIon> ListCompounds7 = compound66.getFragments();
										                                        for(int j7 = 0; j7 < ListCompounds7.size(); j7++){
										                                                ParentIon compound77 = ListCompounds7.get(j7);
										                                                //////////////////PRINTING ///////////////////////////////////////////
										                                            	printingResultTitle(compound77);
										                                         		////////////////// PRINTING ///////////////////////////////////////////
										                                                
										                                            	if(compound77.getFormulaSet() == null)
										                                            		continue;
										                                            	
										                                            	//may be can be smaller but also bigger
										                                                IMolecularFormula formulaMax7 = compound77.getMaximalMF();
										                                                IMolecularFormula formulaMin7 = compound77.getMinimalMF();
										                                                
										                                                MolecularFormulaRange mfRangeNew7 = new MolecularFormulaRange();
										                                                Iterator<IIsotope> ii7 = formulaMax7.isotopes().iterator();
										                                                while(ii7.hasNext()){
										                                                        IIsotope isoto7 = ii7.next();
										                                                        mfRangeNew7.addIsotope(isoto7, formulaMin7.getIsotopeCount(isoto7), formulaMax7.getIsotopeCount(isoto7));
										                                                }
										                                                
										                                                IMolecularFormulaSet newMFC7 = MolecularFormulaSetManipulator.remove(compound77.getFormulaSet(), mfRangeNew7);
	
										                                                diff = compound77.getFormulaSet().size() - newMFC7.size();
										                                               	outPrint += "MF removed "+diff+", from "+compound77.getFormulaSet().size()+"\n";
										                                               	outPrint += "MF remaining "+newMFC7.size()+"\n";
										                                               	printingResultInitialCycle(newMFC7, mfRangeNew7);
										                                               	///////////////////////////////////////////////////
										                                               	
										                                                compound66.addFormulaSet(newMFC7);
										                                                MolecularFormulaRange rangeExtract7 = MolecularFormulaRangeManipulator.getRange(newMFC7);
										
										                                                IMolecularFormula formulaMinE7 = MolecularFormulaRangeManipulator.getMinimalFormula(rangeExtract7,builder);
										                                                IMolecularFormula formulaMaxE7 = MolecularFormulaRangeManipulator.getMaximalFormula(rangeExtract7,builder);
	
										                                                outPrint += "Final Maximal MF Post = "+MolecularFormulaManipulator.getString(formulaMaxE7)+"\n";
										                                                outPrint += "Final Minimal MF Post = "+MolecularFormulaManipulator.getString(formulaMinE7)+"\n";
										                                                
										                                                compound77.addMaximalFormula(formulaMaxE7);
										                                                compound77.addMinimalFormula(formulaMinE7);

										        	                                    // constraining loss
										        	        	                        compound77.addFormulaLossSet(fittingRemovingLoss(compound66, compound77, compound77.getFormulaLossSet()));
										        	        	                        // constraining fragment
										        	        	                        compound66 = fittingRemovingParent1(compound66, compound77, compound77.getFormulaLossSet());
										        	        	                        
										                                        }//listCompounds7.size()
											            	                }//levels[7]
							                                                
							                                        }//listCompounds6.size()
								            	                }//levels[6]
				                                                
				                                        }//listCompounds5.size()
					            	                }//levels[5]
	                                        }//listCompounds4.size()
		            	                }//levels[4]
		            	                
	                            }//listCompounds3.size()
	    	                }//levels[3]
	                        
	                }//listCompounds2.size()
	            }//levels[2]
	            
	            
	                boolean stopFlagF = isMonoResult(compound1);
	                if(stopFlagF){
	//                	printResultFinal(compound1,count);
	                	break;
	                }
	                
	                boolean stopFlag = isFinishCycle2(compound1);
	                
	                if((count == cyclesMax)|(stopFlag)){
	                	
	                	printResultFinal(compound1,count);
	                	break;
	                }
	                    	
	                    
	            }// cycles
	            countRedAccur++;
	            if(countRedAccur == 2)
	            	if(isotopePatternR){
	            	if(resultsMF1 != null)
	            		if(spExp != null){
	            			compound1.addFormulaSet(ruleIsotopePattern(compound1.getFormulaSet()));
	            		}
		           
	            }
	            if(countRedAccur == 3){
	            	outPrint += "\nBreaking because not more variation\n";
	//                PrintTools.printOutPutCompress(compound1,out);
	                break;
	            }
	            
	            boolean stopFlagFinal = isFinishCycle1(compound1);
	            if(stopFlagFinal){
	//            	PrintTools.printOutPutCompress(compound1,out);
	                outPrint += "\nBreaking because fitted in one all not losses\n";
		            break;
	            	
	            }
	        }// while(flagRedAccur)
		} catch (IOException e) {
			System.err.println("IOException");
			System.out.println(outPrint);
			
		}catch (CDKException e) {
			System.err.println("CDKException");
			System.out.println(outPrint);
		}
    }
        private void printing1Step(int countRedAccur) {
        	 outPrint += "\n";
    		 outPrint += "\n";
    		 outPrint += "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
    		 outPrint += "%\n";
    		 outPrint += "%  1er. STEP: Generate EF for all fragments\n";
    		 outPrint += "%\n";
    		 outPrint += "%  Accuracy = "+countRedAccur+", Cycle "+countRedAccur+"\n";
    		 outPrint += "%\n";
    		 outPrint += "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
    		 outPrint += "\n";
		}
		private void printingINITIAL_COND(MolecularFormulaRange mfRangei) throws FileNotFoundException {
        	outPrint += "\n";
	   		outPrint += "\n";
	   		outPrint += "////////////////////////////////////////////////////\n";
	   		outPrint += "////////////////////////////////////////////////////\n";
	  		outPrint += "//                                                \n";
	   		outPrint += "//                  NEW TREE                      \n";
	   		outPrint += "//                                                \n";
	   		outPrint += "//  Steps followed:\n";
	   		outPrint += "//  1: Generate all fragments\n";
	   		outPrint += "//    1.1: Validate fragments according the rules\n";
	   		outPrint += "//  2: Generate Losses\n";
	   		outPrint += "//    2.1: Validate Losses according the rules\n";
	   		outPrint += "//    2.2: Validate relation sum: Parent = fragment + loss\n";
	   		outPrint += "//  3: Reduction fragmetns according relation Parent/children\n";
	   		outPrint += "//\n";
	   		outPrint += "////////////////////////////////////////////////////\n";
	   		outPrint += "////////////////////////////////////////////////////\n";
	  		outPrint += "\n";
	  		outPrint += "INITIAL CONDITIONS\n";
        	outPrint += "Max Levels: "+levels+"\n";
        	outPrint += "Polarity  : "+polarity+"\n";
        	outPrint += "Cycle Max.: "+cyclesMax+"\n";
        	outPrint += "ISO patern: "+isotoPattern+"\n";
        	if(isotoPattern != null){
        		for(int i = 0 ; i < isotoPattern.getNumberOfIsotopes(); i++)
        			outPrint += isotoPattern.getIsotope(i).getMass()+":"+isotoPattern.getIsotope(i).getIntensity()+ " ";
        		outPrint += "\n";
        	}
			outPrint += "RULES     : "+rules.length+"\n";
        	for(int i = 0; i < rules.length; i++){
        		if(rules[i].equals("nitrogenR")){
        			outPrint += "   nitrogenR\n";
        			nitrogenR = true;
        		}else if(rules[i].equals("isotopePatternR")){
        			isotopePatternR = true;
        			outPrint += "   isotopePatternR\n";
        		}else if(rules[i].equals("RDBER")){
        			outPrint += "   RDBER\n";
        			RDBER = true;
        		}
        	}
//        	outPrint += "Fix Parent: "+fixParent+"\n";
        	String ima = "";
        	String imi = "";
//        	if(fixParent){ // when is fixed parent the range is the EF
//				MolecularFormulaRange mfRangeInitialParent = new MolecularFormulaRange();
//	        	for(IIsotope isotop:mfRangei.isotopes()){
//	        		mfRangeInitialParent.addIsotope(isotop, mfRangei.getIsotopeCountMax(isotop), mfRangei.getIsotopeCountMax(isotop));
//	        	}
//	         	ima = MolecularFormulaManipulator.getString(MolecularFormulaRangeManipulator.getMaximalFormula(mfRangeInitialParent,builder));
//	            imi = MolecularFormulaManipulator.getString(MolecularFormulaRangeManipulator.getMinimalFormula(mfRangeInitialParent,builder));
//        	}else{
	        	ima = MolecularFormulaManipulator.getString(MolecularFormulaRangeManipulator.getMaximalFormula(mfRangei,builder));
	            imi = MolecularFormulaManipulator.getString(MolecularFormulaRangeManipulator.getMinimalFormula(mfRangei,builder));
//        	}
        	outPrint += "   "+imi+"\n";
        	outPrint += "   "+ima+"\n";
		}
		
		private IMolecularFormulaSet applyingLossConstraints(
				ParentIon parentIon, ParentIon fragmentIon) throws CDKException, IOException {
        	// If the fragmentIon doesn't contain any formula. The loss neither
			if(fragmentIon.getFormulaSet() == null){
        		   fragmentIon.addFormulaSet(null);
        		   return null;
        	}
        	
        	// Generating EC for neutral loss
            double differe = parentIon.getMass() - fragmentIon.getMass();
            IMolecularFormulaSet resultsMFL = generateMF(parentIon,fragmentIon,differe);
        	
            // Constraints to apply
        	if(resultsMFL != null)
			    resultsMFL = ruleNNeutral(resultsMFL);
            
            // test if someone loss MF is not fitted. Some loss + fragment => Parent
        	if(resultsMFL != null & (fragmentIon.getFormulaSet() != null))
            	resultsMFL = fittingRemovingLoss(parentIon, fragmentIon, resultsMFL);
            
            // test if someone fragment MF is not fitted. Some fragment + loss = Parent
            if(resultsMFL != null)
            	fragmentIon = fittingRemovingFragment(parentIon, fragmentIon, resultsMFL);
                
            // test if someone parent MF is not fitted. Some parent = fragment + loss
            if(resultsMFL != null)
              	parentIon = fittingRemovingParent(parentIon, fragmentIon, resultsMFL);
            
           // It it is not found any neutral loss the compound doesn't exist
       	   if(resultsMFL == null)
       		   fragmentIon.addFormulaSet(null);
       	   
       	   fragmentIon.addFormulaLossSet(resultsMFL);
       	   
       	   return resultsMFL;
		}
        
		private IMolecularFormula getMinimalMF(ParentIon compound) {
        	IMolecularFormula formulaMin = new MolecularFormula();
            
        	for(int f2 = 0 ; f2 < compound.getFragments().size(); f2++){
            	IMolecularFormula fragment2 = compound.getFragments().get(f2).getMinimalMF();
            	if(fragment2 == null)
            		continue;
            	Iterator<IIsotope> ii2 = fragment2.isotopes().iterator();
                while(ii2.hasNext()){
                    IIsotope isoto = ii2.next();
                    if(!formulaMin.contains(isoto)){
                    	formulaMin.addIsotope(isoto, fragment2.getIsotopeCount(isoto));
                    }else{
                    	int occurNew = fragment2.getIsotopeCount(isoto);
                    	int occurOld = formulaMin.getIsotopeCount(isoto);
                    	if(occurNew < occurOld){
                        	formulaMin.removeIsotope(isoto);
                    		formulaMin.addIsotope(isoto, occurNew);
                    	}
                    }
                }
            }
        	if(formulaMin.getIsotopeCount() == 0)
        		return null;
        	return formulaMin;
		}
		private IMolecularFormulaSet ruleNIon(IMolecularFormulaSet mfS) {
			if(!nitrogenR)
				return mfS;
				
			outPrint += "Not accepted due to  Rule_Nitrogen:\n";
			IRule rule  = new NitrogenRule();
        	IMolecularFormulaSet newmfS = new NNMolecularFormulaSet();

        	
    		try {
    			int count = 0;
	        	for(IMolecularFormula formula: mfS.molecularFormulas()){
//		        	if(polarity.equals(Polarity.positive))
//		        		formula.setCharge(1);
//		        	else if(polarity.equals(Polarity.negative))
//		        		formula.setCharge(-1);
//		        	else if(polarity.equals(Polarity.nothing))
//		        		formula.setCharge(0);
		        	
	        		String stringMF = MolecularFormulaManipulator.getString(formula);
//		    		IMolecularFormula newFormula = MolecularFormulaManipulator.getMajorIsotopeMolecularFormula(stringMF, builder);
		    		
					if(rule.validate(formula) == 1.0)
						newmfS.addMolecularFormula(formula);
					else{
						if(count == 0){
							outPrint += " "+stringMF+"\n";
						}else if(count == 7){
							count = 0;
							outPrint += "\n "+stringMF+"\n";
						}else{
							outPrint += ", "+stringMF+"\n";
						}
						count ++;
					}
					
	        	}
    		} catch (CDKException e) {
				e.printStackTrace();
			}
    		if(newmfS.size() == 0){
    			outPrint += "!! removed all MF !!\n";
    			newmfS = null;
    		}
    		outPrint += "---------------------------------------------\n";
            outPrint += "\n";
    		return newmfS;
		}
        
        private IMolecularFormulaSet ruleRDBE(IMolecularFormulaSet mfS) {
        	
        	/*should to inspected before*/
        	if(!RDBER)
				return mfS;
			
        	outPrint += "Not accepted due to  Rule_RDBE:\n";
        	
		    RDBERule rule  = new RDBERule();
        	IMolecularFormulaSet newmfS = new NNMolecularFormulaSet();
    		try {
    			int count = 0;
	        	for(IMolecularFormula formula: mfS.molecularFormulas()){
//		        	if(polarity.equals(Polarity.positive))
//		        		formula.setCharge(1);
//		        	else if(polarity.equals(Polarity.negative))
//		        		formula.setCharge(-1);
		        	
	        		String stringMF = MolecularFormulaManipulator.getString(formula);
//		    		IMolecularFormula newFormula = MolecularFormulaManipulator.getMajorIsotopeMolecularFormula(stringMF, builder);
//		    		System.out.print(stringMF);
	            	if(rule.validate(formula) == 1.0){
						newmfS.addMolecularFormula(formula);
					}else{
						if(count == 0){
							outPrint += " "+stringMF+"\n";
						}else if(count == 7){
							count = 0;
							outPrint += "\n "+stringMF+"\n";
						}else{
							outPrint += ", "+stringMF+"\n";
						}
						count ++;
					}
	        	}
			} catch (CDKException e) {
				e.printStackTrace();
			}

        	if(newmfS.size() == 0){
    			outPrint += "!! removed all MF !!\n";
    			newmfS = null;
    		}
        	outPrint += "---------------------------------------------\n";
            outPrint += "\n";
    		return newmfS;
		}

        private IMolecularFormulaSet ruleIsotopePattern(IMolecularFormulaSet mfS) {
        	if(!isotopePatternR)
				return mfS;
			
        	outPrint += "Not accepted due to  Rule_IsotopePattern:\n";
		   
			IMolecularFormulaSet newmfS = new NNMolecularFormulaSet();
			IsotopePatternSimilarity is = new IsotopePatternSimilarity();
			
			int count = 0;
        	for(IMolecularFormula formula: mfS.molecularFormulas()){
	        	
        		String stringMF = MolecularFormulaManipulator.getString(formula);
	    		IsotopePatternGenerator isotopeGe = new IsotopePatternGenerator(minAbISO);
    			outPrint += "\n";
	    		IsotopePattern patternIsoPredicted = isotopeGe.getIsotopes(formula);
	    		is.seTolerance(toleranceISO);
	    		
//	    		for(int i = 0 ; i < patternIsoPredicted.getNumberOfIsotopes(); i++)
//	        			System.out.println(stringMF+">"+patternIsoPredicted.getIsotope(i).getMass()+":"+patternIsoPredicted.getIsotope(i).getIntensity());
	        	    			
        		double tempScore = is.compare(isotoPattern, patternIsoPredicted);
    			outPrint += "formula: "+stringMF+"- tempScore: "+tempScore+"\n";
    	    	
            	if(tempScore > scoreISO){
					newmfS.addMolecularFormula(formula);
				}else{
					if(count == 0){
						outPrint += " Remved > "+stringMF+"\n";
					}else if(count == 7){
						count = 0;
						outPrint += "\n "+stringMF+"\n";
					}else{
						outPrint += ", "+stringMF+"\n";
					}
					count ++;
				}
        	}

        	if(newmfS.size() == 0){
    			outPrint += "!! removed all MF !!\n";
    			newmfS = null;
    		}
    		outPrint += "\n";
        	outPrint += "---------------------------------------------\n";
    		outPrint += "\n";
    		return newmfS;
		}
       
        private IMolecularFormulaSet ruleNNeutral(IMolecularFormulaSet mfS) {
			if(!nitrogenR)
				return mfS;

			outPrint += "Not accepted due to  Rule_Nitrogen(neutral loss):\n";
			IRule rule  = new NitrogenRule();
        	IMolecularFormulaSet newmfS = new NNMolecularFormulaSet();
    		try {
    			int count = 0;
    			if(mfS == null)
    				return mfS;
        	
    			for(IMolecularFormula formula: mfS.molecularFormulas()){
		        	
	        		String stringMF = MolecularFormulaManipulator.getString(formula);
		    		IMolecularFormula newFormula = MolecularFormulaManipulator.getMajorIsotopeMolecularFormula(stringMF, builder);
		    			if(rule.validate(formula) == 1.0)
							newmfS.addMolecularFormula(newFormula);
						else{
							if(count == 0){
								outPrint += " "+stringMF+"\n";
							}else if(count == 7){
								count = 0;
								outPrint += "\n "+stringMF+"\n";
							}else{
								outPrint += ", "+stringMF+"\n";
							}
							count ++;
						}
	        	}
			} catch (CDKException e) {
				e.printStackTrace();
			}
			if(newmfS.size() == 0){
    			outPrint += "!! removed all MF !!\n";
    			newmfS = null;
    		}
    		outPrint += "\n";
        	outPrint += "---------------------------------------------\n";
    		outPrint += "\n";
    		return newmfS;
		}
		/**
         * Print final result.
         *  
         * @param compound1  The parent ion
         * @param cycle
         */
        private void printResultFinal(ParentIon compound1,int cycle) {

        	int countFragLevel = 0;
        	List<String> stringID = new ArrayList<String>();
        	String printes = "";
//        	if(flagPrintFinal){
//                System.out.println();
//                printes = "################################################################";
//                System.out.println(printes);
//                printOut += printes+ "\n ";
//                double toleranPPM = compound1.getAccuracy()*1000000/compound1.getMass();
//                printes = "RESULT: (FINAL), aflter cycle: "+cycle;
//                System.out.println(printes);
//                printOut += printes+ "\n ";
//                printes = "      "+levels+" levelMax, "+nrElements+" elements, "+toleranPPM+" final ppm, ";
//                System.out.println(printes);
//                printOut += printes+ "\n ";
//                printes = "################################################################";
//                System.out.println(printes);
//                printOut += printes+ "\n ";
//                System.out.println();
//                printes = "";
//                System.out.println(printes);
//                printOut += printes+ "\n ";
//                printes = compound1.getID()+" > NODE [1] ("+compound1.getMass()+", +- "+compound1.getAccuracy()+")-----------------";
//                System.out.println(printes);
//                printOut += printes+ "\n ";
//                printes = "Nr MF: "+compound1.getFormulaSet().size()+", ";
//                System.out.print(printes);
//                printOut += printes;
//            }
        	
//        	if(flagPrintStati){
//        		System.out.print(compound1.getFormulaSet().size()+";");
//        		System.outPrint += cycle+";");
//        	}
        	
        	if(levels == 1)
            	countFragLevel +=  1;
        	
//            if(flagPrintFinal){
//                for(int j = 0; j < compound1.getFormulaSet().size(); j++){
//                	printes = MolecularFormulaManipulator.getString(compound1.getFormulaSet().getMolecularFormula(j));
//                    System.out.println(printes);
//                    printOut += printes+"\n";
//                }
//            }
            // level 2
            if(levels >= 2){
            	if(levels == 2){
            		boolean flagBreaking = false;
            		for(int i2 = 0; i2 < compound1.getFragments().size(); i2++){
                		if(!stringID.contains(compound1.getFragments().get(i2).getID())){
                			stringID.add(compound1.getFragments().get(i2).getID());
                		}else
                			flagBreaking = true;
            		}
            		if(flagBreaking)
            			return;
                	
                	countFragLevel +=  compound1.getFragments().size();
                }
                
                
                for(int i2 = 0; i2 < compound1.getFragments().size(); i2++){
                        
                    ParentIon compound22 = compound1.getFragments().get(i2);

//                    compound22.setAccuracy(getAccuracy(compound22,compound22.getFormulaSet()));
                    
//                	if(flagPrintFinal){
//                        printes = "";
//                        System.out.println(printes);
//                        printOut += printes+"\n";
//                        printes = compound22.getID()+" > NODE [2] ("+compound22.getMass()+", +- "+compound22.getAccuracy()+")-----------------";
//                        System.out.println(printes);
//                        printOut += printes+"\n";
//                        
//                        for(int j2 = 0; j2 < compound22.getFormulaSet().size(); j2++){
//                        	printes = MolecularFormulaManipulator.getString(compound22.getFormulaSet().getMolecularFormula(j2))+" -> ";
//                            if(compound22.getLoss() != null)
//                        		printes += MolecularFormulaManipulator.getString(compound22.getLoss());
//                            else
//                            	printes += "+"+searchCorrespondingLoss(compound1,compound22.getFormulaSet().getMolecularFormula(j2),compound22.getFormulaLossSet());
//                            printOut += printes+"\n";
//                        	System.out.println(printes);
//                            
//                    	}
//                	}
                	
                    // level 3
	                if(levels >= 3){
	                	if(levels == 3){
                        		boolean flagBreaking = false;
                        		for(int i3 = 0; i3 < compound22.getFragments().size(); i3++){
	                        		if(!stringID.contains(compound22.getFragments().get(i3).getID())){
	                        			stringID.add(compound22.getFragments().get(i3).getID());
	                        		}else
	                        			flagBreaking = true;
                        		}
                        		if(flagBreaking)
                        			continue;
                        	countFragLevel +=  compound22.getFragments().size();
                        }
                        
                        for(int i3 = 0; i3 < compound22.getFragments().size(); i3++){
                                
                            ParentIon compound33 = compound22.getFragments().get(i3);
//	                        compound33.setAccuracy(getAccuracy(compound33,compound33.getFormulaSet()));

//    	                	if(flagPrintFinal){
//    	                		 printes = "";
//    	                         System.out.println(printes);
//    	                         printOut += printes+"\n";
//    	                         printes = compound33.getID()+" > NODE [3] ("+compound33.getMass()+", +- "+compound33.getAccuracy()+")-----------------";
//    	                         System.out.println(printes);
//    	                         printOut += printes+"\n";
//		                        for(int j3 = 0; j3 < compound33.getFormulaSet().size(); j3++){
//		                        	printes = MolecularFormulaManipulator.getString(compound33.getFormulaSet().getMolecularFormula(j3))+" -> ";
//		                            if(compound33.getLoss() != null)
//		                        		printes += MolecularFormulaManipulator.getString(compound33.getLoss());
//		                            else
//		                            	printes += "+"+searchCorrespondingLoss(compound22,compound33.getFormulaSet().getMolecularFormula(j3),compound33.getFormulaLossSet());
//		                            printOut += printes+"\n";
//		                        	System.out.println(printes);
//	                        	}
//    	                	}
    	                	
	                        // level 4
        	                if(levels >= 4){
        	                	if(levels == 4){

	                        		boolean flagBreaking = false;
	                        		for(int i4 = 0; i4 < compound33.getFragments().size(); i4++){
		                        		if(!stringID.contains(compound33.getFragments().get(i4).getID())){
		                        			stringID.add(compound33.getFragments().get(i4).getID());
		                        		}else
		                        			flagBreaking = true;
	                        		}
	                        		if(flagBreaking)
	                        			continue;
		                        	countFragLevel +=  compound33.getFragments().size();
		                        }
		                        
		                        for(int i4 = 0; i4 < compound33.getFragments().size(); i4++){
		                                
		                            ParentIon compound44 = compound33.getFragments().get(i4);
		                            if(compound44.getFormulaSet() == null)
		                            	compound44.setAccuracy(0.0);
//		                            else
//		                            	compound44.setAccuracy(getAccuracy(compound44,compound44.getFormulaSet()));
		                            
//	        	                	if(flagPrintFinal){
//	        	                		 printes = "";
//	        	                         System.out.println(printes);
//	        	                         printOut += printes+"\n";
//	        	                         printes = compound44.getID()+" > NODE [4] ("+compound44.getMass()+", +- "+compound44.getAccuracy()+")-----------------";
//	        	                         System.out.println(printes);
//	        	                         printOut += printes+"\n";
//	        	                         if(compound44.getFormulaSet() != null)
//		        	                         for(int j4 = 0; j4 < compound44.getFormulaSet().size(); j4++){
//		        	                        	 printes = MolecularFormulaManipulator.getString(compound44.getFormulaSet().getMolecularFormula(j4))+" -> ";
//		        	                             if(compound44.getLoss() != null)
//		        	                         		printes += MolecularFormulaManipulator.getString(compound44.getLoss());
//		        	                             else
//		     		                            	printes += "+"+searchCorrespondingLoss(compound33,compound44.getFormulaSet().getMolecularFormula(j4),compound44.getFormulaLossSet());
//		        	                             printOut += printes+"\n";
//		        	                         	System.out.println(printes);
//				                        	}
//	        	                         
//	        	                	}
	        	                	
			                        // level 5
		        	                if(levels >= 5){
		        	                	if(levels == 5){

			                        		boolean flagBreaking = false;
			                        		for(int i5 = 0; i5 < compound44.getFragments().size(); i5++){
				                        		if(!stringID.contains(compound44.getFragments().get(i5).getID())){
				                        			stringID.add(compound44.getFragments().get(i5).getID());
				                        		}else
				                        			flagBreaking = true;
			                        		}
			                        		if(flagBreaking)
			                        			continue;
				                        	countFragLevel +=  compound44.getFragments().size();
				                        }
				                        
				                        for(int i5 = 0; i5 < compound44.getFragments().size(); i5++){
				                                
				                            ParentIon compound55 = compound44.getFragments().get(i5);
//					                        compound55.setAccuracy(getAccuracy(compound55,compound55.getFormulaSet()));
				                            

//			        	                	if(flagPrintFinal){
//			        	                		 printes = "";
//			        	                         System.out.println(printes);
//			        	                         printOut += printes+"\n";
//			        	                         printes = compound55.getID()+" > NODE [5] ("+compound55.getMass()+", +- "+compound55.getAccuracy()+")-----------------";
//			        	                         System.out.println(printes);
//			        	                         printOut += printes+"\n";
//			        	                         for(int j5 = 0; j5 < compound55.getFormulaSet().size(); j5++){
//			        	                        	 printes = MolecularFormulaManipulator.getString(compound55.getFormulaSet().getMolecularFormula(j5))+" -> ";
//			        	                             if(compound55.getLoss() != null)
//			        	                         		printes += MolecularFormulaManipulator.getString(compound55.getLoss());
//			        	                             else
//			     		                            	printes += "+"+searchCorrespondingLoss(compound44,compound55.getFormulaSet().getMolecularFormula(j5),compound55.getFormulaLossSet());
//			        	                             printOut += printes+"\n";
//			        	                         	System.out.println(printes);
//					                        	}
//			        	                	}
			        	                	
					                        // level 6
				        	                if(levels >= 6){
				        	                	if(levels == 6){

					                        		boolean flagBreaking = false;
					                        		for(int i6 = 0; i6 < compound55.getFragments().size(); i6++){
						                        		if(!stringID.contains(compound55.getFragments().get(i6).getID())){
						                        			stringID.add(compound55.getFragments().get(i6).getID());
						                        		}else
						                        			flagBreaking = true;
					                        		}
					                        		if(flagBreaking)
					                        			continue;
						                        	countFragLevel +=  compound55.getFragments().size();
						                        }
						                        
						                        for(int i6 = 0; i6 < compound55.getFragments().size(); i6++){
						                                
						                            ParentIon compound66 = compound55.getFragments().get(i6);
//							                        compound66.setAccuracy(getAccuracy(compound66,compound66.getFormulaSet()));
						                            

//					        	                	if(flagPrintFinal){
//					        	                		 printes = "";
//					        	                         System.out.println(printes);
//					        	                         printOut += printes+"\n";
//					        	                         printes = compound66.getID()+" > NODE [6] ("+compound66.getMass()+", +- "+compound66.getAccuracy()+")-----------------";
//					        	                         System.out.println(printes);
//					        	                         printOut += printes+"\n";
//					        	                         for(int j6 = 0; j6 < compound66.getFormulaSet().size(); j6++){
//					        	                        	 printes = MolecularFormulaManipulator.getString(compound66.getFormulaSet().getMolecularFormula(j6))+" -> ";
//					        	                             if(compound66.getLoss() != null)
//					        	                         		printes += MolecularFormulaManipulator.getString(compound66.getLoss());
//					        	                             else
//					     		                            	printes += "+"+searchCorrespondingLoss(compound55,compound66.getFormulaSet().getMolecularFormula(j6),compound66.getFormulaLossSet());
//					        	                             printOut += printes+"\n";
//					        	                         	System.out.println(printes);
//							                        	}
//					        	                	}
					        	                	
					        	                	// level 7
						        	                if(levels >= 7){
						        	                	if(levels == 7){

							                        		boolean flagBreaking = false;
							                        		for(int i7 = 0; i7 < compound66.getFragments().size(); i7++){
								                        		if(!stringID.contains(compound66.getFragments().get(i7).getID())){
								                        			stringID.add(compound66.getFragments().get(i7).getID());
								                        		}else
								                        			flagBreaking = true;
							                        		}
							                        		if(flagBreaking)
							                        			continue;
								                        	countFragLevel +=  compound66.getFragments().size();
								                        }
								                        
								                        for(int i7 = 0; i7 < compound66.getFragments().size(); i7++){
								                                
								                            ParentIon compound77 = compound66.getFragments().get(i7);
//									                        compound77.setAccuracy(getAccuracy(compound77,compound77.getFormulaSet()));
								                            

//							        	                	if(flagPrintFinal){
//							        	                		 printes = "";
//							        	                         System.out.println(printes);
//							        	                         printOut += printes+"\n";
//							        	                         printes = compound77.getID()+" > NODE [7] ("+compound77.getMass()+", +- "+compound77.getAccuracy()+")-----------------";
//							        	                         System.out.println(printes);
//							        	                         printOut += printes+"\n";
//							        	                         for(int j7 = 0; j7 < compound77.getFormulaSet().size(); j7++){
//							        	                        	 printes = MolecularFormulaManipulator.getString(compound77.getFormulaSet().getMolecularFormula(j7))+" -> ";
//							        	                             if(compound77.getLoss() != null)
//							        	                         		printes += MolecularFormulaManipulator.getString(compound77.getLoss());
//							        	                             else
//							     		                            	printes += "+"+searchCorrespondingLoss(compound66,compound77.getFormulaSet().getMolecularFormula(j7),compound77.getFormulaLossSet());
//							        	                             printOut += printes+"\n";
//							        	                         	System.out.println(printes);
//									                        	}
//							        	                	}
								                        }
						        	                }// level 7
					        	                	
						                        }
				        	                }// level 6
				                        }
		        	                }// level 5
		                        }
        	                }// level 4    
                        }
	                }// level 3
                }
            }// level 2

//        	if(flagPrintStati)
//        		System.out.println(countFragLevel);
        	
//        	if(countRedAccur == 0)
//        		this.numberCFrag = countFragLevel;
//        	break;
		
		}
        /**
         * search for the corresponding loss for a specific compound.
         * 
         * @param compound22
         * @param molecularFormula
         * @param formulaLossSet
         * 
         * @return The loss as String
         */
		private String searchCorrespondingLoss(ParentIon parent,
				IMolecularFormula fragment,
				IMolecularFormulaSet formulaLossSet) {
			for(int ls = 0; ls < formulaLossSet.size(); ls ++){ //fragment
				IMolecularFormula fSum = new MolecularFormula();
	        	fSum.add(fragment);
	        	fSum.add(formulaLossSet.getMolecularFormula(ls));
	    		
	    		// compare if the sum is some of the some compounds1
	    		boolean same = MolecularFormulaSetManipulator.contains(parent.getFormulaSet(), fSum);
	    		if(same)
	    			return MolecularFormulaManipulator.getString(formulaLossSet.getMolecularFormula(ls));
        	}
			return "";
		}
		/**
         * Check if the parent is correctly fitted. If not it is removed from the
         * set.
         * 
         * @param parent      The parent
         * @param fragment    The fragment
         * @param lossSet     The loss set
         * @return            The loss set without the not fitted loss  
         */
        private ParentIon fittingRemovingParent(
        			ParentIon parent,
        			ParentIon fragment, 
        			IMolecularFormulaSet lossSet) {

        	IMolecularFormulaSet fSetFrag = fragment.getFormulaSet();
        	IMolecularFormulaSet fSetPare = parent.getFormulaSet();
        	
        	outPrint += "Inspection PARENT:\n";
        	for(int ls = 0; ls < fSetPare.size(); ls ++){ // Parent
        		IMolecularFormula parFormula = fSetPare.getMolecularFormula(ls);
        		outPrint += "  for "+MolecularFormulaManipulator.getString(parFormula)+"\n";
            	boolean found2 = false;
        		for(IMolecularFormula fragFormula: fSetFrag.molecularFormulas()){ // fragment
            		for(IMolecularFormula lossFormula: lossSet.molecularFormulas()){ // loss
            			
            			IMolecularFormula fSum = new MolecularFormula();
            			fSum.setCharge(0);
            			fSum.add(fragFormula);
                		fSum.add(lossFormula);
                		boolean same = MolecularFormulaManipulator.compare(parFormula, fSum);

                    	//////////////////////////////////////////////////////////////////
                		outPrint += MolecularFormulaManipulator.getString(parFormula)+" = "+MolecularFormulaManipulator.getString(fragFormula)
                    				+" + "+MolecularFormulaManipulator.getString(lossFormula)+" -> "+same+"\n";
                		//////////////////////////////////////////////////////////////////
                		
                		if(same){ // if it is found go out of the look
                    		outPrint += "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
                        	found2 = true;
                        	break;
                		}
            		}
            		if(found2) // if it is found go out of the look
            			break;
            	}
            	if(!found2){ // if it is not found remove
        			//////////////////////////////////////////////////////////////////
        			outPrint += "Removing Parent: "+MolecularFormulaManipulator.getString(parFormula)+"\n";
        			//////////////////////////////////////////////////////////////////
        			fSetPare.removeMolecularFormula(parFormula);
        			parent.addFormulaSet(fSetPare);
        		}
             	outPrint += "\n";
        	}
         	outPrint += "\n";
// 			if(fSet22.size() == 0){
//				outPrint += "removed all MF\n");
//				fSet22 = null;
//			}

			return parent;
		}
        /**
         * Check if the parent is correctly fitted when you have one fragment and one loss. If not it is removed from the
         * set.
         * 
         * @param parent      The parent
         * @param fragment    The fragment
         * @param lossSet     The loss set
         * @return            The loss set without the not fitted loss  
         */
        private ParentIon fittingRemovingParent1(
        			ParentIon parent,
        			ParentIon fragment, 
        			IMolecularFormulaSet lossSet) {

        	IMolecularFormulaSet fSetFrag = fragment.getFormulaSet();
        	IMolecularFormulaSet fSetPare = parent.getFormulaSet();
//        	if(lossSet == null){
//				outPrint += "return the same lossSet is null\n";
//    			return parent;
//        	}
        	if(fSetFrag.size() != 1 || lossSet.size() != 1)
        		return parent;
        	
        	if(fSetPare.size() <= 1)
        		return parent;
        	

        	outPrint += "Inspection PARENT 1_TO_1:\n";
        	for(int ls = 0; ls < fSetPare.size(); ls ++){ // Parent
        		IMolecularFormula parFormula = fSetPare.getMolecularFormula(ls);
        		outPrint += "  for "+MolecularFormulaManipulator.getString(parFormula)+"\n";
            	boolean found2 = false;
        		for(IMolecularFormula fragFormula: fSetFrag.molecularFormulas()){ // fragment
            		for(IMolecularFormula lossFormula: lossSet.molecularFormulas()){ // loss
            			
            			IMolecularFormula fSum = new MolecularFormula();
            			fSum.setCharge(0);
            			fSum.add(fragFormula);
                		fSum.add(lossFormula);
                		boolean same = MolecularFormulaManipulator.compare(parFormula, fSum);

                    	//////////////////////////////////////////////////////////////////
                		outPrint += MolecularFormulaManipulator.getString(parFormula)+" = "+MolecularFormulaManipulator.getString(fragFormula)
                    				+" + "+MolecularFormulaManipulator.getString(lossFormula)+" -> "+same+"\n";
                		//////////////////////////////////////////////////////////////////
                		
                		if(same){ // if it is found go out of the look
                    		outPrint += "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
                        	found2 = true;
                        	break;
                		}
            		}
            		if(found2) // if it is found go out of the look
            			break;
            	}
            	if(!found2){ // if it is not found remove
        			//////////////////////////////////////////////////////////////////
        			outPrint += "Removing Parent: "+MolecularFormulaManipulator.getString(parFormula)+"\n";
        			//////////////////////////////////////////////////////////////////
        			fSetPare.removeMolecularFormula(parFormula);
        			parent.addFormulaSet(fSetPare);
        		}
             	outPrint += "\n";
        	}
         	outPrint += "\n";
// 			if(fSet22.size() == 0){
//				outPrint += "removed all MF\n");
//				fSet22 = null;
//			}

			return parent;
		}
		/**
         * Check if the Fragment is correctly fitted. If not it is removed from the
         * set.
         * 
         * @param parent      The parent
         * @param fragment    The fragment
         * @param lossSet     The loss set
         * @return            The loss set without the not fitted loss  
         */
        private ParentIon fittingRemovingFragment(
        			ParentIon parent,
        			ParentIon fragment, 
        			IMolecularFormulaSet lossSet) {
        	
        	IMolecularFormulaSet fSet22 = fragment.getFormulaSet();
        	outPrint += "Inspection FRAGMENT:\n";
        	for(int ls = 0; ls < fSet22.size(); ls ++){ // fragment
        		boolean flagFoundSum = false;
            	IMolecularFormula fVT = fSet22.getMolecularFormula(ls);
            	outPrint += "  for "+MolecularFormulaManipulator.getString(fVT)+"\n";
            	for(IMolecularFormula lossFormula: lossSet.molecularFormulas()){ // loss
            		IMolecularFormula fSum = new MolecularFormula();
                	fSum.setCharge(0);
            		fSum.add(fVT);
            		fSum.add(lossFormula);
            		
            		// compare if the sum is in someone combination of the parents
            		boolean same = MolecularFormulaSetManipulator.contains(parent.getFormulaSet(), fSum);

        			//////////////////////////////////////////////////////////////////
            			outPrint += MolecularFormulaManipulator.getString(fSum)+" = "+MolecularFormulaManipulator.getString(fVT)
                				+" + "+MolecularFormulaManipulator.getString(lossFormula)+" -> "+same+"\n";
            		//////////////////////////////////////////////////////////////////
            		
                    if(same){
                		outPrint += "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
                    	flagFoundSum = true;
                    	break;
                    }
            	}
            	outPrint += "\n";
            	if(!flagFoundSum){
            		//////////////////////////////////////////////////////////////////
            		outPrint += "Fragments to be removed:\n";
                	outPrint += "     "+MolecularFormulaManipulator.getString(fVT)+"\n";
            		//////////////////////////////////////////////////////////////////
                	fSet22.removeMolecularFormula(fVT);
                 	outPrint += "\n";
            	}
        	}
			if(fSet22.size() == 0){
				outPrint += "removed all MF\n";
				fSet22 = null;
			}
    		outPrint += "-------------------------------------\n";
        	outPrint += "\n";
			fragment.addFormulaSet(fSet22);
			return fragment;
		}

		/**
         * Check if the Loss is correctly fitted. If not it is removed from the
         * set.
         * 
         * @param parent      The parent
         * @param fragment    The fragment
         * @param lossSet     The loss set
         * @return            The loss set without the not fitted loss  
         */
        private IMolecularFormulaSet fittingRemovingLoss(
        			ParentIon parent,
        			ParentIon fragment, 
        			IMolecularFormulaSet lossSet) {
        	
        	outPrint += "Inspection LOSS: \n";
        	IMolecularFormulaSet formulasToRemove = new MolecularFormulaSet();
//        	if(lossSet == null){
//				outPrint += "removed all neutral MF because lossSet is null\n";
//    			return null;
//        	}
        	for(IMolecularFormula lossFormula: lossSet.molecularFormulas()){
        		
        		outPrint += "  for "+MolecularFormulaManipulator.getString(lossFormula)+"\n";
            	boolean flagSame = false;
        		IMolecularFormulaSet fragmentSet = fragment.getFormulaSet();
        		for(int ls = 0; ls < fragmentSet.size(); ls ++){ //fragment
            		IMolecularFormula fSum = new MolecularFormula();
            		fSum.setCharge(0);
                	IMolecularFormula fVT = fragmentSet.getMolecularFormula(ls);
                	fSum.add(fVT);
                	
                	fSum.add(lossFormula);
                	
                	// compare if the sum is some of the some compounds1
                	boolean same = MolecularFormulaSetManipulator.contains(parent.getFormulaSet(), fSum);
            		//////////////////////////////////////////////////////////////////
           			outPrint += "     "+MolecularFormulaManipulator.getString(fSum)+"["+fSum.getCharge()+"/"+parent.getFormulaSet().getMolecularFormula(0).getCharge()+"]"
                				+" ("+MolecularFormulaManipulator.getString(lossFormula)+") contained -> "+same+"\n";
                	//////////////////////////////////////////////////////////////////
                	
            		if(same){
            			flagSame = true;
                    	break;
                    }
            	}
            	if(!flagSame){
            		outPrint += "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
                	formulasToRemove.addMolecularFormula(lossFormula);
            	}
        	}
        	
        	outPrint += "\n";
        	outPrint += "Neutral loss being removed:\n";
        	for(IMolecularFormula lossFormula: formulasToRemove.molecularFormulas()){
        		//////////////////////////////////////////////////////////////////
        		outPrint += "     "+MolecularFormulaManipulator.getString(lossFormula)+"\n";
        		//////////////////////////////////////////////////////////////////
        		lossSet.removeMolecularFormula(lossFormula);
        	}
        	
        	if(lossSet.size() == 0){
				outPrint += "removed all neutral MF\n";
    			return null;
			}
        	outPrint += "\n";
    		outPrint += "-------------------------------------\n";
        	outPrint += "\n";
        	
        	return lossSet;
		}

        /**
         * Generate the possible EC for a fragment from a range of MolecularFormulaRange.
         * 
         * @param fragment    The Fragment ion to generate EC
         * @param mfRange     Range containing the EC constraint.
         * @param accuracy    Accuracy to take account
         * @return            A IMolecularFormulaSet containing the possible EC
         * 
         * @throws CDKException
         * @throws IOException 
         */
        private IMolecularFormulaSet generateMF(ParentIon fragment,
				MolecularFormulaRange mfRange,
				double accuracy) throws CDKException, IOException {

            ToleranceRangeRule ruleToleran = new ToleranceRangeRule();
            Object[] paramsT = new Object[2];
            paramsT[0] = 133.0;
            paramsT[1] = accuracy;
            ruleToleran.setParameters(paramsT);
            
            MassToFormulaTool mToF = new MassToFormulaTool(builder);
            List<IRule> rules = new ArrayList<IRule>();
            rules.add(ruleToleran);
            
            ElementRule elementRule  = new ElementRule();
            Object[] params = new Object[1];
            params[0] = mfRange;
            elementRule.setParameters(params);
            rules.add(elementRule);
            
            mToF.setRestrictions(rules);
            double mass = fragment.getMass();
            IMolecularFormulaSet resultsMF = mToF.generate(mass);
            if(resultsMF != null)
	            for(IMolecularFormula formula:resultsMF.molecularFormulas()){
	                if(polarity.equals(Polarity.positive))
	            		formula.setCharge(1);
	            	else if(polarity.equals(Polarity.negative))
	            		formula.setCharge(-1);
	            	else
	            		formula.setCharge(0);
	            }
            ////////////////// PRINTING ///////////////////////////////////////////
        	printingResultInitialCycle(resultsMF, mfRange);
     		////////////////// PRINTING ///////////////////////////////////////////
        	
            return resultsMF;
		}
        
        /**
         * Generate the possible EC for the loss total from his mass.
         * 
         * @param compoundPa    The parent ion. Given the EC constraints
         * @param compoundCh    The fragment ion
         * @param massLoss      Mass of the loss total
         * @return              A IMolecularFormulaSet containing the possible EC for the loss
         * @throws CDKException
         * @throws IOException
         */
		private IMolecularFormulaSet generateMF(ParentIon parent,
				ParentIon fragment, double massLoss) throws CDKException, IOException {
        	
        	List<IRule> rulesLoss = new ArrayList<IRule>();
             
            Object[] paramsTL = new Object[2];
            paramsTL[0] = 133.0;
            // checking which a
            double acc1 = fragment.getAccuracy();
            double acc2 = parent.getAccuracy();
            if(acc1 > acc2)
            	paramsTL[1] = acc1;
            else
            	paramsTL[1] = acc2;
            
            ToleranceRangeRule ruleToleranLoss = new ToleranceRangeRule();
            ruleToleranLoss.setParameters(paramsTL);
            
            MolecularFormulaRange mfRangeLoss = createRangeFromMax(parent.getMaximalMF());

            ElementRule elementRuleLoss  = new ElementRule();
            Object[] paramsLoss = new Object[1];
            paramsLoss[0] = mfRangeLoss;
            elementRuleLoss.setParameters(paramsLoss);
            
            rulesLoss.add(elementRuleLoss);
            rulesLoss.add(ruleToleranLoss);
            
            MassToFormulaTool mToFLoss = new MassToFormulaTool(builder); 
            mToFLoss.setRestrictions(rulesLoss);
            IMolecularFormulaSet resultsMFL2 = mToFLoss.generate(massLoss);
            
            if(resultsMFL2 != null)
	            for(IMolecularFormula formula:resultsMFL2.molecularFormulas()){
	            	formula.setCharge(0);
//	                if(polarity.equals(Polarity.positive))
//	            		formula.setCharge(1);
//	            	else if(polarity.equals(Polarity.negative))
//	            		formula.setCharge(-1);
	            }
            //////////////////PRINTING ///////////////////////////////////////////
        	printingResultInitialCycle(resultsMFL2, mfRangeLoss);
     		////////////////// PRINTING ///////////////////////////////////////////
			return resultsMFL2;
		}

		/**
         * Creates a minimal IMolecularFormula which contains for all IElements
         * with a occurrence of zero from other IMolecularFormula. 
         * Exception for Carbon which will be one occurrence.
         * 
         * @param mfMax IMolecularFormula to extract the IElements
         * @return      The minimal IMolecularFormula
         */
		private IMolecularFormula creatingMinimal0(IMolecularFormula mfMax) {
			IMolecularFormula newMin = new MolecularFormula();
//			double mass = MolecularFormulaManipulator.getMajorIsotopeMass(mfMax);
            Iterator<IIsotope> itIsotope2 = mfMax.isotopes().iterator();
            while(itIsotope2.hasNext()){
            	IIsotope isotope2 = itIsotope2.next();
            	if(isotope2.getSymbol().equals("C")){
//            		if(mass < 200.00){
            			newMin.addIsotope(isotope2, 0);
//            		}else{
//            			newMin.addIsotope(isotope2, 1);
//                	}
            	}else
            		newMin.addIsotope(isotope2, 0);
            }
			return newMin;
		}
		private void printingResultTitle(ParentIon compound){
			outPrint += "\n";
			outPrint += "+++++++++++++++++++++++++++++++++++++++++++++\n";
			outPrint += "\n";
            outPrint += "NODE ID = "+compound.getID()+" \n";
			outPrint += "\n";
            outPrint += "Mass: "+compound.getMass()+"\n";
            outPrint += "Acur:"+compound.getAccuracy()+"\n";
			outPrint += "\n";
		}
		private void printingResultTitleLoss(ParentIon parent, ParentIon fragment){
			outPrint += "\n";
			outPrint += "+++++++++++++++++++++++++++++++++++++++++++++\n";
			outPrint += "\n";
            outPrint += "LOSS ID = "+fragment.getID()+"loss \n";
			outPrint += "\n";
			double diff = parent.getMass() - fragment.getMass();
            outPrint += "Mass: "+diff+"\n";
            outPrint += "Acur: "+fragment.getAccuracy()+"\n";
			outPrint += "\n";
		}
				

		private void printingResultInitialCycle(IMolecularFormulaSet resultsMF, MolecularFormulaRange mfRangei) throws IOException {
		    
            String ima = MolecularFormulaManipulator.getString(MolecularFormulaRangeManipulator.getMaximalFormula(mfRangei,builder));
            String imi = MolecularFormulaManipulator.getString(MolecularFormulaRangeManipulator.getMinimalFormula(mfRangei,builder));
            
            outPrint += "MF Max. Range = "+ima+"\n";
            outPrint += "MF Min. Range = "+imi+"\n";
            outPrint += "---------------------------------------------\n";
            outPrint += "\n";
			if(resultsMF == null){
            	outPrint += "     Nr. MF: 0\n";
            	outPrint += "\n";
            	return;
            }else
            	outPrint += "     Nr. MF: "+resultsMF.size()+"\n";
			outPrint += "\n";
			
            int count = 0;
            for(int i = 0; i < resultsMF.size(); i++){
            	String stringMF = MolecularFormulaManipulator.getString(resultsMF.getMolecularFormula(i));
    			if(count == 0){
    				outPrint += " "+stringMF+"\n";
    			}else if(count == 7){
    				count = 0;
    				outPrint += "\n "+stringMF+"\n";
    			}else{
    				outPrint += ", "+stringMF+"\n";
    			}
    			count ++;
    		}outPrint += "\n";
            	
            MolecularFormulaRange mfRange1 = MolecularFormulaRangeManipulator.getRange(resultsMF);
	        IMolecularFormula maximalEle2T = MolecularFormulaRangeManipulator.getMaximalFormula(mfRange1,builder);
	        IMolecularFormula minimalEle11 = MolecularFormulaRangeManipulator.getMinimalFormula(mfRange1,builder);
	        
            outPrint += "MF Max. Range Post = "+MolecularFormulaManipulator.getString(maximalEle2T)+"\n";
            outPrint += "MF Min. Range Post = "+MolecularFormulaManipulator.getString(minimalEle11)+"\n";
            outPrint += "---------------------------------------------\n";
            
            outPrint += "\n";
		}
		/**
		 * Creating an ElementRule for the Elemental composition of the Loss total
		 * 
		 * @param mfMax   The IMolecularFormula to take as representation
		 * @return        The ElementRule object
		 * @throws IOException
		 * @throws CDKException
		 */
		private MolecularFormulaRange createRangeFromMax(IMolecularFormula mfMax) throws IOException, CDKException {
            
            IsotopeFactory ifac = IsotopeFactory.getInstance(DefaultChemObjectBuilder.getInstance());
            MolecularFormulaRange mfRangeLoss = new MolecularFormulaRange();
            Iterator<IIsotope> it = mfMax.isotopes().iterator();
            while(it.hasNext()){
            	IIsotope isotope = it.next();
            	mfRangeLoss.addIsotope(ifac.getMajorIsotope(isotope.getSymbol()),0,mfMax.getIsotopeCount(isotope));
            }
            
			return mfRangeLoss;
		}
		/**
		 * Get a new accuracy extracted from the parent and a group of IMolecularFormulaSet.
		 * 
		 * @param parent      The Parent ion
		 * @param formulaSet  The IMolecularFormulaSet
		 * @return            The new accuracy
		 */
		private double getAccuracy(ParentIon parent, IMolecularFormulaSet formulaSet) {
			double newAccuracy = 0;
			for (IMolecularFormula formula : formulaSet.molecularFormulas()){
				double diff = Math.abs(MolecularFormulaManipulator.getTotalExactMass(formula) - parent.getMass());
				if(diff > newAccuracy)
					newAccuracy = diff;
			}
			return newAccuracy+0.001;
		}

		/**
		 * The controller checks if the cycles is finished. It is based on 
		 * looking if in all nodes there is one only possibility combination
		 * 
		 * @param compound1  The initial parent ion 
		 * @return           True, if the cycle is finish
		 */
		private boolean isMonoResult(ParentIon compound1) {
			
			if(compound1.getFormulaSet().size() != 1)
				return false;

			if(levels >= 2)
            for(int i2 = 0; i2 < compound1.getFragments().size(); i2++){
            	ParentIon compound2 = compound1.getFragments().get(i2);
            	if(compound2.getFormulaSet() != null && compound2.getFormulaSet().size() != 1)
            		if(compound2.getFormulaLossSet().size() != 1)
					return false;
    			if(levels >= 3)
            	for(int i3 = 0; i3 < compound2.getFragments().size(); i3++){
                	ParentIon compound3 = compound2.getFragments().get(i3);
                	if(compound3.getFormulaSet() != null && compound3.getFormulaSet().size() != 1)
    					return false;
        			if(levels >= 4)
                	for(int i4 = 0; i4 < compound3.getFragments().size(); i4++){
                    	ParentIon compound4 = compound3.getFragments().get(i4);
                    	if(compound4.getFormulaSet() != null && compound4.getFormulaSet().size() != 1)
        					return false;
            			if(levels >= 5)
                    	for(int i5 = 0; i5 < compound4.getFragments().size(); i5++){
                        	ParentIon compound5 = compound4.getFragments().get(i5);
                        	if(compound5.getFormulaSet() != null && compound5.getFormulaSet().size() != 1)
            					return false;
                			if(levels >= 6)
                        	for(int i6 = 0; i6 < compound5.getFragments().size(); i6++){
                            	ParentIon compound6 = compound5.getFragments().get(i6);
                            	if(compound6.getFormulaSet() != null && compound6.getFormulaSet().size() != 1)
                					return false;
                    			if(levels >= 7)
                            	for(int i7 = 0; i7 < compound6.getFragments().size(); i7++){
                                	ParentIon compound7 = compound6.getFragments().get(i7);
                                	if(compound7.getFormulaSet() != null && compound7.getFormulaSet().size() != 1)
                    					return false;
                                    
                                }
                            }
                        }
                    }
                }
            }
            return true;
		}
		/**
		 * The controller checks if the cycles is finished. It is based on 
		 * looking if in all nodes there is one only possibility combination
		 * 
		 * @param compound1  The initial parent ion 
		 * @return           True, if the cycle is finish
		 */
		private boolean isFinishCycle1(ParentIon compound1) {
			
			
            if(levels >= 2){// levels 2
                for(int i2 = 0; i2 < compound1.getFragments().size(); i2++){
                    ParentIon compound22 = compound1.getFragments().get(i2);
                    if(compound22.getFormulaSet() != null){
                    		if(compound22.getFormulaSet().size() != 1)
                    			return false;
                    		if(compound22.getFormulaLossSet().size() != 1)
                    			return false;
                    }
	                if(levels >= 3){// levels 3
                        for(int i3 = 0; i3 < compound22.getFragments().size(); i3++){
                                
                            ParentIon compound33 = compound22.getFragments().get(i3);
                            if(compound33.getFormulaSet() != null){
                        		if(compound33.getFormulaSet().size() != 1)
                        			return false;
                        		if(compound33.getFormulaLossSet().size() != 1)
                        			return false;
                            }
                            
        	                if(levels >= 4){// levels 4
		                        for(int i4 = 0; i4 < compound33.getFragments().size(); i4++){
		                                
		                            ParentIon compound44 = compound33.getFragments().get(i4);
		                            if(compound44.getFormulaSet() != null){
		                        		if(compound44.getFormulaSet() == null)
		                        			continue;
		                        		if(compound44.getFormulaLossSet().size() != 1)
		                        			return false;
		                            }
		                            
		        	                if(levels >= 5){// levels 5
				                        for(int i5 = 0; i5 < compound44.getFragments().size(); i5++){
				                                
				                            ParentIon compound55 = compound44.getFragments().get(i5);
				                            if(compound55.getFormulaSet() != null ){
				                        		if(compound55.getFormulaSet().size() != 1)
				                        			return false;

				                        		if(compound55.getFormulaLossSet().size() != 1)
				                        			return false;
				                            }
				        	                if(levels >= 6){// levels 6
						                        for(int i6 = 0; i6 < compound55.getFragments().size(); i6++){
						                                
						                            ParentIon compound66 = compound55.getFragments().get(i6);
						                            if(compound66.getFormulaSet() != null){
						                        		if(compound66.getFormulaSet().size() != 1)
						                        			return false;

						                        		if(compound66.getFormulaLossSet().size() != 1)
						                        			return false;
						                            }
						                            
						        	                if(levels >= 7){// levels 7
								                        for(int i7 = 0; i7 < compound66.getFragments().size(); i7++){
								                                
								                            ParentIon compound77 = compound66.getFragments().get(i7);
								                            if(compound77.getFormulaSet() != null){
								                        		if(compound77.getFormulaSet().size() != 1)
								                        			return false;

								                        		if(compound77.getFormulaLossSet().size() != 1)
								                        			return false;
								                            }
								                            
								                        }
						        	                }// levels 7
						                        }
				        	                }// levels 6
				                        }
		        	                }// levels 5
		                        }
        	                }// levels 4    
                        }
	                }// levels 3
                }
            }// levels 2
            
            return true;
		}
		

		/**
		 * The controller checks if the cycles is finished. It is based on 
		 * looking if the occurrence of each node remains equal
		 * 
		 * @param compound1  The initial parent ion 
		 * @return           True, if the cycle is finish
		 */
		private boolean isFinishCycle2(ParentIon compound1) {
			boolean flag = true;
			
            if(levels >= 2){// levels 2
                for(int i2 = 0; i2 < compound1.getFragments().size(); i2++){
                    
                    ParentIon compound22 = compound1.getFragments().get(i2);
                    if(compound22.getFormulaSet() == null)
                    	continue;
                    int numNew2 = compound22.getFormulaSet().size();
                    Integer numOld2I = mapOccur.get(compound22.getID());
                    if(numOld2I == null){
                    	mapOccur.put(compound22.getID(), numNew2);
                    	flag = false;
                    }else{
                    
	                    int numOld2 = numOld2I.intValue();
	                    if(numNew2 != numOld2){
	                    	mapOccur.put(compound22.getID(), numNew2);
	                    	flag = false;
	                    }
                    }
                    // levels 3
	                if(levels >= 3){
                        for(int i3 = 0; i3 < compound22.getFragments().size(); i3++){
                                
                            ParentIon compound33 = compound22.getFragments().get(i3);
                            if(compound33.getFormulaSet() == null)
                            	continue;
                            int numNew3 = compound33.getFormulaSet().size();
                            Integer numOld3I = mapOccur.get(compound33.getID());
                            if(numOld3I == null){
                            	mapOccur.put(compound33.getID(), numNew3);
                            	flag = false;
                            }else{
                            
	                            int numOld3 = numOld3I.intValue();
	                            if(numNew3 != numOld3){
	                            	mapOccur.put(compound33.getID(), numNew3);
	                            	flag = false;
	                            }
                            }
	                        // levels 4
        	                if(levels >= 4){
		                        for(int i4 = 0; i4 < compound33.getFragments().size(); i4++){
		                                
		                            ParentIon compound44 = compound33.getFragments().get(i4);
		                            
		                            if(compound44.getFormulaSet() == null)
		                            	continue;
		                            int numNew4 = compound44.getFormulaSet().size();
		                            Integer numOld4I = mapOccur.get(compound44.getID());
		                            if(numOld4I == null){
		                            	mapOccur.put(compound44.getID(), numNew4);
		                            	flag = false;
		                            }else{
		                            
			                            int numOld4 = numOld4I.intValue();
			                            if(numNew4 != numOld4){
			                            	mapOccur.put(compound44.getID(), numNew4);
			                            	flag = false;
			                            }
		                            }
			                        // levels 5
		        	                if(levels >= 5){
				                        for(int i5 = 0; i5 < compound44.getFragments().size(); i5++){
				                                
				                            ParentIon compound55 = compound44.getFragments().get(i5);
				                            if(compound55.getFormulaSet() == null)
				                            	continue;
				                            int numNew5 = compound55.getFormulaSet().size();
				                            Integer numOld5I = mapOccur.get(compound55.getID());
				                            if(numOld5I == null){
				                            	mapOccur.put(compound55.getID(), numNew5);
				                            	flag = false;
				                            }else{
				                            
					                            int numOld5 = numOld5I.intValue();
					                            if(numNew5 != numOld5){
					                            	mapOccur.put(compound55.getID(), numNew5);
					                            	flag = false;
					                            }
				                            }
					                        // levels 6
				        	                if(levels >= 6){
						                        for(int i6 = 0; i6 < compound55.getFragments().size(); i6++){
						                                
						                            ParentIon compound66 = compound55.getFragments().get(i6);
						                            if(compound66.getFormulaSet() == null)
						                            	continue;
						                            int numNew6 = compound66.getFormulaSet().size();
						                            Integer numOld6I = mapOccur.get(compound66.getID());
						                            if(numOld6I == null){
						                            	mapOccur.put(compound66.getID(), numNew6);
						                            	flag = false;
						                            }else{
							                            
							                            int numOld6 = numOld6I.intValue();
							                            if(numNew6 != numOld6){
							                            	mapOccur.put(compound66.getID(), numNew6);
							                            	flag = false;
							                            }
						                            }
						                            // levels 7
						        	                if(levels >= 7){
								                        for(int i7 = 0; i7 < compound66.getFragments().size(); i7++){
								                                
								                            ParentIon compound77 = compound66.getFragments().get(i7);
								                            if(compound77.getFormulaSet() == null)
								                            	continue;
								                            int numNew7 = compound77.getFormulaSet().size();
								                            Integer numOld7I = mapOccur.get(compound77.getID());
								                            if(numOld7I == null){
								                            	mapOccur.put(compound77.getID(), numNew7);
								                            	flag = false;
								                            }else{
									                            
									                            int numOld7 = numOld7I.intValue();
									                            if(numNew7 != numOld7){
									                            	mapOccur.put(compound77.getID(), numNew7);
									                            	flag = false;
									                            }
								                            }
								                        }
						        	                }// levels 7
						                        }
				        	                }// levels 6
				                        }
		        	                }// levels 5
		                        }
        	                }// levels 4    
                        }
	                }// levels 3
                }
            }// levelss 2
            
            if(!flag){
            	this.countInt = 0;
            	return false;
            }else if(this.countInt == 3)
            	return true;
            else{
            	this.countInt++;
            	return false;
            }
		}
		public String getOutprint() {
			return outPrint;
		}
}