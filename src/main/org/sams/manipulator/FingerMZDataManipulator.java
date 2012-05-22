package org.sams.manipulator;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.sams.FingerMZData;
import org.sams.MZData;

public class FingerMZDataManipulator {
	static DefaultChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
	
	/**
	 * Compare if the Value1 is contained in value2
	 * @param values1 The values to check
	 * @param values2 The values container
	 * @return        Percentatge of similarity
	 */
	public static double compare(List<String> values1, List<String> values2){
		
		double count = 0;
		for(String value1 : values1){
			if(values2.contains(value1)){
				count++;
			}
		}
		double value = count / values2.size();
		
		return value;
	}
	/**
	 * Extract the tanimoto coefficient from two list of strings.
	 * 
	 * @param values1 A list containing strings
	 * @param values2 A list containing strings
	 * @return        Tanimoto coefficient
	 */
	public static double tanimotto(List<String> values1, List<String> values2){
		
		double countS = 0.0;
		for(String value1 : values1){
			if(values2.contains(value1)){
				countS++;
			}
		}
		double value = countS / (values2.size() + values1.size() - countS);
		if(Double.isNaN(value))
			return 0.0;
		return value;
	}
	/**
	 * Extract the tanimotto coefficient from two list of strings. But only use those fingerprint that 
	 * are not repetitive.
	 * 
	 * @param values1 A list containing strings
	 * @param values2 A list containing strings
	 * @return        Tanimotto coefficient
	 */
	public static double tanimotto2(List<String> values1, List<String> values2){
		
		double countS = 0.0;
		
		List<String> values1R = new ArrayList<String>();
		List<String> values1S = new ArrayList<String>();
		int count1R = 0;
		int count1S = 0;
		for(String value1 : values1){
//			System.out.println(value1+" "+values2.contains(value1));
			for(String v1R : values1R){
				if(!value1.contains("[")){
					if(value1.startsWith(v1R)){
//						System.out.println("predec: "+v1R);
						count1R++;
						break;
					}
				}else{
//					System.out.println("value1: "+value1);
					if(value1.split("@").length > 2)
						continue;
					String p1 = value1.split("@")[0];
					String p2 = value1.split("@")[1].replace("]", "").replace("[", "");
					String p21 = p2.split("\\|\\|")[0];
					String p22 = p2.split("\\|\\|")[1];
					String p1F = p1+"@"+p21+"@";
					String p2F = p1+"@"+p22+"@";
//					System.out.println(value1+" > "+p1+" , "+p21+" "+p22);
					if(p1F.startsWith(v1R) || p2F.startsWith(v1R)){
//						System.out.println("predec: "+v1R);
						count1R++;
						break;
					}
//					break;
				}
			}
			if(values2.contains(value1)){
				countS++;
				values1S.add(value1+"@");
				for(String v1S : values1S){
					if(!value1.contains("[")){
						if(value1.startsWith(v1S)){
//							System.out.println(v1S+" st_ "+value1);
//							countS--;
							values1S.remove(v1S);
							count1S++;
							break;
						}
					}else{
//						System.out.println("value1: "+value1);
						if(value1.split("@").length > 2)
							continue;
						String p1 = value1.split("@")[0];
						String p2 = value1.split("@")[1].replace("]", "").replace("[", "");
						String p21 = p2.split("\\|\\|")[0];
						String p22 = p2.split("\\|\\|")[1];
						String p1F = p1+"@"+p21+"@";
						String p2F = p1+"@"+p22+"@";
//						System.out.println(value1+" > "+p1+" , "+p21+" "+p22);
						if(p1F.startsWith(v1S) || p2F.startsWith(v1S)){
//							System.out.println("predec: "+v1R);
							values1S.remove(v1S);
							count1S++;
							break;
						}
//						break;
					}
				}
						
			}else{
//				System.out.println("add "+value1+"@");
				values1R.add(value1+"@");
			}
		}

		List<String> values2R = new ArrayList<String>();
		List<String> values2S = new ArrayList<String>();
		int count2R = 0;
		int count2S = 0;
		for(String value2 : values2){
//			System.out.println(value1+" "+values2.contains(value1));
			for(String v2R : values2R){
				if(!value2.contains("[")){
					if(value2.startsWith(v2R)){
	//					System.out.println("predec: "+v2R);
						count2R++;
						break;
					}
				}else{
//					System.out.println("value1: "+value2);
					if(value2.split("@").length > 2 )
						continue;
					String p1 = value2.split("@")[0];
					String p2 = value2.split("@")[1].replace("]", "").replace("[", "");
					String p21 = p2.split("\\|\\|")[0];
					String p22 = p2.split("\\|\\|")[1];
					String p1F = p1+"@"+p21+"@";
					String p2F = p1+"@"+p22+"@";
//					System.out.println(value2+" > "+p1+" , "+p21+" "+p22);
					if(p1F.startsWith(v2R) || p2F.startsWith(v2R)){
//						System.out.println("predec: "+v2R);
						count2R++;
						break;
					}
//					break;
				}
			}
			if(values1.contains(value2)){
//				countS++;
				values2S.add(value2+"@");
				for(String v2S : values2S){
					if(!value2.contains("[")){
						if(value2.startsWith(v2S)){
//							System.out.println(v2S+" st "+value2);
//							countS--;
							values2S.remove(v2S);
//							values2.remove(value2);
							count2S++;
							break;
						}
					}
				}
			}else{
				values2R.add(value2+"@");
			}
		}
		double value = (countS-count1S) / ((values1.size() -count1R - count1S) + (values2.size() - count2R - count1S) - (countS - count1S));
//		System.out.println(value+" = ("+countS+" - "+count1S+") / (("+values1.size()+" - "+count1R+" - "+count1S+") + ("+values2.size()+" - "+count2R+" - "+count2S+") - ("+countS+" - "+count1S+"))");
		if(Double.isNaN(value))
			return 0.0;
		return value;
	}

	/**
	 * Extract the tanimotto coefficient from two list of strings. But only use those fingerprint that 
	 * are not repetitive.
	 * 
	 * @param values1 A list containing strings
	 * @param values2 A list containing strings
	 * @return        Tanimotto coefficient value
	 */
	public static double tanimotto3(List<String> values1, List<String> values2){
		
		double countC = 0.0;// number in common features C
		int countCP = 0; // number of overlapped features C'
		
		List<String> values1R = new ArrayList<String>();
		List<String> values1S = new ArrayList<String>();
		List<String> values1Sb = new ArrayList<String>();
		int count1R = 0; // number of B'
//		System.out.println(values1.size() +" == "+values2.size());

		for(String value1 : values1){
			if(!values2.contains(value1))
				values1R.add(value1+"@");
		}

		for(String value1 : values1){
			for(String v1R : values1R){
				if(!value1.contains("[")){
//					System.out.println(value1+" starts "+v1R+" "+value1.startsWith(v1R));
					if(value1.startsWith(v1R)){
//						System.out.println("predec: "+v1R);
						count1R++;
						break;
					}
				}else{
//					System.out.println("value1: "+value1);
					if(value1.split("@").length > 2)
						continue;
					String p1 = value1.split("@")[0];
					String p2 = value1.split("@")[1].replace("]", "").replace("[", "");
					String p21 = p2.split("\\|\\|")[0];
					String p22 = p2.split("\\|\\|")[1];
					String p1F = p1+"@"+p21+"@";
					String p2F = p1+"@"+p22+"@";
//					System.out.println(value1+" > "+p1+" , "+p21+" "+p22);
					if(p1F.startsWith(v1R) || p2F.startsWith(v1R)){
//						System.out.println("predec[]: "+v1R);
						count1R++;
						break;
					}
				}
			}

//			System.out.println(value1+" "+values2.contains(value1));
			/////////////////////////////////
			if(values2.contains(value1)){
				countC++; 
				values1S.add(value1+"@");
//				System.out.println(countC+", val1: "+value1+"@");
				for(String v1S : values1S){
					if(!value1.contains("[")){
//						System.out.println(v1S+" aa ");
						if(value1.startsWith(v1S)){
//							System.out.println(v1S+" st_ "+value1);
							if(!values1Sb.contains(v1S)){
								values1Sb.add(v1S);
								countCP++;
							}
						}else if((value1+"@").endsWith(v1S)){
							int nv = (value1+"@").split("@").length;
							int nvs = v1S.split("@").length;
							if(nv != nvs){
//								System.out.println(v1S+" end_ "+value1);
								if(!values1Sb.contains(v1S)){
									values1Sb.add(v1S);
									countCP++;
								}
							}
							
						}
					}else{
//						System.out.println("value1: "+value1);
						if(value1.split("@").length > 2)
							continue;
						String p1 = value1.split("@")[0];
						String p2 = value1.split("@")[1].replace("]", "").replace("[", "");
						String p21 = p2.split("\\|\\|")[0];
						String p22 = p2.split("\\|\\|")[1];
						String p1F = p1+"@"+p21+"@";
						String p2F = p1+"@"+p22+"@";
//						System.out.println(value1+" > "+p1+" , "+p21+" "+p22);
						if(p1F.startsWith(v1S) || p2F.startsWith(v1S)){
//							System.out.println("predec: "+v1R);
							if(!values1Sb.contains(v1S)){
								values1Sb.add(v1S);
								countCP++;
							}
							break;
						}
					}
				}	
			}
		}

		List<String> values2R = new ArrayList<String>();
		List<String> values2S = new ArrayList<String>();
		int count2R = 0; // number of A'
		for(String value2 : values2){
			if(!values1.contains(value2))
				values2R.add(value2+"@");
		}
		for(String value2 : values2){
//			System.out.println(value1+" "+values2.contains(value1));
			for(String v2R : values2R){
				if(!value2.contains("[")){
					if(value2.startsWith(v2R)){
	//					System.out.println("predec: "+v2R);
						count2R++;
						break;
					}
				}else{
//					System.out.println("value1: "+value2);
					if(value2.split("@").length > 2 )
						continue;
					String p1 = value2.split("@")[0];
					String p2 = value2.split("@")[1].replace("]", "").replace("[", "");
					String p21 = p2.split("\\|\\|")[0];
					String p22 = p2.split("\\|\\|")[1];
					String p1F = p1+"@"+p21+"@";
					String p2F = p1+"@"+p22+"@";
//					System.out.println(value2+" > "+p1+" , "+p21+" "+p22);
					if(p1F.startsWith(v2R) || p2F.startsWith(v2R)){
//						System.out.println("predec: "+v2R);
						count2R++;
						break;
					}
				}
			}
			if(values1.contains(value2)){
				values2S.add(value2+"@");
				for(String v2S : values2S){
					if(!value2.contains("[")){
						if(value2.startsWith(v2S)){
//							System.out.println(v2S+" st "+value2);
							values2S.remove(v2S);
							break;
						}
					}
				}
			}
		}
		double value = (countC-countCP) / ((values1.size() - count1R - countCP) + (values2.size() - count2R - countCP) - (countC - countCP));
//		System.out.println(value+" = ("+countC+" - "+countCP+") / (("+values1.size()+" - "+count1R+" - "+countCP+") + ("+values2.size()+" - "+count2R+" - "+ countCP+") - ("+countC+" - "+countCP+"))");
		if(Double.isNaN(value))
			return 0.0;
		return value;
	}

	/**
	 * Convert the Formula Path to Nominal Path identity
	 * 
	 * @param formulapathL A list containing the formula path
	 * @return A list containing the nominal path
	 */
	public static List<String> getNominalPath(List<String> bs151L) {
		return getNominalPath(bs151L, "#0.0");
	}
	/**
	 * Convert the Formula Path to Nominal Path identity
	 * 
	 * @param bs151L A list containing the formula path
	 * @return A list containing the nominal path
	 */
	public static List<String> getNominalPath(List<String> bs151L, String format) {
		List<String> formulapath = new ArrayList<String>();
		for(Object ob:bs151L){
			String path = (String)ob;
			if(path.length() < 1 || path.equals("||"))
				continue;
			// if only lineal AB||BC||DE
			if(!path.contains("@[")){
//				String sub[] = path.split("\\|\\|");
				String sub[] = path.split("@");

				String nomPath = "";
				for(int i = 0 ; i < sub.length; i++){
					String fm = sub[i];
					if(fm.length() < 1)
						continue;
					IMolecularFormula formula = MolecularFormulaManipulator.getMajorIsotopeMolecularFormula(fm,builder);
			    	formula.setCharge(+1);
			    	double mass = MolecularFormulaManipulator.getTotalExactMass(formula);
			    	DecimalFormat Currency = new DecimalFormat(format);
			        String formated = Currency.format(mass);
//			    	nomPath += formated+"||";
			    	nomPath += formated+"@";
				}
//				nomPath = nomPath.substring(0, nomPath.length()-2);
				nomPath = nomPath.substring(0, nomPath.length()-1);
				formulapath.add(nomPath);
			}else{
				// if only parallel AB@[BC||DE]
				DecimalFormat Currency = new DecimalFormat(format);
		        if(path.startsWith("@"))
		        	continue;
				String sub[] = path.split("\\@\\[");
				IMolecularFormula formula = MolecularFormulaManipulator.getMajorIsotopeMolecularFormula(sub[0],builder);
		    	formula.setCharge(+1);
		    	double mass = MolecularFormulaManipulator.getTotalExactMass(formula);
		    	String nmP = Currency.format(mass);
				String subSub[] = sub[1].replace("]", "").split("\\|\\|");
				String nomPath = "";
				formula = MolecularFormulaManipulator.getMajorIsotopeMolecularFormula(subSub[0],builder);
		    	formula.setCharge(+1);
		    	mass = MolecularFormulaManipulator.getTotalExactMass(formula);
		    	String nomValue = Currency.format(mass);
				String sG = "["+nomValue;
				for(int j = 1 ; j < subSub.length; j++){
					formula = MolecularFormulaManipulator.getMajorIsotopeMolecularFormula(subSub[j],builder);
			    	formula.setCharge(+1);
			    	mass = MolecularFormulaManipulator.getTotalExactMass(formula);
			    	nomValue = Currency.format(mass);
					if(j == subSub.length-1){
						sG = sG+"||"+nomValue+"]";
					}else{
						sG += "||"+nomValue;
					}
	        	}
				nomPath = nmP+"@"+sG;
				
				formulapath.add(nomPath);
			}
		}
		return formulapath;
	}

	/**
	 * Get the bit string for bs15
	 * 
	 * @param inchiList1_L The List in the MF
	 * @return A list with the bits
	 */
	public static List<String> getbs24_24_1Nom(MZData mzdata){

		FingerMZData fpMZD = new FingerMZData();
		
		List<String> bs11_L = fpMZD.getFingerprintNomL2(mzdata);
		List<String> bs21_L = fpMZD.getFingerprintNomL2_Loss(mzdata);
    	List<String> bs31_L = new ArrayList<String>();
		for(String b:bs11_L)
			bs31_L.add(b);
		for(String b:bs21_L)
			bs31_L.add(b);
		
		List<String> bs41_L = fpMZD.getFingerprintNomL3(mzdata);
		List<String> bs51_L = fpMZD.getFingerprintNomL3_Loss(mzdata);
		List<String> bs61_L = new ArrayList<String>();
		for(String b:bs41_L)
			bs61_L.add(b);
		for(String b:bs51_L)
			bs61_L.add(b);
		
		List<String> bs71_L = fpMZD.getFingerprintNomL4(mzdata);
		List<String> bs81_L = fpMZD.getFingerprintNomL4_Loss(mzdata);
		List<String> bs91_L = new ArrayList<String>();
		for(String b:bs71_L)
			bs91_L.add(b);
		for(String b:bs81_L)
			bs91_L.add(b);
		
		List<String> bs71_V = fpMZD.getFingerprintNomV2(mzdata);
		List<String> bs81_V = fpMZD.getFingerprintNomV2_Loss(mzdata);
		List<String> bs91_V = new ArrayList<String>();
		for(String b:bs71_V)
			bs91_V.add(b);
		for(String b:bs81_V)
			bs91_V.add(b);
		
//		List<String> bs11_V = fpMZD.getFingerprintV3(mzdata);
//		List<String> bs21_V = fpMZD.getFingerprintV3(mzdata);
//		List<String> bs31_V = new ArrayList<String>();
//		for(String b:bs11_V)
//			bs31_V.add(b);
//		for(String b:bs21_V)
//			bs31_V.add(b);
		
		List<String> bs151_L = new ArrayList<String>();
		for(String b:bs31_L)
			bs151_L.add(b);
		for(String b:bs61_L)
			bs151_L.add(b);
		for(String b:bs91_L)
			bs151_L.add(b);
//		for(String b:bs31_V)
//			bs151_L.add(b);
		for(String b:bs91_V)
			bs151_L.add(b);
		
		return FingerMZDataManipulator.removeDuplicates(bs151_L);
	}
	/**
	 * Get the bit string for bs15
	 * 
	 * @param inchiList1_L The List in the MF
	 * @return A list with the bits
	 */
	public static List<String> getbs24_24_1(MZData mzdata){

		FingerMZData fpMZD = new FingerMZData();
		
		List<String> bs11_L = fpMZD.getFingerprintL2(mzdata);
		List<String> bs21_L = fpMZD.getFingerprintL2_Loss(mzdata);
    	List<String> bs31_L = new ArrayList<String>();
		for(String b:bs11_L)
			bs31_L.add(b);
		for(String b:bs21_L)
			bs31_L.add(b);
		
		List<String> bs41_L = fpMZD.getFingerprintL3(mzdata);
		List<String> bs51_L = fpMZD.getFingerprintL3_Loss(mzdata);
		List<String> bs61_L = new ArrayList<String>();
		for(String b:bs41_L)
			bs61_L.add(b);
		for(String b:bs51_L)
			bs61_L.add(b);
		
		List<String> bs71_L = fpMZD.getFingerprintL4(mzdata);
		List<String> bs81_L = fpMZD.getFingerprintL4_Loss(mzdata);
		List<String> bs91_L = new ArrayList<String>();
		for(String b:bs71_L)
			bs91_L.add(b);
		for(String b:bs81_L)
			bs91_L.add(b);
		
		List<String> bs71_V = fpMZD.getFingerprintV2(mzdata);
		List<String> bs81_V = fpMZD.getFingerprintV2_Loss(mzdata);
		List<String> bs91_V = new ArrayList<String>();
		for(String b:bs71_V)
			bs91_V.add(b);
		for(String b:bs81_V)
			bs91_V.add(b);
		
//		List<String> bs11_V = fpMZD.getFingerprintV3(mzdata);
//		List<String> bs21_V = fpMZD.getFingerprintV3(mzdata);
//		List<String> bs31_V = new ArrayList<String>();
//		for(String b:bs11_V)
//			bs31_V.add(b);
//		for(String b:bs21_V)
//			bs31_V.add(b);
		
		List<String> bs151_L = new ArrayList<String>();
		for(String b:bs31_L)
			bs151_L.add(b);
		for(String b:bs61_L)
			bs151_L.add(b);
		for(String b:bs91_L)
			bs151_L.add(b);
//		for(String b:bs31_V)
//			bs151_L.add(b);
		for(String b:bs91_V)
			bs151_L.add(b);
		
		return bs151_L;
	}
	/**
	 * Get the bit string for bs15
	 * 
	 * @param inchiList1_L The List in the MF
	 * @return A list with the bits
	 */
	public static List<String> getbs24(List<String> inchiList1_L,List<String> inchiList1L_L){
		List<String> bs11_L = getFingerprintL(inchiList1_L, 2);
		List<String> bs21_L = getFingerprintL(inchiList1L_L, 2);
    	List<String> bs31_L = new ArrayList<String>();
		for(String b:bs11_L)
			bs31_L.add(b);
		for(String b:bs21_L)
			bs31_L.add(b);
		
		List<String> bs41_L = getFingerprintL(inchiList1_L, 3);
		List<String> bs51_L = getFingerprintL(inchiList1L_L, 3);
		List<String> bs61_L = new ArrayList<String>();
		for(String b:bs41_L)
			bs61_L.add(b);
		for(String b:bs51_L)
			bs61_L.add(b);
		
		List<String> bs71_L = getFingerprintL(inchiList1_L, 4);
		List<String> bs81_L = getFingerprintL(inchiList1L_L, 4);
		List<String> bs91_L = new ArrayList<String>();
		for(String b:bs71_L)
			bs91_L.add(b);
		for(String b:bs81_L)
			bs91_L.add(b);
		
		List<String> bs151_L = new ArrayList<String>();
		for(String b:bs31_L)
			bs151_L.add(b);
		for(String b:bs61_L)
			bs151_L.add(b);
		for(String b:bs91_L)
			bs151_L.add(b);
		
		return bs151_L;
	}

	/**
	 * Get the bit string for bs25
	 * 
	 * @param inchiList1_L The List in the MF
	 * @return A list with the bits
	 */
	public static List<String> getbs25(List<String> inchiList1_L,List<String> inchiList1L_L){
		List<String> bs11_L = getFingerprintL(inchiList1_L, 1);
		List<String> bs21_L = getFingerprintL(inchiList1L_L, 1);
    	List<String> bs31_L = new ArrayList<String>();
		for(String b:bs11_L)
			bs31_L.add(b);
		for(String b:bs21_L)
			bs31_L.add(b);
		
		List<String> bs41_L = getFingerprintL(inchiList1_L, 2);
		List<String> bs51_L = getFingerprintL(inchiList1L_L, 2);
		List<String> bs61_L = new ArrayList<String>();
		for(String b:bs41_L)
			bs61_L.add(b);
		for(String b:bs51_L)
			bs61_L.add(b);
		
		List<String> bs71_L = getFingerprintL(inchiList1_L, 3);
		List<String> bs81_L = getFingerprintL(inchiList1L_L, 3);
		List<String> bs91_L = new ArrayList<String>();
		for(String b:bs71_L)
			bs91_L.add(b);
		for(String b:bs81_L)
			bs91_L.add(b);
		
		List<String> bs151_L = new ArrayList<String>();
		for(String b:bs31_L)
			bs151_L.add(b);
		for(String b:bs61_L)
			bs151_L.add(b);
		for(String b:bs91_L)
			bs151_L.add(b);
		
		return bs151_L;
	}
	public static List<String> getFingerprintL(List<String> listPath, int sign) {
		List<String> fingerList = new ArrayList<String>();
		for(String path:listPath){
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
		}
		return fingerList;
	}
	public static List<String> removeDuplicates(List<String> listI){
		List<String> listR = new ArrayList<String>();
		for(String fp:listI){
			if(!listR.contains(fp)){
				listR.add(fp);
			}
		}
		return listR;
	}
}
