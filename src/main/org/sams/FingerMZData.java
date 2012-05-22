package org.sams;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

import org.sams.manipulator.MZDataManipulator;

public class FingerMZData {

	private int size = 1024;
	/**
	 * Signatures of lineal level 0
	 * 
	 * @param mzData
	 * @param listfingers
	 * @return
	 */
	public BitSet getFingerprint(MZData mzData,List<String> fingerprintList) {
		List<String> fingInt = MZDataManipulator.getListPath(mzData);
//		BitSet bitSet = new BitSet(fingerprintList.size());
		int count = 0;
		BitSet bitSet = new BitSet(size);


//        Set<String> cleanPath = new HashSet<String>();
//        for (StringBuffer s : allPaths) {
//            if (cleanPath.contains( s.toString() )) continue;
//            String s2 = s.reverse().toString();
//            if (cleanPath.contains(s2)) continue;
//            cleanPath.add(s2);
//        }
//
//        // convert paths to hashes
//        int[] hashes = new int[fingInt.size()];
//        int i= 0;
//        for (String s: fingInt) hashes[i++] = s.hashCode();
//
//        for (int hash : hashes) {
//            int position = new java.util.Random(hash).nextInt(size);
//            bitSet.set(position);
//        }for (String finger : fingInt) {
//            if(fingerprintList.contains(finger))
//            	bitSet.set(count);
//            count ++;
//        }
        

//		BitSet bitSet = new BitSet(size);
//
//        // convert paths to hashes
//        int[] hashes = new int[fingInt.size()];
//        int i= 0;
//        for (String s: fingInt) hashes[i++] = s.hashCode();
//
//        for (int hash : hashes) {
//            int position = new java.util.Random(hash).nextInt(size);
//            bitSet.set(position);
//        }
        
		return bitSet;
	}
	/**
	 * Signatures of fragments lineal level 0. What we know as EFP.
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintL0(MZData mzData) {
		return MZDataManipulator.getListPath(mzData);
	}

	/**
	 * Signatures of losses lineal level 0. What we know as EFP.
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintL0_loss(MZData mzData) {
		return MZDataManipulator.getListPathLoss(mzData);
	}
	/**
	 * Signatures of fragments lineal level 1. It will contains only EC.
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintL1(MZData mzData) {
		return getFingerprintL(mzData,1);
	}
	/**
	 * Signatures of Losses lineal level 1 It will contains only EC.
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintL1_Loss(MZData mzData) {
		return getFingerprintL_Loss(mzData,1);
	}
	/**
	 * Signatures of fragments lineal level 2. Connected two EC.
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintL2(MZData mzData) {
		return getFingerprintL(mzData,2);
	}
	/**
	 * Signatures of losses lineal level 2. Connected two EC.
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintL2_Loss(MZData mzData) {
		return getFingerprintL_Loss(mzData,2);
	}
	/**
	 * Signatures of lineal fragments level 3. Connected two EC.
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintL3(MZData mzData) {
		return getFingerprintL(mzData,3);
	}
	/**
	 * Signatures of lineal losses level 3. Connected two EC.
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintL3_Loss(MZData mzData) {
		return getFingerprintL_Loss(mzData,3);
	}
	/**
	 * Signatures of lineal fragments level 4. Connected two EC.
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintL4(MZData mzData) {
		return getFingerprintL(mzData,4);
	}
	/**
	 * Signatures of lineal losses level 4. Connected two EC.
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintL4_Loss(MZData mzData) {
		return getFingerprintL_Loss(mzData,4);
	}
	/**
	 * Signatures of lineal fragments level 5. Connected two EC.
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintL5(MZData mzData) {
		return getFingerprintL(mzData,5);
	}
	/**
	 * Signatures of lineal losses level 5. Connected two EC.
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintL5_Loss(MZData mzData) {
		return getFingerprintL_Loss(mzData,5);
	}
	
	////////////////////////
	/**
	 * Signatures of fragments lineal level 0. What we know as EFP.
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintNL0(MZData mzData) {
		return MZDataManipulator.getListPathNominal(mzData);
	}

	/**
	 * Signatures of losses lineal level 0. What we know as EFP.
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintNL0_loss(MZData mzData) {
		return MZDataManipulator.getListPathNominal_loss(mzData);
	}

	/**
	 * Signatures of fragments lineal level 1. It will contains only Nominal.
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintNomL1(MZData mzData) {
		return getFingerprintNomL(mzData,1);
	}
	/**
	 * Signatures of Losses lineal level 1 It will contains only Nominal.
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintNomL1_Loss(MZData mzData) {
		return getFingerprintNomL_Loss(mzData,1);
	}
	/**
	 * Signatures of fragments lineal level 2. Connected two Nominal.
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintNomL2(MZData mzData) {
		return getFingerprintNomL(mzData,2);
	}
	/**
	 * Signatures of losses lineal level 2. Connected two Nominal.
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintNomL2_Loss(MZData mzData) {
		return getFingerprintNomL_Loss(mzData,2);
	}
	/**
	 * Signatures of lineal fragments level 3. Connected two Nominal.
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintNomL3(MZData mzData) {
		return getFingerprintNomL(mzData,3);
	}
	/**
	 * Signatures of lineal losses level 3. Connected two Nominal.
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintNomL3_Loss(MZData mzData) {
		return getFingerprintNomL_Loss(mzData,3);
	}
	/**
	 * Signatures of lineal fragments level 4. Connected two Nominal.
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintNomL4(MZData mzData) {
		return getFingerprintNomL(mzData,4);
	}
	/**
	 * Signatures of lineal losses level 4. Connected two Nominal.
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintNomL4_Loss(MZData mzData) {
		return getFingerprintNomL_Loss(mzData,4);
	}
	/**
	 * Signatures of lineal fragments level 5. Connected two Nominal.
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintNomL5(MZData mzData) {
		return getFingerprintNomL(mzData,5);
	}
	/**
	 * Signatures of lineal losses level 5. Connected two Nominal.
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintNomL5_Loss(MZData mzData) {
		return getFingerprintNomL_Loss(mzData,5);
	}
	///////////////////////

	private List<String> getFingerprintL(MZData mzData, int sign) {
		List<String> listPath = MZDataManipulator.getListPath(mzData);
		return getFingerprintL(listPath,sign);
	}

	private List<String> getFingerprintNomL(MZData mzData, int sign) {
		List<String> listPath = MZDataManipulator.getListPathNominal(mzData);
		return getFingerprintL(listPath,sign);
	}
	
	private List<String> getFingerprintL(List<String> listPath, int sign) {
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
		return removeDuplicate(fingerList);
	}

	private List<String> getFingerprintL_Loss(MZData mzData, int sign) {
		List<String> listPath = MZDataManipulator.getListPathLoss(mzData);
		return getFingerprintL_Loss(listPath,sign);
	}

	private List<String> getFingerprintNomL_Loss(MZData mzData, int sign) {
		List<String> listPath = MZDataManipulator.getListPathNominal_loss(mzData);
		return getFingerprintL_Loss(listPath,sign);
	}
	
	private List<String> getFingerprintL_Loss(List<String> listPath, int sign) {
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
		return removeDuplicate(fingerList);
	}

	/**
	 * Signatures of fragments parallel level 0
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintV0(MZData mzData) {
		return MZDataManipulator.getListPathV(mzData);
	}

	/**
	 * Signatures of losses parallel level 0
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintV0_Loss(MZData mzData) {
		return MZDataManipulator.getListPathVLoss(mzData);
	}
	/**
	 * Signatures of fragments Vertical level 2
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintV2(MZData mzData) {
		return getFingerprintV(mzData,2);
	}
	/**
	 * Signatures of losses Vertical level 2
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintV2_Loss(MZData mzData) {
		return getFingerprintV_Loss(mzData,2);
	}
	/**
	 * Signatures of fragments Vertical level 3
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintV3(MZData mzData) {
		return getFingerprintV(mzData,3);
	}
	/**
	 * Signatures of losses Vertical level 3
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintV3_Loss(MZData mzData) {
		return getFingerprintV_Loss(mzData,3);
	}
	
	////////////////////////////////////////
	/**
	 * Signatures of fragments Vertical level 2
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintNomV2(MZData mzData) {
		return getFingerprintNomV(mzData,2);
	}
	/**
	 * Signatures of losses Vertical level 2
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintNomV2_Loss(MZData mzData) {
		return getFingerprintNomV_Loss(mzData,2);
	}
	/**
	 * Signatures of fragments Vertical level 3
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintNomV3(MZData mzData) {
		return getFingerprintNomV(mzData,3);
	}
	/**
	 * Signatures of losses Vertical level 3
	 * 
	 * @param listPath
	 * @return
	 */
	public List<String> getFingerprintNomV3_Loss(MZData mzData) {
		return getFingerprintNomV_Loss(mzData,3);
	}

	////////////////////////////////////////

	private List<String> getFingerprintNomV(MZData mzData, int sign) {
		List<String> listPath = MZDataManipulator.getListPathNominal(mzData);
		return getFingerprintV(listPath,sign);
	}

	private List<String> getFingerprintV(MZData mzData, int sign) {
		List<String> listPath = MZDataManipulator.getListPath(mzData);
		return getFingerprintV(listPath,sign);
	}
	
	private List<String> getFingerprintV(List<String> listPath, int sign) {
		List<String> fingerList = new ArrayList<String>();
		for(String pathParent:listPath){
			List<String> listNei = new ArrayList<String>();
			String[] fm = pathParent.split("\\|\\|");
			
			for(String pathKid:listPath){

				String parent = "";
				if(!pathKid.contains("||"))
					continue;
				else 
					parent = pathKid.substring(0, pathKid.lastIndexOf("||"));
				
				if(parent.equals(pathParent)){
					String[] fmK = pathKid.split("\\|\\|");
					listNei.add(fmK[fmK.length-1]);
				}
					
			}
			if(listNei.size() < 2)
				continue;
			List<String> listPackNei = getPack(listNei,sign);
			for(String pathPack:listPackNei){
				fingerList.add(fm[fm.length-1]+"@"+pathPack);
			}
			
		}
		return removeDuplicate(fingerList);
	}

	private List<String> getFingerprintV_Loss(MZData mzData, int sign) {
		List<String> listPath = MZDataManipulator.getListPathLoss(mzData);
		return getFingerprintV_Loss(listPath,sign);
	}

	private List<String> getFingerprintNomV_Loss(MZData mzData, int sign) {
		List<String> listPath = MZDataManipulator.getListPathNominal_loss(mzData);
		return getFingerprintV_Loss(listPath,sign);
	}

	private List<String> getFingerprintV_Loss(List<String> listPath, int sign) {
		List<String> fingerList = new ArrayList<String>();
		for(String pathParent:listPath){
			List<String> listNei = new ArrayList<String>();
			String[] fm = pathParent.split("\\|\\|");
			
			for(String pathKid:listPath){

				String parent = "";
				if(!pathKid.contains("||"))
					continue;
				else 
					parent = pathKid.substring(0, pathKid.lastIndexOf("||"));
				
				if(parent.equals(pathParent)){
					String[] fmK = pathKid.split("\\|\\|");
					listNei.add(fmK[fmK.length-1]);
				}
					
			}
			if(listNei.size() < 2)
				continue;
			List<String> listPackNei = getPack(listNei,sign);
			for(String pathPack:listPackNei){
				fingerList.add(fm[fm.length-1]+"@"+pathPack);
			}
			
		}
		return removeDuplicate(fingerList);
	}

	private List<String> removeDuplicate(List<String> fingerList) {
		List<String> newfingerList = new ArrayList<String>();
		HashSet<String> set = new HashSet<String>();
		 for (int i = 0; i < fingerList.size(); i++) {
		  boolean val = set.add(fingerList.get(i));
		  if (val == true) {
		  	newfingerList.add(fingerList.get(i));
		  }
		 }
		return newfingerList;
	}
	public List<String> getPackList(String pack, int nNei) {
    	String parent = pack.split("@")[0];
    	String rest = pack.split("@")[1];
    	rest = rest.substring(1, rest.length()-1);
    	String[] rests = rest.split("\\|\\|");
    	List<String> listFP = new ArrayList<String>(); 
    	List<String> listNei = Arrays.asList(rests);
    	List<String> listPackNei = getPack(listNei,nNei);
		for(String pathPack:listPackNei){
			listFP.add(parent+"@"+pathPack);
		}
		return listFP;
	}
	private List<String> getPack(List<String> listNei, int nNei) {
    	Collections.sort(listNei);
    	List<String> listPackNei =  new ArrayList<String>();
		for(int i = 0 ; i < listNei.size(); i++){
			String sG = "["+listNei.get(i);
			for(int j = i+1 ; j < listNei.size(); j++){
				if(j - i >= nNei-1){
					listPackNei.add(sG+"||"+listNei.get(j)+"]");
				}else{
					sG += "||"+listNei.get(j);
				}
				
        	}
		}
		return listPackNei;
	}
}
