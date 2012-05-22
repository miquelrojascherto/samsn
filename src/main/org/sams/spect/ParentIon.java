package org.sams.spect;

import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.formula.MolecularFormulaSet;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.interfaces.IMolecularFormulaSet;
import org.openscience.cdk.interfaces.IMolecule;
import org.xmlcml.cml.element.CMLPeak;
/**
 * Class which define a parent ion in a multi-stage mass spectrometry.
 * An example use would be a set of templates ...
 * 
 * @author Miguel Rojas-Cherto
 *
 */
public class ParentIon 	{
	/** Mass of the ion*/
	private double mass;
	/** Identification of the ion*/
	private String id;
	/** Level which is found the ion onto the spectral tree*/
	private int level;
	/** List of fragment ions which descents of this ion*/
	private List<ParentIon> fragments = new ArrayList<ParentIon>();
	/** Set of IMolecularFormula defining this ion*/
	private IMolecularFormulaSet formulaSet = new MolecularFormulaSet();
	private IMolecularFormulaSet formulaLossSet = new MolecularFormulaSet();
	private IMolecularFormula maximalMF;
	private IMolecularFormula minimalMF;
	private double accuracy = 0.0;
	private List<Double> massP;
	private List<Double> inte;
	private int position;
	private Double int_Abs;
	private Double int_Nor;
	private Double int_Nor2;
	private Double int_SN;
	private ParentIon parentIon;
	private String id_XCMS;
	private String rt;
	private String formulaPath;
	private IMolecule molecule;
	private Integer numGroup;
	private CMLPeak cmlPeak;
	private IMolecule moleculeLoss;
	/**
	 * constructor ParentIon Object.
	 * 
	 * @param mass  The mass of this ion
	 * @param intensity 
	 * @param id    The ion identification
	 * @param level The level of this ion in the spectral tree
	 */
	public ParentIon(double mass, double intensity, String id, int level){
		this.mass = mass;
		this.int_Abs = intensity;
		this.id = id;
		this.level = level;
	}
	/**
	 * Add the experimental Isotope pattern of this compound.
	 * 
	 * @param isotopes The isotope pattern
	 */
	public void addIsotopePattern(List<Double> mass,List<Double> inte){
		this.massP = mass;
		this.inte = inte;
	}

	/**
	 * Get the experimental Isotope pattern of this ion.
	 * 
	 * @return The isotope pattern
	 */
	public List<Double> getIsotopePatternMass(){
		return this.massP;
	}
	/**
	 * Get the experimental Isotope pattern of this ion.
	 * 
	 * @return The isotope pattern
	 */
	public List<Double> getIsotopePatternInte(){
		return this.inte;
	}
	/**
	 * Return the Level onto the spectral tree of this ion.
	 * @return
	 */
	public int getLevel(){
		return level;
	}
	/**
	 * Add a new fragment ion.
	 * 
	 * @param fragment
	 */
	public void addFragment(ParentIon fragment){
		fragments.add(fragment);
	}
	/**
	 * Add possible IMolecularFormulaSet for this ion.
	 * 
	 * @param formulaSet The IMolecularFormulaSet
	 */
	public void addFormulaSet(IMolecularFormulaSet formulaSet){
		this.formulaSet  = formulaSet;
	}

	/**
	 * Add possible IMolecularFormulaSet for this compound
	 * 
	 * @param formulaLossSet The IMolecularFormulaSet
	 */
	public void addFormulaLossSet(IMolecularFormulaSet formulaLossSet){
		this.formulaLossSet  = formulaLossSet;
	}
	
	/**
	 * Add maximal 
	 * @param maximalMF
	 */
	public void addMaximalFormula(IMolecularFormula maximalMF) {
		this.maximalMF = maximalMF;
		
	}
	

	/**
	 * Add minimal 
	 * @param minimalMF
	 */
	public void addMinimalFormula(IMolecularFormula minimalMF) {
		this.minimalMF = minimalMF;
		
	}

	/**
	 * Get the experimental mass of this Compound.
	 * 
	 * @return The experimental of this Compound.
	 */
	public double getMass(){
		return mass;
	}
	/**
	 * Set the mass of the ion.
	 * 
	 * @param mass The mass of the ion
	 */
	public void setMass(double mass){
		this.mass = mass;
	}
	
	/**
	 * Return the fragments.
	 * 
	 * @return A List with the fragments
	 */
	public List<ParentIon> getFragments(){
		return fragments;
	}
	/**
	 * Return IMolecularFormulaSet.
	 * 
	 * @return The IMolecularFormulaSet object
	 */
	public IMolecularFormulaSet getFormulaSet(){
		return formulaSet;
	}
	/**
	 * Return IMolecularFormulaSet.
	 * 
	 * @return The IMolecularFormulaSet object
	 */
	public IMolecularFormulaSet getFormulaLossSet(){
		return formulaLossSet;
	}
	/**
	 * Return Maximal MF.
	 * 
	 * @return The IMolecularFormula object
	 */
	public IMolecularFormula getMaximalMF(){
		return maximalMF;
	}

	/**
	 * Return minimal MF.
	 * 
	 * @return The IMolecularFormula object
	 */
	public IMolecularFormula getMinimalMF(){
		return minimalMF;
	}
	/**
	 * Get the Identification of this compound.
	 * 
	 * @return The identification value
	 */
	public String getID(){
		return id;
	}
	/**
	 * Set the Accuracy of this compound.
	 * 
	 * @return The identification value
	 */
	public double getAccuracy(){
		return accuracy;
	}
	/**
	 * Set the Accuracy of this compound.
	 * 
	 * @return The identification value
	 */
	public void setAccuracy(double accuracy){
		this.accuracy = accuracy;
	}
	/**
	 * Position in the acquisition
	 * 
	 * @param i position
	 */
	public void setPosition(int i) {
		this.position = i;
		
	}
	public int getPosition(){
		return position;
	}
	public void setIntensity_Abs(Double int_Abs) {
		this.int_Abs = int_Abs;
	}
	public Double getIntensity_Abs(){
		return int_Abs;
	}
	public void setIntensity_Nor(Double int_Nor) {
		this.int_Nor = int_Nor;
	}
	public Double getIntensity_Nor(){
		return int_Nor;
	}
	public void setIntensity_Nor2(Double int_Nor2) {
		this.int_Nor2 = int_Nor2;
	}
	public Double getIntensity_Nor2(){
		return int_Nor2;
	}
	public void setIntensity_SN(Double int_SN) {
		this.int_SN = int_SN;
	}
	public Double getIntensity_SN(){
		return int_SN;
	}
	public void addParent(ParentIon parentIon) {
		this.parentIon = parentIon;
	}
	public ParentIon getParent(){
		return parentIon;
	}
	public void addID_XCMS(String id_XCMS) {
		this.id_XCMS = id_XCMS;
	}
	public String getID_XCMS(){
		return id_XCMS;
	}
	public void setScan(String rt) {
		this.rt = rt;
	}
	public String getScan(){
		return rt;
	}
	public void addFormulaPath(String formulaPath) {
		this.formulaPath = formulaPath;
	}
	public String getFormulaPath() {
		return formulaPath;
	}
	public void setMolecule(IMolecule molecule) {
		this.molecule = molecule;
	}
	public IMolecule getMolecule() {
		return molecule;
	}
	public void setGroup(Integer numGroup) {
		this.numGroup = numGroup;
	}
	public Integer getGroup() {
		return numGroup;
	}
	public void addCMLPeak(CMLPeak cmlPeak) {
		this.cmlPeak = cmlPeak;
	}

	public CMLPeak getCMLPeak() {
		return cmlPeak;
	}
	public void setMoleculeLoss(IMolecule moleculeLoss) {
		this.moleculeLoss = moleculeLoss;
	}

	public IMolecule getMoleculeLoss() {
		return this.moleculeLoss;
	}
}
