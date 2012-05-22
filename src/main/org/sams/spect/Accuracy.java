package org.sams.spect;

/**
 * Class which define the accuracy value.
 * 
 * @author Miguel Rojas-Cherto
 *
 */
public class Accuracy 	{
	/** The accuracy of the ion in Da*/
	private double accuDa;
	/** The accuracy of the ion in PPM*/
	private double accuPPM;

	/**
	 * Constructor of the Accuracy object.
	 */
	public Accuracy(){
	}
	/**
	 * Set the accuracy.
	 * 
	 * @param massIon  The mass of the ion
	 * @param accuracyDa The accuracy in Da
	 */
	public void setAccuracy(double massIon, double accuracyDa){
		accuDa = accuracyDa;
		accuPPM = accuracyDa*1000000/massIon;
	}
	/**
	 * Set the accuracy.
	 * 
	 * @param massIon     The mass of the ion
	 * @param accuracyPPM The accuracy in ppm
	 */
	public void setAccuracyPPM(double massIon, double accuracyPPM){
		accuPPM = accuracyPPM;
		accuDa = massIon*accuracyPPM/1000000;
	}
	/**
	 * get the accuracy in PPM.
	 * 
	 * @return The accuracy in PPM
	 */
	public double getAccuracyPPM(){
		return accuPPM;
	}
	/**
	 * get the accuracy in Da.
	 * 
	 * @return The accuracy in Da
	 */
	public double getAccuracyDa(){
		return accuDa;
	}

}
