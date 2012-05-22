package org.sams.spect;

/**
 * This class express a scan of MZXML format
 * 
 * @author Miguel Rojas-Cherto
 */
public class ScanMZXML {

	private String num;
	private String polarity;
	private String collisionEnergy;

	/**
	 * Set the parameters of a scan
	 * 
	 * @param num        The number of the scan
	 * @param polarity   The polarity
	 * @param collisionEnergy The collisionEnergy
	 */
	public ScanMZXML(String num, String polarity, String collisionEnergy) {
		this.num = num;
		this.polarity = polarity;
		this.collisionEnergy = collisionEnergy;
	}
	/**
	 * Get the number of the scan
	 * 
	 * @return The Number of the scan
	 */
	public String getNum(){
		return num;
	}

	/**
	 * Get the polarity of the scan
	 * 
	 * @return The polarity of the scan
	 */
	public String getpolarity(){
		return polarity;
	}
	/**
	 * Get the collision Energy of the scan
	 * 
	 * @return The collision Energy of the scan
	 */
	public String getcollisionEnergy(){
		return collisionEnergy;
	}

}
