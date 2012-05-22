package org.sams;
/**
 * An interface providing predefined values for a number of
 * constants used throughout the NMC.
 * 
 * @author Miguel Rojas-Cherto
 */
public class MZDataConstants {

    /** Description of the fragmentation tree type. */
	public static final String FRAGMENTS = "nmc:fragments";
    public static final String LOSSES = "nmc:losses";

    /** Attributes of CML*/
    public static final String INSTRUMENT = "nmc:instrumentRef";
    public static final String PARENT_FILE = "nmc:parentFile";
    public static final String PROTOCOL = "nmc:protocol";
    public static final String NUM_GROUPS = "nmc:numGroups";
    public static final String SCANS = "nmc:countScans";
    public static final String CREATED = "nmc:created";
    public static final String COMMENT = "nmc:comment";

    /** Attributes of CML spectrum*/
    public static final String SCAN_NUM = "nmc:scanNum";
    public static final String MS_LEVEL = "nmc:msLevel";
    public static final String RT = "nmc:rt";
    public static final String GROUP_PEAK_MSN = "nmc:groupPeakMSn";
    public static final String PRECURSOR_MZ = "nmc:precursorMZ";
    public static final String PRECURSOR_SCAN = "nmc:precursorScan";
    public static final String WINDOW = "nmc:window";
    public static final String COLLISION_ENERGY = "nmc:collisionEnergy";
    public static final String ACTIVATION_METHOD = "nmc:activationMethod";
    public static final String POLARITY = "nmc:polarity";
    public static final String MASS = "nmc:mass";
    public static final String INTENSITIY = "nmc:intensity";
    
    /** Parameters processing MZXML*/
    public static final String MZGAP = "nmc:mzgap";
    public static final String SNTHRESH = "nmc:snthresh";
    public static final String RINT = "nmc:rint";
    public static final String VERSION_XCMS = "nmc:XCMSversion";
    
    /** Parameters enriching CML*/
    public static final String ACCURACY = "nmc:accuracy";
    public static final String RULES = "nmc:rules";
    public static final String ELEMENTS = "nmc:elements";
    public static final String VERSION_SAMS = "nmc:SAMSversion";
    public static final String OCCURRENCE = "nmc:occurrence";

    /** Compound Identifier CML*/
    public static final String PARENT_COMPOUND = "nmc:parentCompounds";
    public static final String PRECURSOR_ID = "nmc:precursorId";
    
    /** Compound Attributes*/
    public static final String COLOR = "nmc:color";
    
}
