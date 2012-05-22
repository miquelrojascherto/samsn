package org.sams.main;

import java.io.IOException;
import java.io.InputStream;
import java.util.jar.Attributes;
import java.util.jar.Manifest;

/**
 * Class that return package version and more.
 * 
 * @author Miguel Rojas-Cherto
 *
 */
public class PackageVersion {
	
	Attributes attributes;
	public static String ATTRIBUTE_IMPLEMENTATIONVERSION = "Implementation-Version";
	public static String ATTRIBUTE_BUILTDATE = "Built-Date";
	
	public PackageVersion(){
		try{         
	        InputStream stream = ProcessMZData.class.getResourceAsStream("/META-INF/MANIFEST.MF2");
	
		        if (stream == null)
		            System.out.println("Couldn't find manifest.");
		        else{
			        Manifest manifest = new Manifest(stream);
			        attributes = manifest.getMainAttributes();
		        }
        }catch (IOException e){            
            System.out.println("Couldn't read manifest.");
        }     
	}
	public void printVersion()
    {
		if(attributes != null){
            String impVersion = attributes.getValue("Implementation-Version");
            String impBuildDate = attributes.getValue("Built-Date");
            
            if (impVersion != null)
            {
                System.out.println("Implementation-Version: " + impVersion);
            }
            if (impBuildDate != null)
            {
                System.out.println("Built-Date: " + impBuildDate);
            }
		}
    }
	
	public String getAttribute(String attribute){
		return attributes.getValue(attribute);
	}
	/**
	 * Check if it exist the manifest
	 * @return
	 */
	public boolean existManifest() {
		if(attributes != null)
			return true;
		else
			return false;
	}
}