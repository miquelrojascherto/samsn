MEF - The Molecular Elemental Formula
 
Copyright 2011 The Netherland Metabolomics Center Team License: GPL v3, see
doc/gpl.license

1. Introduction
----------------------------------------------------------------------------------
You are currently reading the README file for the MEF Project.This project will
be hosted very soon under http://sourceforge.net/.

The MEF is an open-source tool for the assignment of chemical structure information to
mass spectrometry data, implemented in the programming language Java(tm). The
library is published under terms of the standard The GNU General Public License 
(GPL) v3. This has implications on what you can do with sources and binaries of the
MEF library. For details, please refer to the file LICENSE, which should have
been provided with this distribution.

PLEASE NOTE:  is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.

2. Prerequisites
----------------------------------------------------------------------------------
Before to be able to run successful the MEF tool you need a working installation
of the following other tools:

 - java 1.5 or higher

 - R/Bioconductor Project for Statiscal Computing. http://www.r-project.org/
 
 - Rserve R library. More information you can find here: 
 http://www.rforge.net/Rserve/doc.html
  
 - XCMS library. More information you can find here: http://metlin.scripps.edu/xcms/. 
 You may install bioconductor's XCMS library just issuing the following within your
 R console:
 
		source("http://bioconductor.org/biocLite.R")
		biocLite("xcms")

 - Graphviz library
 		
 		biocLite("Rgraphviz")
 		sudo apt-get install libgraphviz-dev


3. Using MEF tool
----------------------------------------------------------------------------------
In order to use the MEF tool in your program, you need to run the jar file. An
example using the command line could be that: 
USER/$ R CMD Rserve
USER/$ java -jar sams-0.1.git.jar -occurr=0.4 -sn=1 -mzgap=0.5 -rint=0 -acc=5 -e=[C1..15,H1..9,O0..4,N0..2] -rules=[RDBER,nitrogenR] -imzXML  F002169_C15H9O4.mzXML -ocml  F002169_C15H9O4.cml process

After the process is finished it will be created a cml enriched with chemical
information together with the spectra.

Command to create a pdf file which show the fragmentation tree. 
USER/$ java -jar sams-0.1.git.jar -icml F002169_C15H9O4.cml -opdf F002169_C15H9O4.pdf convert

Command to compare two fragmentation trees. 
USER/$ java -jar sams-0.1.git.jar -i1cml F002169_C15H9O4.cml -i2cml  F002169_C15H9O4.cml compare


4. Help
----------------------------------------------------------------------------------
If you need help don't hesitate to contact us (m.rojas@lacdr.leidenuniv.nl)

----------------------------------------------------------------------------------
Enjoy!  Comments and feedback are appreciated!

Miguel Rojas Cherto
m.rojas@lacdr.leidenuniv.nl
Leiden, Netherlands
http://abs.lacdr.gorlaeus.net/people/rojas-cherto
