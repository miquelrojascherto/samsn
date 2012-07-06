#!/bin/sh

# **********************************************
# script to install MEF/SAMSn dependencies
# **********************************************

lowercase(){
    echo "$1" | sed "y/ABCDEFGHIJKLMNOPQRSTUVWXYZ/abcdefghijklmnopqrstuvwxyz/"
}

######### Store the current working directory to be able to return here
CWD=`pwd`

######### Detect OS
OS=`lowercase \`uname\``
echo "${OS}"
if [ "{$OS}" == "windowsnt" ]; then
	#better switch to Mac or Linux!
	echo "\nYou OS is not supported"
	exit 1
fi

echo "\nYou are running $OS". 


######### Check if Curl is installed
which curl &> /dev/null
if [ $? -eq 1 ]; then
    echo >&2 "\nCurl not found. Make sure Curl (http://curl.haxx.se/download.html) is installed"
    exit 1
fi


######### Check if Java is installed
which java &> /dev/null
if [ $? -eq 1 ]; then
    echo >&2 "\nJava not found. Make sure Java (http://java.com) is installed"
    exit 1
fi

######### Check if R is installed
which R &> /dev/null
if [ $? -eq 1 ]; then
	echo >&2 "\nR not found. Trying to build R (http://r-project.org) from source"

	#Download source (2.15)
	curl -O http://cran-mirror.cs.uu.nl/src/base/R-2/R-2.15.0.tar.gz
	tar -zxvf R-2.15.0.tar.gz
	cd R-2.15.0
	
	# configure for mac or linux
	if [ "{$OS}" == "darwin"]; then
		./configure --with-blas='-framework vecLib' --with-lapack --with-aqua --enable-R-framework
	else
		./configure
	fi
	
	# run make and install
	make
	make install
	cd $CWD
fi

######### Check if Graphviz is installed
which dot &> /dev/null
if [ $? -eq 1 ]; then
    echo >&2 "\nGraphviz not found. Make sure Graphviz (http://graphviz.org) is installed"
    
    #Download source (2.28)
    curl -O http://graphviz.org/pub/graphviz/stable/SOURCES/graphviz-2.28.0.tar.gz
    tar -zxvf graphviz-2.28.0.tar.gz
    cd graphviz-2.28.0
    ./configure
    make
    make install
	cd $CWD
fi

######### Now to make sure we have the required R packages installed
echo "\nInstalling xcms and Rgraphviz"
R --slave <<EOF
	source("http://bioconductor.org/biocLite.R")
	biocLite("xcms")
	biocLite("Rgraphviz")
EOF

######### And now Rserve
echo "\nInstalling Rserve"
curl -O http://www.rforge.net/Rserve/snapshot/Rserve_1.7-0.tar.gz
R CMD INSTALL Rserve_1.7-0.tar.gz








