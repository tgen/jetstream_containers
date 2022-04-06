#!/bin/bash

Rver=3.6.1
INPUT=R_modules.csv
IFS=','
if [ -f modInst.txt ]; then
    rm modInst.txt
fi
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read modName modVer
do
   if [ $modVer == "version" ]; then
      modVer=$Rver
   fi
   modFile=$modName"_"$modVer".tar.gz"
   modURL="https://cran.r-project.org/src/contrib/"$modFile
   if curl --output /dev/null --silent --fail -r 0-0 "$modURL"; then
      echo "Installing $modName version $modVer which is current version of module."
      wget $modURL
      Rscript --no-save install_packages_or_die.R $modFile
      rm $modFile
   else
      modURL="https://cran.r-project.org/src/contrib/Archive/"$modName"/"$modFile
      if curl --output /dev/null --silent --fail -r 0-0 "$modURL"; then
         echo "Installing $modName version $modVer from archive."
         wget $modURL
         Rscript --no-save install_packages_or_die.R $modFile
         rm $modFile
      else
         echo "NO VERSION FOUND!!!!!"
         echo "$modName - $modVer"
         exit 99
      fi
   fi
done < $INPUT
