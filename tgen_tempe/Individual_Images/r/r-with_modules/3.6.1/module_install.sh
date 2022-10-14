#!/bin/bash

function checkCRAN() {
    modURLcheck="https://cran.r-project.org/src/contrib/"$modFile
    if curl --output /dev/null --silent --fail -r 0-0 "$modURLcheck"; then
        echo "$modURLcheck"
    else
        modURLcheck="https://cran.r-project.org/src/contrib/Archive/"$modName"/"$modFile
        if curl --output /dev/null --silent --fail -r 0-0 "$modURLcheck"; then
            echo "$modURLcheck"
        else
            echo "none"
        fi
    fi
}

function installCranMod() {
    wget $modURL
    if [ $modArgs ]; then
        R CMD INSTALL ${modFile} --configure-args="${modArgs}" || { echo "Install failed for module: ${modFile}"; rm ${modFile}; exit 1; }
    else
        R CMD INSTALL ${modFile} || { echo "Install failed for module: ${modFile}"; rm ${modFile}; exit 1; }
    fi
    rm ${modFile}
}

function installBioConductor() {
    Rscript -e 'BiocManager::install('\""${modName}"\"')' || { echo "Install failed for module: ${modName}"; rm ${modFile}; exit 1; }
    echo "Installed BioConductor Module ${modName}"
}

Rver=3.6.1
INPUT=R_modules.csv
IFS=','
if [ -f modInst.txt ]; then
    rm modInst.txt
fi
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read modName modVer modArgs
do
    if [ $modVer == "version" ]; then
        modVer=$Rver
    fi
    modFile=$modName"_"$modVer".tar.gz"
    modURL=$(checkCRAN "$modFile")
    if [ $modURL == "none" ]; then
        if [ -f ${modFile} ]; then
            installCranMod "$modURL" "$modFile" "$modArgs"
            echo "Installed mod ${modName} from local file."
        else
            installBioConductor "$modName"
        fi
    else
        echo "Installing from CRAN."
        installCranMod "$modURL" "$modFile" "$modArgs"
    fi
done < $INPUT