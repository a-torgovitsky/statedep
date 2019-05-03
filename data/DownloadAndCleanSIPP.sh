#!/bin/bash
echo "Downloading, extracting and converting raw data from NBER..."
for (( wave = 1; wave <= 16; wave ++ ))
do
    if test -f "sippl08puw${wave}.dta"; then
        echo "Wave ${wave} .dta already exists; skipping."
    else
        wget -nc "http://www.nber.org/sipp/2008/l08puw${wave}.zip"
        unzip -o l08puw${wave}.zip
        rm l08puw${wave}.zip
        echo "Running Stata to convert wave ${wave}..."
        stata -b sippl08puw${wave}.do
        echo "Done"
        rm l08puw${wave}.dat sippl08puw${wave}.log
    fi
done

echo "Cleaning"
stata -b CleanSIPP.do

if [ "$1" != "nocleanup" ]; then
    echo "Removing SIPP data"
    rm sipp08_w*.dta sippl08puw*.dta sipp08_merged.dta sipp08_clean.dta
fi

echo "Done"
