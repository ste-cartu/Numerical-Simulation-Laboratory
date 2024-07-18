#!/bin/bash


if [[ $# -ne 1 && $# -ne 2 ]]; then
    echo -e "ERROR! This script needs one or two arguments (metro/gibbs)"
    exit 1
fi

for type in "$@"
do

    if [ $type == "metro" ]; then
        INPATH="Metropolis/INPUT/"
        OUTPATH="Metropolis/OUTPUT/"
    elif [ $type == "gibbs" ]; then
        INPATH="Gibbs/INPUT/"
        OUTPATH="Gibbs/OUTPUT/"
    else 
        echo "ERROR! Unknown argument!"
        exit 1
    fi
    INFILE="${INPATH}input.dat"
    if [[ ! -f "$INFILE" ]]; then
        echo "ERROR! File $INFILE not found!"
        exit 1
    fi



    H=$(printf "%.2f" 0)
    echo
    echo "MAGNETIC FIELD: $H"
    OUTFILE="${OUTPATH}summary_H=${H}.csv"
    touch $OUTFILE
    echo "TEMP,ENERGY,ENERGY_ERR,CV,CV_ERR,CHI,CHI_ERR" > $OUTFILE

    SEARCH="^SIMULATION_TYPE.*$"
    REPLACE="SIMULATION_TYPE        2 1 ${H}"
    sed -i.bak "s|$SEARCH|$REPLACE|g" $INFILE

    mv ${INPATH}properties_H=0.00.dat ${INPATH}properties.dat

    start=0.2
    step=0.04
    end=3.0

    for i in $(seq $start $step $end)
    do    
        TEMP=$(printf "%.3f" "$i")
        echo "TEMPERATURE: $TEMP"
        SEARCH="^TEMP.*$"
        REPLACE="TEMP                   $TEMP"
        sed -i.bak "s|$SEARCH|$REPLACE|g" $INFILE
        if [[ $? -ne 0 ]]; then
            echo "ERROR! Impossible to replace the temperature value in file input.csv!"
        fi
        make $type

        ENERGY=$(cut -d ',' -f 3 ${OUTPATH}SIMULATION/total_energy_H=${H}_t=${TEMP}.csv | tail -n 1)
        ENERGY_ERR=$(cut -d ',' -f 4 ${OUTPATH}SIMULATION/total_energy_H=${H}_t=${TEMP}.csv | tail -n 1)
        CV=$(cut -d ',' -f 3 ${OUTPATH}SIMULATION/specific_heat_H=${H}_t=${TEMP}.csv | tail -n 1)
        CV_ERR=$(cut -d ',' -f 4 ${OUTPATH}SIMULATION/specific_heat_H=${H}_t=${TEMP}.csv | tail -n 1)
        CHI=$(cut -d ',' -f 3 ${OUTPATH}SIMULATION/susceptibility_H=${H}_t=${TEMP}.csv | tail -n 1)
        CHI_ERR=$(cut -d ',' -f 4 ${OUTPATH}SIMULATION/susceptibility_H=${H}_t=${TEMP}.csv | tail -n 1)

        echo "${TEMP},${ENERGY},${ENERGY_ERR},${CV},${CV_ERR},${CHI},${CHI_ERR}" >> $OUTFILE
    
    done

    mv ${INPATH}properties.dat ${INPATH}properties_H=${H}.dat

    H=$(printf "%.2f" 0.02)
    echo
    echo
    echo "MAGNETIC FIELD: $H"
    OUTFILE="${OUTPATH}summary_H=${H}.csv"
    touch $OUTFILE
    echo "TEMP,ENERGY,ENERGY_ERR,MAGNET,MAGNET_ERR" > $OUTFILE

    SEARCH="^SIMULATION_TYPE.*$"
    REPLACE="SIMULATION_TYPE        2 1 ${H}"
    sed -i.bak "s|$SEARCH|$REPLACE|g" $INFILE

    mv ${INPATH}properties_H=0.02.dat ${INPATH}properties.dat

    for i in $(seq $start $step $end)
    do
        echo 
        TEMP=$(printf "%.3f" "$i")
        echo "TEMPERATURE: $TEMP"
        SEARCH="^TEMP.*$"
        REPLACE="TEMP                   $TEMP"
        sed -i.bak "s|$SEARCH|$REPLACE|g" $INFILE
        if [[ $? -ne 0 ]]; then
            echo "ERROR! Impossible to replace the temperature value in file input.csv!"
        fi
        make $type

        ENERGY=$(cut -d ',' -f 3 ${OUTPATH}SIMULATION/total_energy_H=${H}_t=${TEMP}.csv | tail -n 1)
        ENERGY_ERR=$(cut -d ',' -f 4 ${OUTPATH}SIMULATION/total_energy_H=${H}_t=${TEMP}.csv | tail -n 1)
        MAGNET=$(cut -d ',' -f 3 ${OUTPATH}SIMULATION/magnetization_H=${H}_t=${TEMP}.csv | tail -n 1)
        MAGNET_ERR=$(cut -d ',' -f 4 ${OUTPATH}SIMULATION/magnetization_H=${H}_t=${TEMP}.csv | tail -n 1)

        echo "${TEMP},${ENERGY},${ENERGY_ERR},${MAGNET},${MAGNET_ERR}" >> $OUTFILE
    
    done

    mv ${INPATH}properties.dat ${INPATH}properties_H=${H}.dat

done
