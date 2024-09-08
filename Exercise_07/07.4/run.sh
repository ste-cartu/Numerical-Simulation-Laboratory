#!/bin/bash


if [[ $# -ne 1 && $# -ne 2 ]]; then
    echo -e "ERROR! This script needs one or two arguments (NVE/NVT)"
    exit 1
fi

phases=("Solid" "Liquid" "Gas")

for type in "$@"
do
    if [[ "$type" != "NVE" && "$type" != "NVT" ]]; then
        echo "ERROR! Unknown argument!"
        exit 1
    fi
    for phase in "${phases[@]}"; do
        echo
        echo "––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
        echo "${phase} - ${type}"
        make
        ./main ${type} ${phase}
        echo
    done
done