#!/bin/bash

# ANSI
GREEN="\033[0;92m"
RED="\033[0;31m"
NC="\033[0m" # No color

# (Remove lines startign with "#" as these are comments?)
# Sort both files and check diff
for path_to_cif in Resources/KnownCIFs/*.cif; do
    cif="$(basename $path_to_cif | sed 's/.cif//g')"
    bin/sbu $path_to_cif Output/$cif
    echo "\n"
done
for path_to_cif in Resources/KnownCIFs/*.cif; do
    cif="$(basename $path_to_cif | sed 's/.cif//g')"
    echo -e "${GREEN}CHECKING${NC} $cif..."
    output_path="Output/$cif"
    known_path="Resources/KnownCIFs/Outputs/$cif"
    # Check orig_mol.cif
    diff <(sort $known_path/orig_mol.cif) <(sort $output_path/orig_mol.cif) > /dev/null
    if [[ $? -ne 0 ]]; then
        echo -e "${RED}WARNING:${NC} orig_mol.cif is different" 
    fi
    for path_to_dir in $(find $known_path/* -type d); do
        dir="$(basename $path_to_dir)"
        for path_to_file in $path_to_dir/*; do
            file="$(basename $path_to_file)"
            #echo $output_path/$dir/$file
            diff <(sort $path_to_file) <(sort $output_path/$dir/$file) > /dev/null
            if [[ $? -ne 0 ]]; then
                echo -e "${RED}WARNING:${NC} $dir/$file is different"
                diff -y <(sort $path_to_file) <(sort $output_path/$dir/$file)
            fi
        done
    done
done

