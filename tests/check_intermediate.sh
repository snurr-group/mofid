#!/bin/bash

# ANSI
GREEN="\033[0;92m"
RED="\033[0;31m"
NC="\033[0m" # No color

for path_to_cif in Resources/KnownCIFs/*.cif; do
    cif="$(basename $path_to_cif | sed 's/.cif//g')"
    echo -e "${GREEN}RUNNING${NC} bin/sbu $path_to_cif Output/$cif..."
    # TODO: make sbu output condition on verbose flag
    bin/sbu $path_to_cif Output/$cif &> /dev/null
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
        # In order, ELEMENT# -> ELEMENT, NODE# -> "", Lines with # -> "", # -> 1 decimal place
        pattern="s/([a-zA-Z]+)[0-9]+/\1/g; s/NODE [0-9]+/NODE/g; s/^[ \t]*#.*//g; s/([0-9]+)(.[0-9])?[0-9]+/\1\2/g"
        for path_to_file in $path_to_dir/*; do
            file="$(basename $path_to_file)"
            diff -b <(sed -E "$pattern" $path_to_file | sort) <(sed -E "$pattern" $output_path/$dir/$file | sort) > /dev/null
            if [[ $? -ne 0 ]]; then
                echo -e "${RED}WARNING:${NC} $dir/$file is different"
                diff -yb --suppress-common-lines <(sed -E "$pattern" $path_to_file | sort) <(sed -E "$pattern" $output_path/$dir/$file | sort) | colordiff
            fi
        done
    done
done
# TODO: remove this once the script is complete
exit 0
