#!/bin/bash

# ANSI
GREEN="\033[0;92m"
RED="\033[0;31m"
NC="\033[0m" # No color

# FLAGS
verbose=false

while getopts "v" flag; do
    case $flag in
        v) verbose=true ;;
    esac
done

for path_to_cif in Resources/KnownCIFs/*.cif; do
    cif="$(basename $path_to_cif | sed 's/.cif//g')"
    echo -e "${GREEN}RUNNING${NC} bin/sbu $path_to_cif Output/$cif..."
    if [[ $verbose = true ]]; then
        bin/sbu $path_to_cif Output/$cif 
        echo "\n"
    else
        bin/sbu $path_to_cif Output/$cif &> /dev/null
    fi
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
        # In order, a-zA-Z# -> a-zA-Z, NODE# -> "", Lines with # -> "", # -> 1 decimal place
        pattern="s/([a-zA-Z]+)[0-9]+/\1/g; s/NODE [0-9]+/NODE/g; s/^[ \t]*#.*//g; s/([0-9]+)(.[0-9])?[0-9]+/\1\2/g"
        for path_to_file in $path_to_dir/*; do
            file="$(basename $path_to_file)"
            diff -b <(sed -E "$pattern" $path_to_file | sort) <(sed -E "$pattern" $output_path/$dir/$file | sort) > /dev/null
            exit_code="$?"
            if [[ $exit_code -ne 0 ]] && [[ $verbose = true ]] ; then
                echo -e "${RED}WARNING:${NC} $cif/$dir/$file is different"
                diff -yb --suppress-common-lines <(sed -E "$pattern" $path_to_file | sort) <(sed -E "$pattern" $output_path/$dir/$file | sort) | colordiff
            elif [[ $exit_code -ne 0 ]]; then
                # This pattern tests simpler equality
                simple_pattern="s/[ \t0-9.\-]e?//g; s/^[ \t]*#.*//g;"
                diff -b <(sed -E "$simple_pattern" $path_to_file | sort) <(sed -E "$simple_pattern" $output_path/$dir/$file | sort) > /dev/null
                if [[ $? -ne 0 ]]; then
                    diff -yb --suppress-common-lines <(sed -E "$simple_pattern" $path_to_file | sort) <(sed -E "$simple_pattern" $output_path/$dir/$file | sort) | colordiff
                else
                    echo -e "${GREEN}SUCCESS:${NC} $cif/$dir/$file is identical"
                fi
            else
                echo -e "${GREEN}SUCCESS:${NC} $cif/$dir/$file is identical"
            fi
        done
    done
done
exit 0
