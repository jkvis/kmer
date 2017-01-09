#!/bin/bash

getopt --test > /dev/null
if [[ $? -ne 4 ]]; then
    echo "`getopt --test` failed in this environment."
    exit 1
fi

SHORT=?hk:
LONG=help,kmer:

# -temporarily store output to be able to check for errors
# -activate advanced mode getopt quoting e.g. via “--options”
# -pass arguments only via   -- "$@"   to separate them correctly
PARSED=`getopt --options $SHORT --longoptions $LONG --name "$0" -- "$@"`
if [[ $? -ne 0 ]]; then
    # e.g. $? == 1
    #  then getopt has complained about wrong arguments to stdout
    exit 2
fi
# use eval with "$PARSED" to properly handle the quoting
eval set -- "$PARSED"

k=2

# now enjoy the options in order and nicely split until we see --
while true; do
    case "$1" in
        -h|--help)
            echo "usage: $0"
            exit 0
            ;;
        -k|--kmer)
            k="$2"
            shift 2
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Programming error"
            exit 3
            ;;
    esac
done

# handle non-option arguments
if [[ $# -ne 1 ]]; then
    echo "$0: A single input file is required."
    exit 4
fi
extension="${1##*.}"

if [[ "$extension" == "gz" ]]; then
    zcat $1 | awk '{if(++count%4==2) print $0;}' | tr -c 'ACGT.' 'ACGTN' | tr -s 'N' '\n' | awk '{if(length($0)>'$k') print $0;}'
else
    cat $1 | awk '{if(++count%4==2) print $0;}' | tr -c 'ACGT.' 'ACGTN' | tr -s 'N' '\n' | awk '{if(length($0)>'$k') print $0;}'
fi

