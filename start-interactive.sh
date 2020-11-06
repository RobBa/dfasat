#!/bin/bash

if [ $# -ne 2 ]; then
    echo "usage: $0 ini-file inputfile"
    exit 1
fi

#echo $param

cp tests/template.pdf tests/pre_0.dot.pdf
cp tests/template.pdf tests/pre_1.dot.pdf

evince tests/pre_0.dot.pdf &>/dev/null &
evince tests/pre_1.dot.pdf &>/dev/null &

sleep 5

./interactive.sh &>/dev/null &

echo "Starting interactive mode"

./flexfringe --ini="$1" "$2"

evince tests/final.dot.pdf &>/dev/null

trap "trap - SIGTERM && kill -- -$$" SIGINT SIGTERM EXIT

