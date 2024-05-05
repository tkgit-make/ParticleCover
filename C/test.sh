#!/bin/bash

# accepts one argument, which controls profiling options:
# - gprof: compiles profiler options into command and uses gprof to build a report
# - time: uses `time` program for profiling (apt install time)
# - stderr: redirects stdout to /dev/null to show only stderr
# compiles all C files into the ProcessInput executable. assumes presence of main function  
# ./test.sh gprof, or ./test.sh time

OPTS=""
INVOKER=""

if [[ $1 == "gprof" ]]; then
	OPTS="-pg"
elif [[ $1 == "time" ]]; then
	INVOKER="/usr/bin/time -v -o bin/time.results.txt"
elif [[ $1 == "stderr" ]]; then
    OUTPUT="/dev/null"
fi

mkdir -p bin
rm -f bin/ProcessInput
if ! gcc -g $OPTS C/*.c -o bin/ProcessInput -lm; then
	echo "compile fail"
	exit
fi
#zcat C/wedgedata.txt.gz | $INVOKER bin/ProcessInput > $OUTPUT
cat wedgeData_v3_128.txt | $INVOKER bin/ProcessInput > $OUTPUT
if [[ $1 == "gprof" ]]; then
	date > bin/gprof.results.txt
	gprof bin/ProcessInput gmon.out >> bin/gprof.results.txt
fi
