#!/bin/bash
#./C/test.sh 2>&1 1>/dev/null (only stderr)
#./C/test.sh 2>/dev/null (only stdout)
#./C/test.sh (both)
# accepts one argument, which controls profiling options:
# - gprof: compiles profiler options into command and uses gprof to build a report
# - time: uses `time` program for profiling (apt install time)
# compiles all C files into the ProcessInput executable. assumes presence of main function  
# ./test.sh gprof, or ./test.sh time
# Diff output file in main and dev branch: git diff main..devFPGA -- C/cOutput.txt

OPTS=""
INVOKER=""

if [[ $1 == "gprof" ]]; then
	OPTS="-pg"
elif [[ $1 == "time" ]]; then
	INVOKER="/usr/bin/time -v -o bin/time.results.txt"
fi

mkdir -p bin
rm -f bin/ProcessInput
if ! gcc -g $OPTS C/*.c -o bin/ProcessInput -lm; then
	echo "compile fail"
	exit
fi
gunzip -c C/wedgedata.txt.gz | $INVOKER bin/ProcessInput
# zcat C/wedgedata.txt.gz | $INVOKER bin/ProcessInput
# cat wedgedata-sample.txt | $INVOKER bin/ProcessInput
if [[ $1 == "gprof" ]]; then
	date > bin/gprof.results.txt
	gprof bin/ProcessInput gmon.out >> bin/gprof.results.txt
fi