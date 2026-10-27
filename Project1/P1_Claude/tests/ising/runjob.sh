#!/bin/bash
# usage: runjob.sh "name|namelist body"
cd /Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/ising
name="${1%%|*}"; body="${1#*|}"
if [ -f "runs/$name.dat" ] && grep -q "cpu_s" "runs/$name.dat"; then exit 0; fi
start=$(date +%s)
echo "&run $body outfile='runs/$name.dat' /" | ./ising > "runs/$name.log" 2>&1
echo "$name done in $(( $(date +%s) - start )) s"
