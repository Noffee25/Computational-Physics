#!/bin/bash
# worker: each line of the queue file is "tag gen seed single plan"
cd /Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/spacing
while read tag gen seed single plan; do
  ./spacing_scan $gen $seed $single $plan scan/$tag.txt >> scan/$tag.log 2>&1
done < $1
