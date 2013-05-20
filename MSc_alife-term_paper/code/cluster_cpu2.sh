#!/bin/bash
for i in $(seq 6 8) ; do 
	./daisies -v 2 --mut_rate 0.$i --num_runs 100 1> /dev/null 2> /tmp/lg80-daisies/100runs_mut0.$i 
done 
echo "third core done" | mail -s "sim report" lorenzo.grespan@gmail.com
