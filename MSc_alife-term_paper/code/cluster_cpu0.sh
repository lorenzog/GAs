#!/bin/bash
for i in $(seq 0 2) ; do 
	./daisies -v 2 --mut_rate 0.$i --num_runs 100 1> /dev/null 2> /tmp/lg80-daisies/100runs_mut0.$i 
done 
echo "first core done" | mail -s "sim report" lorenzo.grespan@gmail.com
