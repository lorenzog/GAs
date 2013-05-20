#!/bin/bash
./daisies -v 2 --mut_rate 0.9 --num_runs 100 1> /dev/null 2> /tmp/lg80-daisies/100runs_mut0.9
./daisies -v 2 --mut_rate 1.0 --num_runs 100 1> /dev/null 2> /tmp/lg80-daisies/100runs_mut1.0

echo "fourth core done" | mail -s "sim report" lorenzo.grespan@gmail.com
