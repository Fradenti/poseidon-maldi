#!/bin/bash

set -e  # stop if any job fails

count=0

for i_s in {1..4}; do
    for i in {1..40}; do
        i_sim="$i" i_str="$i_s" Rscript E1_run_BayesSpace.R &

        ((++count))

        # Every 5 jobs, wait
        if (( count % 5 == 0 )); then
            wait
        fi
    done
done

wait

# In terminal, first make it executable with:
# chmod +x E2_run_BayesSpace_simulations.sh
# Then run it with:
# ./E2_run_BayesSpace_simulations.sh
