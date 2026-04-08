#!/bin/bash

set -e  # stop if any job fails

source .venv/bin/activate

count=0

for i_str in {1..4}; do
    for i in {1..40}; do # it was 25
        python3 D1_run_UTAG.py --i_sim $i --i_str $i_str &

        ((++count))

        # Every 5 jobs, wait
        if (( count % 5 == 0 )); then
            wait
        fi
    done
done

wait

# In terminal, first make it executable with:
# chmod +x D2_run_python_simulations.sh
# Then run it with:
# ./D2_run_python_simulations.sh
