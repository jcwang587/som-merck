#!/bin/bash

# Start a local job server
/apps/prod/COMPCHEM/schrodinger/schrodinger/jsc local-server-start

# cd ../data/merck

# Get all the files in the current directory
for file in $(ls); do

    # Check if the file is a .mae file
    if [[ $file == *.mae ]]; then

        # Run the command and wait for it to finish before proceeding to the next file
        /apps/prod/COMPCHEM/schrodinger/schrodinger/epik -cyp -imae $file -omae 3a4_reactivity-${file}

    fi
done

# Stop the local job server
# Do not stop the server because the job server may still be running
# /apps/prod/COMPCHEM/schrodinger/schrodinger/jsc local-server-stop
