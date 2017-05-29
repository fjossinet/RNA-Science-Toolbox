#!/bin/bash


for f in $(more $1 | sed -r 's/.+jobs\/job_(.+).jdl .+/\1/')
do
	more script_$f.sh | grep "#RFAM families processed" | sed -r 's/.+processed: (.+)/\1,/' | tr -d '\n'
done
