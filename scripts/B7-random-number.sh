#!/bin/bash

input_file=$1

output_directory=data/neutral
mkdir -p "$output_directory"

file_name=${output_directory}/$(basename "${input_file%.*}" | sed 's/-dataset//')
output_file="$file_name"-random.csv 


### RANDOM NUMBER ###

total="$(wc -l < ${input_file})"

echo 'Random' > "$output_file" 

i=1
while [ "$i" -lt "$total" ]; do 
  echo $RANDOM >> "$output_file"
  (( i++))
done 

