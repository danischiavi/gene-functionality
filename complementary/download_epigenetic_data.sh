#!/bin/bash
set -uex

echo "Starting Epigenetic marks download..."

for file in $(find ../data/epigenetic_data/ -type f -print | grep .txt); do
    mkdir -p "${file%.txt}"  # Create directory if not present

    for link in $(cat "$file" | grep bed.gz); do
        code=$(basename "$link" .bed.gz)
        if [ ! -f "${file%.txt}"/"$code".bed.gz ]; then
            wget -q "$link" -O "${file%.txt}"/"$code".bed.gz
            echo "Download successful: " "$code"
        else
            echo "File present: $code"
        fi
    done
    
done
echo "Finished downloading Epigenetic marks."

