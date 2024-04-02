#!/bin/bash

set -uex

# Check for correct usage
if [ $# -ne 1 ]; then
    echo "Usage: $0 listOfCodes"
    exit 1
fi

# Store input files with descriptive names
LIST_OF_CODES=$1


echo "Starting Epigenetic marks download..."

for file in $(cat "$LIST_OF_CODES"); do
    mkdir -p ../"${file%/*}"  # Create directory if not present
    code=$(basename "$file" .bed)  # Strip the .bed extension
    echo "Downloading " "$code"
    if wget -q https://www.encodeproject.org/files/"$code"/@@download/"$code".bed.gz -O ../"$file".gz; then
        echo "Download successful: " "$code"
    else
        echo "Download failed." "$code"
    fi
    
done
echo "Finished downloading Epigenetic marks"

