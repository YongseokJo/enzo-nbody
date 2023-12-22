#!/bin/bash

# Specify the directory where your source code files are located
directory="."

# Use the find command to locate all relevant files
files=$(find "$directory" -name "*.for")

# Loop through each file and use sed to replace the text
for file in $files
do
    # Use sed to replace "include common6.h" with "USE MODULE"
    sed -i "s/INCLUDE 'common6.h'/USE POINTERS\n      INCLUDE 'common6.h'/g" "$file"

done

echo "Replacement complete"
