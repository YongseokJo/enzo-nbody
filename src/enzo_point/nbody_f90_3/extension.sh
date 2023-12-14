#!/bin/bash

# Set the absolute source and destination directories
src_dir="."

# Iterate through each .for file in the source directory
for file in "$src_dir"/*.for; do
    # Get the base filename without extension
    base_name=$(basename "$file" .for)
    
    # Create the new filename with .f90 extension
    new_file="$src_dir/$base_name.f90"
    
    # Use f2f90 to convert the Fortran 77 code to Fortran 90
    mv $file $new_file
done

echo "Conversion completed!"

