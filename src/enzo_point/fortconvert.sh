#!/bin/bash

# Set the source and destination directories
src_dir="./nbody"
dest_dir="./nbodyf90"

# Ensure the destination directory exists
mkdir -p "$dest_dir"

# Iterate through each .for file in the source directory
for file in "$src_dir"/*.for; do
    # Get the base filename without extension


    base_name=$(basename "$file" .for)
    
    # Create the new filename with .f90 extension
    new_file="$dest_dir/$base_name.f90"
	
    echo $file
    echo $new_file
    
    # Use f2f90 to convert the Fortran 77 code to Fortran 90
    f2f $file $new_file
done

echo "Conversion completed!"

