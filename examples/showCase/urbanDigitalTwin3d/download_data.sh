#!/usr/bin/env bash

# Base URL of the files
base_url="https://openlb.net/data/city_showcase_osm"

# List of files to download
files=(
    "winddataReutlingen6_14.CSV"
    "map.osm"
    "map_small.osm"
    "station_2024-11-06-2024-11-13_lederStrasse.csv"
    "station_2024-11-06-2024-11-13_alteBurgStrasse.csv"
)

# Target directory
#target_dir=""
#mkdir -p "$target_dir"

# Download each file
for file in "${files[@]}"; do
    echo "Downloading $file..."
    wget --user-agent="Mozilla/5.0" -O "$file" "$base_url/$file"
    
    if [ $? -eq 0 ]; then
        echo "✓ $file downloaded successfully."
    else
        echo "✗ Error downloading $file"
    fi

    sleep 2 # ⏳ Wait 1 second between requests to avoid rate limiting
done

