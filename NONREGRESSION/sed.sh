#!/bin/bash

# Loop over directories starting with CAS_
for dir in CAS_*; do
    # Ensure it's a directory and that Donnees.arc exists
    if [[ -d "$dir" && -f "$dir/Donnees.arc" ]]; then
        echo "Processing $dir/Donnees.arc"

        # Comment out the <period>...</period> line
        sed -i 's|^\s*<output-period>.*</output-period>|<!-- & -->|' "$dir/Donnees.arc"
    fi
done

echo "Done."
