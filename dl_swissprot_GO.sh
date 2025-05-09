#!/bin/bash

set -e

for YEAR in $(seq 2010 2024); do
    REL="release-${YEAR}_01"
    FILE="uniprot_sprot-only${YEAR}_01.tar.gz"
    URL="https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/${REL}/knowledgebase/${FILE}"
    OUTDIR="${YEAR}_01"
    TSV="swissprot_${YEAR}_01_annotations.tsv"

    # Download if not already present
    if [ ! -f "$FILE" ]; then
        echo "Downloading $FILE..."
        wget -q "$URL" || { echo "Failed to download $URL"; continue; }
    fi

    # Extract tar.gz
    if [ ! -d "$OUTDIR" ]; then
        mkdir "$OUTDIR"
        tar -xzf "$FILE" -C "$OUTDIR"
    fi

    # Find and gunzip uniprot_sprot.dat.gz
    DATGZ=$(find "$OUTDIR" -name "uniprot_sprot.dat.gz" | head -n 1)
    if [ -z "$DATGZ" ]; then
        echo "No uniprot_sprot.dat.gz found for $YEAR"
        continue
    fi
    DAT="${DATGZ%.gz}"
    if [ ! -f "$DAT" ]; then
        gunzip -k "$DATGZ"
    fi

    # Parse uniprot_sprot.dat
    awk '
    # Set the record separator to "//\n" so each UniProt entry is a record
    BEGIN {
        RS = "//\n";      # Each entry ends with "//"
        OFS = "\t";        # Output fields separated by tab
    }
    {
        id = "";           # Variable to store the UniProt ID
        go = "";           # Variable to store GO terms (semicolon-separated)
        # Loop through each line of the current record
        for (i = 1; i <= NF; i++) {
            # If the line starts with "ID   ", extract the protein ID (second field)
            if ($i ~ /^ID   /) {
                split($i, a, " ");  # Split line by spaces
                id = a[2];          # The second field is the UniProt ID
            }
            # If the line starts with "DR   GO;", extract the GO term
            if ($i ~ /^DR   GO;/) {
                match($i, /GO:[0-9]+/, m);  # Find GO term using regex
                if (m[0] != "") {
                    # If go already has a value, append with semicolon
                    if (go != "") go = go ";" m[0];
                    else go = m[0];
                }
            }
        }
        # If both ID and at least one GO term were found, print them
        if (id != "" && go != "") print id, go;
    }
    ' "$DAT" > "$TSV"
    
    rm *.gz *.dat
    echo "Processed $YEAR: $TSV"
done
