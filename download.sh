#!/bin/bash

# Download everything
wget -r -l 1 http://irimodel.org/IRI-2016/00_iri2016.tar
wget -r -l 1 http://irimodel.org/IRI-2016/00_iri2016-License.txt
wget -r -l 1 http://irimodel.org/IRI-2016/00_readme_IRI-2016.txt
wget -r -l 1 http://irimodel.org/IRI-2016/00_update_history.txt

# Copy all the text files
mv irimodel.org/IRI-2016/*.txt ./

# unpack the archive into appropriate locations
tar xf irimodel.org/IRI-2016/00_iri2016.tar -C data/ --wildcards "*.dat"
tar xf irimodel.org/IRI-2016/00_iri2016.tar -C src/ --wildcards "*.for"
tar xf irimodel.org/IRI-2016/00_iri2016.tar -C src/ --wildcards "*.txt"

# clean up the download directory
rm -rf irimodel.org

# record the date
date >> downloadDate.txt
