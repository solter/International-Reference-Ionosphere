#!/bin/bash

wget -r -l 1 http://irimodel.org/indices/
mv irimodel.org/indices/*.dat indices/
rm -rf irimodel.org
date >> indexDownloadDate.txt
