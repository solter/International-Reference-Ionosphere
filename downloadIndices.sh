#!/bin/bash

wget -r -l 1 http://irimodel.org/indices/
date >> indexDownloadDate.txt
