#!/bin/bash

# This script will compress figure 1 and 2, and replace the files.

gs -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -sOutputFile=figure_scripts/fig1_compressed.pdf figure_scripts/fig1.pdf

gs -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -sOutputFile=figure_scripts/fig2_compressed.pdf figure_scripts/fig2.pdf
