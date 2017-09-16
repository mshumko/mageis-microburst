#!/bin/bash

# This script will compress figure 1 and 2, and replace the files.

gs -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -sOutputFile=fig1_compressed.pdf fig1.pdf

gs -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -sOutputFile=fig2_compressed.pdf fig2.pdf
