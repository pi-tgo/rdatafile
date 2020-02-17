# rdatafile
Simple minimalistic offline plotting of CADI datafiles

This code is written for python >= 3.6 and require matplotlib and numpy.

Usage: rdatafile.py inputfile outputfile
inputfile = file name of the CADI data input file
outputfile = file name of the resulting plot

Outputfile should have extension .png to produce png-files or .pdf to produce pdf-files.
Supported formats eps, pdf, pgf, png, ps, raw, rgba, svg, svgz (depends on version of matplotlib)

The code is written based on IDL code provided by Chris Meek
IDL code originally written by Ian Grant and modified by the same and JWM

An example datafile (testdata.md49) is included in the repository.

Example: python rdatafile.py testdata.md49 output.png

Will result in the output.png file.
