# rdatafile
Simple minimalistic offline plotting of CADI datafiles

this code is written for python >= 3.6 and require matplotlib and numpy
usage: rdatafile.py inputfile outputfile
inputfile = file name of the CADI data input file
outputfile = file name of the resulting plot

outputfile should have extension .png to produce png-files or .pdf to produce pdf-files.
supported formats eps, pdf, pgf, png, ps, raw, rgba, svg, svgz (depends on version of matplotlib)

the code is written based on IDL code provided by Chris Meek
IDL code originally written by Ian Grant and modified by the same and JWM
