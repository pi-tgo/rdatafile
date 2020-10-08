# rdatafile
Simple minimalistic offline plotting of CADI datafiles

The code is written for python3 and require matplotlib and numpy.

Usage: rdatafile.py inputfile outputfile
Usage: rdatafile-mdx.py inputfile prefix

inputfile = file name of the CADI data input file
outputfile = file name of the resulting plot
prefix = prefix to the output files (when plotting several ionograms from a single inputfile)

Outputfile should have extension .png to produce png-files or .pdf to produce pdf-files.
Supported formats eps, pdf, pgf, png, ps, raw, rgba, svg, svgz (depends on version of matplotlib)

The code is written based on IDL code obtained from the Canadian High Arctic Ionospheric Network (CHAIN) web pages; http://chain.physics.unb.ca/chain/.
IDL code originally written by Ian Grant and modified by the same and JWM"

Example datafile (testdata.md4) is included in the repository. This file contain a single ionogram.
Another example datafile (testdata.md2) is also included which contain 60 ionograms (only noise).

Example: **python rdatafile.py testdata.md4 output.png**

will result in a output.png file.

Example: **python rdatafile-mdx.py testdata.md2 prefix**

will result in multiple prefixX.png files. The number of ionograms is defined in the code
(default 5).
