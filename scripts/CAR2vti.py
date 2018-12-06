# Read VASP CAR format (CHGCAR, AECCAR, LOCPOT) and write the field as a vti file!
import numpy
import argparse

from VASPutils import readCAR
from pyevtk.hl import imageToVTK

# ---------------------------------------------
parser = argparse.ArgumentParser(description='Read a CAR file and write as vti')
parser.add_argument('--infile', metavar='(infile)', required=True, nargs=1, help='Input CAR file')

args = parser.parse_args()
infilename = args.infile[0]

# ---------------------------------------------
(sysname, scalingfactor, lattice, species, pos, grid, data) = readCAR(infilename, True)
spacings = (scalingfactor*lattice[0][0]/(grid[0]-1),
            scalingfactor*lattice[1][1]/(grid[1]-1),
            scalingfactor*lattice[2][2]/(grid[2]-1))

# ---------------------------------------------
# now, write the data
print 'Writing', infilename+'.vti'
imageToVTK(infilename, spacing = spacings, pointData = {"charge" : data} )
