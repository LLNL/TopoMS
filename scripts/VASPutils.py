import os
import numpy
import math

# Read VASP CAR format (CHGCAR, AECCAR, LOCPOT)
def readCAR(infilename, read_field = True):

    if not os.path.isfile(infilename):
        raise Exception('File ({}) not found'.format(infilename))

    # --------------------------------------------
    # start reading CAR file
    f = open(infilename, 'r')
    print " --> Reading file ", f

    # --------------------------------------------
    print " --> Reading header"
    sysname = f.readline()[:-1]
    scalingfactor = float(f.readline())

    lattice = []
    lattice.append( [float(x) for x in f.readline().split()] )
    lattice.append( [float(y) for y in f.readline().split()] )
    lattice.append( [float(z) for z in f.readline().split()] )

    species = []
    specieCounts = []

    l = f.readline().split()

    # if this is not a digit, then these must be symbols
    if(l[0].isdigit() == False):
        species = l
        specieCounts = [int(cnt) for cnt in f.readline().split()]

    # else, this must already be the count. create fake symbols
    else:
        specieCounts = [int(cnt) for cnt in l]
        species = ['-' for x in specieCounts]

    numAtoms = sum(specieCounts)
    species = zip(species, specieCounts)

    # --------------------------------------------
    print " --> Reading atom positions"
    coordType = f.readline()[:-1]

    pos = []
    for i in range(0, numAtoms):
        pos.append( [float(x) for x in f.readline().split()] )

    print " --> Reading grid information"
    t = f.readline()
    while len(t.split()) == 0:
        t = f.readline()

    grid = [int(x) for x in t.split()]
    sz = grid[0]*grid[1]*grid[2]

    # --------------------------------------------
    print " --> Reading data "
    data = numpy.zeros(sz, dtype=numpy.float32)
    if(read_field == True):

        lines = f.readlines()
        curr_count = 0
        for line in lines:

            if line.startswith('augmentation'):
                #print curr_count, data[curr_count-1]
                break

            for token in line.split():
                if curr_count < sz-1:
                    data[curr_count] = numpy.float32(token)
                curr_count = curr_count+1

            if(curr_count >= sz):
                break

    f.close()

    # ----------------------------------------------------------
    # ---- 09.13.2016
    # I want to index data as (x,y,z), so reshape in fortran order
    # this makes it easier to convert to vti, etc.

    data = numpy.reshape(data, [grid[0], grid[1], grid[2]], order = 'F')

    print '\t Data.shape = ', data.shape, '(in Fortran order)'
    return (sysname, scalingfactor, lattice, species, pos, grid, data)


def writeCAR(outfilename, sysname, scalingfactor, lattice, species, positions, grid, data, write_data = False):

    # --------------------------------------------
    # start reading CAR file
    f = open(outfilename, 'w')
    print " --> Writing file ", f

    # --------------------------------------------
    print " --> Writing header"
    print >> f, sysname
    print >> f, '%16.12f' % scalingfactor

    for lat in lattice:
        print >> f, '%12.6f %12.6f %12.6f' % (lat[0], lat[1], lat[2])

    for sp in species:
        f.write('  '+str(sp[0]))
    f.write('\n')
    for sp in species:
        f.write('  '+str(sp[1]))
    f.write('\n')

    print " --> Writing positions"
    f.write('Direct\n')

    for pos in positions:
        print >> f, '%10.6f %10.6f %10.6f' % (pos[0], pos[1], pos[2])

    f.write('\n'+str(grid[0])+' '+str(grid[1])+' '+str(grid[2])+'\n')

    if (write_data):

        print " --> Writing data"
        data = data.flatten('F')
        #data = numpy.reshape(data, grid[0] * grid[1] * grid[2])

        # write 5 values in a line
        nlines = data.shape[0] / 5
        for i in range(0, nlines):
            #print >> f, '%10.6f %10.6f %10.6f %10.6f %10.6f' % (data[5*i], data[5*i+1], data[5*i+2], data[5*i+3], data[5*i+4])
            print >> f, '%0.11E %0.11E %0.11E %0.11E %0.11E' % (data[5*i], data[5*i+1], data[5*i+2], data[5*i+3], data[5*i+4])

        data = numpy.reshape(data, [grid[0], grid[1], grid[2]], order = 'F')

    f.close()

def readPOSCAR(infilename):

    # --------------------------------------------
    # start reading CAR file
    f = open(infilename, 'r')
    print " --> Reading file ", f

    # --------------------------------------------
    print " --> Reading header"
    sysname = f.readline()[:-1]
    scalingfactor = float(f.readline())

    lattice = []
    lattice.append( [float(x) for x in f.readline().split()] )
    lattice.append( [float(y) for y in f.readline().split()] )
    lattice.append( [float(z) for z in f.readline().split()] )

    species = []
    specieCounts = []

    l = f.readline().split()

    # if this is not a digit, then these must be symbols
    if(l[0].isdigit() == False):
        species = l
        specieCounts = [int(cnt) for cnt in f.readline().split()]

    # else, this must already be the count. create fake symbols
    else:
        specieCounts = [int(cnt) for cnt in l]
        species = ['-' for x in specieCounts]

    numAtoms = sum(specieCounts)
    species = zip(species, specieCounts)

    # --------------------------------------------
    print " --> Reading atom positions"
    coordType = f.readline()[:-1]

    pos = []
    for i in range(0, numAtoms):
        pos.append( [float(x) for x in f.readline().split()] )

    # --------------------------------------------
    print " --> Reading atom velocties"

    vel = []
    for i in range(0, numAtoms):
        l = f.readline()
        if len(l) == 0:
            break;
        vel.append( [float(x) for x in l.split()] )

    return (sysname, scalingfactor, lattice, species, pos, vel)

def writePOSCAR(outfilename, sysname, scalingfactor, lattice, species, positions, velocities):

    # --------------------------------------------
    # start reading CAR file
    f = open(outfilename, 'w')
    print " --> Writing file ", f

    # --------------------------------------------
    print " --> Writing header"
    print >> f, sysname
    print >> f, '%16.12f' % scalingfactor

    for lat in lattice:
        print >> f, '%12.6f %12.6f %12.6f' % (lat[0], lat[1], lat[2])

    for sp in species:
        f.write('  '+str(sp[0]))
    f.write('\n')
    for sp in species:
        f.write('  '+str(sp[1]))
    f.write('\n')

    print " --> Writing positions"
    f.write('Direct\n')

    for pos in positions:
        print >> f, '%11.7f %11.7f %11.7f' % (pos[0], pos[1], pos[2])

    print >> f, ''
    for vel in velocities:
        print >> f, '%11.7f %11.7f %11.7f' % (vel[0], vel[1], vel[2])

# Read VASP XDATCAR
def readXDATCAR(infilename, read_only = None):

    # --------------------------------------------
    # start reading CAR file
    f = open(infilename, 'r')
    print " --> Reading file ", f

    # --------------------------------------------
    print " --> Reading header"
    t = f.readline().split()

    numatoms = int(t[0])
    tmp = int(t[1])
    numts = int(t[2])

    t = f.readline().split()

    avol = float(t[0])
    lattice = [float(t[1]), float(t[2]), float(t[3])]
    dt = float(t[4])

    temp = float(f.readline())

    # --- correct the units
    dt = dt*math.pow(10,12)
    lattice = [ k*math.pow(10,10) for k in lattice ]

    f.readline()
    f.readline()
    # --------------------------------------------

    if(read_only != None):
        numts = read_only

    print " --> Allocating space for", numts, "timesteps"
    P = numpy.ndarray( (numatoms, numts, 3), dtype = numpy.float32 )

    print " --> Reading positions"

    # for each time-step and each atom
    for t in range(numts):

        f.readline()  # direct!
        for a in range(numatoms):
            P[a, t] = [float(x) for x in f.readline().split()]

    f.close()
    return (lattice, dt, P)
