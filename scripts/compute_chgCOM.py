import os
import sys
import numpy as np
import argparse

from VASPutils import readCAR
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------
# VTK I/O
# ------------------------------------------------------------------------------
def read_vti(filename, fieldname):

    if not os.path.isfile(filename):
        raise Exception('File ({}) not found'.format(filename))

    import vtk
    from vtk.util import numpy_support

    print ' --> Reading ({}) from ({})'.format(fieldname, filename)

    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(filename)
    reader.Update()
    dims = reader.GetOutput().GetDimensions()
    field = reader.GetOutput().GetPointData().GetArray(fieldname)
    field = numpy_support.vtk_to_numpy(field).reshape(dims[2], dims[1], dims[0])
    field = np.swapaxes(field, 0,2)
    return field

def write_vtp(filename, data):

    import vtk
    from vtk.util import numpy_support

    print ' --> Writing ({})'.format(filename)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(filename)

    polydata = vtk.vtkPolyData()

    points = vtk.vtkPoints()
    cells = vtk.vtkCellArray()
    labels = []

    # --------------------------------------------------------------------------
    pid = 0
    kid = 0
    for key in data.keys():

        pdata = data[key]
        #print key, '-->', kid, len(pdata)

        for p in pdata:
            points.InsertNextPoint(p[0],p[1],p[2])
            pid += 1

            cell = vtk.vtkVertex()
            cell.GetPointIds().SetId(0, pid-1)
            cells.InsertNextCell(cell)

            labels.append(kid)
        kid += 1

    polydata.SetPoints(points)
    polydata.SetPolys(cells)

    VTK_data = numpy_support.numpy_to_vtk(num_array=np.array(labels))
    VTK_data.SetName('plabels')
    polydata.GetPointData().AddArray(VTK_data)

    VTK_data2 = numpy_support.numpy_to_vtk(num_array=np.array(labels))
    VTK_data2.SetName('clabels')
    polydata.GetCellData().AddArray(VTK_data2)

    writer.SetInputData(polydata)
    writer.Write()

# ------------------------------------------------------------------------------
# text output of cocg
# ------------------------------------------------------------------------------
def write_chgcom(filename, apos, cpos):

    print 'Writing', filename, '...',
    sys.stdout.flush()
    fp = open(filename, 'w')

    fp.write(" %4s %11s %11s %11s %13s %13s %13s %13s\n" %
                ("#", "X", "Y", "Z", "cX", "cY", "cZ", "magn(c)"));
    fp.write("-------------------------------------------------------------------------------------------------\n")

    for k in xrange(natoms):
        a = apos[k]
        c = cpos[k]
        fp.write(" %4s %11.6f %11.6f %11.6f %+1.6e %+1.6e %+1.6e %+1.6e\n" %
                    (k, a[0], a[1], a[2], c[0], c[1], c[2], np.linalg.norm(np.array(c))));

        #0:+1.2f

    fp.write("-------------------------------------------------------------------------------------------------\n")
    fp.close()
    print 'Done!'

# ------------------------------------------------------------------------------
def periodic_displacement(a, b, dims):

    debug = False
    disp = a-b
    hdims = 0.5*dims
    if debug:
        print 'compute pdisp', a, b, dims, hdims
        print disp
    for i in [0,1,2]:
        if disp[i] > hdims[i]:
            if debug:
                print 'case 1: ', disp[i],
            disp[i] = disp[i] - dims[i]
            if debug:
                print '--->', disp[i],

        elif disp[i] < -hdims[i]:
            if debug:
                print 'case 2: ', disp[i],
            disp[i] = disp[i] + dims[i]
            if debug:
                print '--->', disp[i],
            #exit()
    if debug:
        print disp
    return disp

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if __name__ == '__main__':

    # --------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description='Compute charge center of mass using Bader analysis done through TopoMS')
    parser.add_argument('--infile', metavar='(infile)', required=True, nargs=1, help='Input CHGCAR file')
    args = parser.parse_args()

    chgVASPfilename = args.infile[0]            # input VASP file
    chgVTIfilename  = chgVASPfilename+'.vti'    # this contains the VASP data in raw format
    volfilename = chgVASPfilename+'-Atoms2Vol.vti'
    periodic = True

    # --------------------------------------------------------------------------
    # read chgcar file to get data about atomic position
    # we read function values from vti, hence pass False at the end
    (sysname, scalingfactor, lattice, species, atompos, grid, data) = readCAR(chgVASPfilename, False)

    def direct2grid(p):
        return np.array([int(round(grid[d]*p[d])) for d in xrange(3)])

    def grid2direct(p):
        return np.array([p[d]/float(grid[d]) for d in xrange(3)])

    def direct2phys(p):
        return np.array([scalingfactor*lattice[d][d]*p[d] for d in xrange(3)])

    def phys2direct(p):
        return np.array([p[d]/(scalingfactor*lattice[d][d]) for d in xrange(3)])

    def grid2phys(p):
        return np.array([scalingfactor*lattice[d][d]*p[d]/float(grid[d]) for d in xrange(3)])

    def phys2grid(p):
        return np.array([p[d]*float(grid[d])/(scalingfactor*lattice[d][d]) for d in xrange(3)])

    grid = np.array(grid)
    pdims = direct2phys([1.0,1.0,1.0])
    spacings = grid2phys([1,1,1])

    voxvol = spacings.prod()
    chgDens_fileUnit2e = 1.0 / float(grid[0]*grid[1]*grid[2]);

    # atoms in CHGCAR files are given in direct coordinates ([0,1])
    # convert them to physical coordinates
    natoms = len(atompos)
    atompos = np.array([direct2phys(apos) for apos in atompos])

    # volumes and charge
    chg = read_vti(chgVTIfilename, 'charge') * chgDens_fileUnit2e   # convert from density to charge!
    vol = read_vti(volfilename, 'atom_labeling')

    # print info about the data
    print '\tLattice  =', pdims
    print '\tGrid     =', grid
    print '\tSpacings =', spacings
    print '\tVox vol  =', voxvol
    print '\tSpecies  =', species
    print '\tSharge   =', chg.shape, chg.dtype, chg.min(), chg.max()
    print '\tLabels   =', vol.shape, vol.dtype, np.unique(vol)
    print '\tAtom pos =\n', atompos

    assert(natoms == len(np.unique(vol)))
    assert(chg.shape[0] == grid[0] and chg.shape[1] == grid[1] and chg.shape[2] == grid[2])
    assert(vol.shape[0] == grid[0] and vol.shape[1] == grid[1] and vol.shape[2] == grid[2])

    # --------------------------------------------------------------------------
    # for each Bader region (each atom)

    # total number of pixels
    anpixels = np.zeros(natoms, np.int)

    # total charge
    ascharges = np.zeros(natoms)

    # net displacement of all pixels (with respect to the atom)
    asdisps   = np.zeros([natoms, 3])

    # center of charge gravity
        # weighted sum of positions of all pixels (with respect to the atom)
    acocg    = np.zeros([natoms, 3])

    # go over all pixels
    for z in xrange(grid[2]):
        for y in xrange(grid[1]):
            for x in xrange(grid[0]):

                l = vol[x,y,z]-1
                c = chg[x,y,z]

                apos = atompos[l]
                ppix = grid2phys([x,y,z])

                # i cant simply add pixels, because of periodic boundary
                # so, compute the displacement with respect to the atom
                disp = periodic_displacement(ppix, apos, pdims)

                anpixels[l]  += 1
                ascharges[l] += c
                asdisps[l]   += disp
                acocg[l]     += c*disp

    '''
    print '   bcnts         :\n', anpixels
    print '   ascharges     :\n', ascharges
    print '   asdisps       :\n', asdisps
    print '   coc-gravity   :\n', acocg
    '''

    anpixels = anpixels.astype(np.float32)
    for k in xrange(natoms):
        acocg[k]   /= anpixels[k]
        asdisps[k] /= anpixels[k]
        anpixels[k] *= voxvol       # convert to physical volume

    for i in xrange(natoms):
        print 'Atom {} at {}: chg = {}, vol = {}, avg_disp {}, cocg = {}'.format(
                        i+1, atompos[i], ascharges[i], anpixels[i], asdisps[i], acocg[i])

    fname = volfilename[:-4] + '-chgcog.txt'
    write_chgcom(fname, atompos, acocg)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
