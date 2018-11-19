import re
import argparse
import numpy as np
import operator
import matplotlib.pyplot as plt

class BondInfo:

    def __init__(self):

        self.id = 0
        self.pos = [0.,0.,0.]

        self.atomid = [0,0]
        self.atomsymb = ['.','.']
        self.atompos = [[0.,0.,0.], [0.,0.,0.]]

        self.area = 0
        self.chg = 0

        self.p = []
        self.val = []
        self.ichg = []
        self.iarea = []

    def title(self):
        return '[{}-{} == {}-{}]'.format(self.atomsymb[0], self.atomid[0], self.atomsymb[1], self.atomid[1])

    def organize(self, l, r):

        if self.atomsymb[0] in l and self.atomsymb[1] in r:
            return

        self.atomsymb[0],self.atomsymb[1] = self.atomsymb[1],self.atomsymb[0]
        self.atomid[0],self.atomid[1] = self.atomid[1],self.atomid[0]
        self.atompos[0],self.atompos[1] = self.atompos[1],self.atompos[0]

        self.p = list(reversed(self.p))
        self.val = list(reversed(self.val))

    def disp(self):
        print 'bond', self.id, 'connects to', self.atomid, ':', self.atomsymb, ':',
        print 'area =', self.area, 'chg =', self.chg,

        l = len(self.p)
        print ' #vals =', l,
        if l > 0:
            print ': [', min(self.p), max(self.p), ']',
            print '[', min(self.val), max(self.val), ']',
            #print '[', min(self.ichg), max(self.ichg), ']',
            #print '[', min(self.iarea), max(self.iarea), ']'
        else:
            print ''

def replace(str, delims):
    for d in delims:
        str = str.replace(d, ' ')
    return str

def read_bondStats(infilename):

    text_file = open(infilename, 'r')
    lines = text_file.read().split('\n')
    text_file.close()

    bonds = []
    s = None
    nlines = len(lines)

    lidx = 0
    while(True):

        if lidx >= nlines:
            break

        l = lines[lidx]
        #print lidx, l, len(l)

        if 'bond_critical_point' in l:

            if s != None:
                s.organize(['Ca'], ['Cl'])
                s.disp()
                bonds.append(s)

            s = BondInfo()

            # read the saddle info from the current line
            tokens = replace(l, '[]()').split(' ')
            s.id = int(tokens[4])
            for i in xrange(3):
                s.pos[i] = float(tokens[8+i])

            # next two lines will give info about atoms
            for a in xrange(2):
                lidx += 1
                l = lines[lidx]
                tokens = replace(l, '[]()').split(' ')
                s.atomid[a] = int(tokens[7])
                s.atomsymb[a] = tokens[9]
                for i in xrange(3):
                    s.atompos[a][i] = float(tokens[11+i])

            # next line will give charge and area
            lidx += 1
            l = lines[lidx]
            tokens = replace(l,',').split(' ')
            s.area = float(tokens[4])
            s.chg = float(tokens[7])

        elif len(l) > 0:
            l = [float(k) for k in l.split(',')]
            s.p.append(l[0])
            s.val.append(l[1])
            #s.ichg.append(l[2])
            #s.iarea.append(l[3])

        lidx += 1

    s.organize(['Ca'], ['Cl'])
    s.disp()
    bonds.append(s)

    #print 'read', len(bonds), 'bonds'
    return bonds


def plot_bonds(bonds, type):

    assert type in [0,1,2]

    labels = ['Ca','Ca','Cl','Cl','Cl','Cl']
    colors = ['black', 'red', 'green', 'blue', 'yellow', 'cyan']
    atoms = []

    #plt.figure(figsize=(8,6))
    plt.figure()
    for s in bonds:

        x = s.p
        if type == 0:
            y = s.val
        elif type == 1:
            y = s.ichg
        elif type == 2:
            y = s.iarea

        plt.semilogy(x, y, linewidth=1, alpha=0.5) #linestyle='--', label=s.title())

        atoms.append([s.atomid[0], x[0],y[0]])
        atoms.append([s.atomid[1], x[-1],y[-1]])

    atoms.sort(key = operator.itemgetter(0))
    atoms = np.array(atoms)
    atoms[:,1] = atoms[:,1].round(decimals=2)
    atoms[:,2] = atoms[:,2].round(decimals=6)
    atoms = np.unique(atoms, axis=0)


    for id in xrange(1,7):
        ids = np.where(atoms[:,0] == id)[0]
        plt.scatter(atoms[ids,1],atoms[ids,2], c=colors[id-1], label=labels[id-1], alpha=0.5)


    plt.legend(loc=3)
    plt.xlim([-2,2])
    plt.xlabel('distance from saddle along bond path (in Angstroms)')

    if type == 0:
        plt.ylabel('charge (number of electrons)')
    elif type == 1:
        plt.ylabel('integrated charge (number of electrons)')
    elif type == 2:
        plt.ylabel('integrated area')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Plot bond statistics output from TopoMS')

    parser.add_argument('--infile', metavar='(infile)', required=True, nargs=1, help='Input stats file')
    args = parser.parse_args()

    bonds = read_bondStats(args.infile[0])

    # --------------------------------------------------------------------------
    # plot all bonds
    outfile = args.infile[0][:-4]

    plot_bonds(bonds, 0)
    plt.savefig(outfile+'_chg.pdf')
    print 'Saved plot:', outfile+'_chg.pdf'

    '''
    plot_bonds(bonds, 1)
    plt.savefig(outfile+'_ichg.pdf')
    print 'Saved plot:', outfile+'_ichg.pdf'

    plot_bonds(bonds, 2)
    plt.savefig(outfile+'_iarea.pdf')
    print 'Saved plot:', outfile+'_iarea.pdf'
    '''
    # --------------------------------------------------------------------------
