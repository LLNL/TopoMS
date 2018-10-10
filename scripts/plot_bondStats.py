import argparse
import re
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
        self.y = []

    def title(self):
        return '[{}-{} == {}-{}]'.format(self.atomsymb[0], self.atomid[0], self.atomsymb[1], self.atomid[1])

    def organize(self, l, r):

        if self.atomsymb[0] in l and self.atomsymb[1] in r:
            return

        self.atomsymb[0],self.atomsymb[1] = self.atomsymb[1],self.atomsymb[0]
        self.atomid[0],self.atomid[1] = self.atomid[1],self.atomid[0]
        self.atompos[0],self.atompos[1] = self.atompos[1],self.atompos[0]

        self.p = list(reversed(self.p))
        self.y = list(reversed(self.y))

    def disp(self):
        print 'bond', self.id, 'connects to', self.atomid, ':', self.atomsymb, ':',
        print 'area =', self.area, 'chg =', self.chg,

        l = len(self.p)
        print ' #vals =', l,
        if l > 0:
            print ': [', min(self.p), max(self.p), ']',
            print '[', min(self.y), max(self.y), ']'
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
            s.y.append(l[1])

        lidx += 1

    s.organize(['Ca'], ['Cl'])
    s.disp()
    bonds.append(s)


    #print 'read', len(bonds), 'bonds'
    return bonds


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Plot bond statistics output from TopoMS')

    parser.add_argument('--infile', metavar='(infile)', required=True, nargs=1, help='Input stats file')
    args = parser.parse_args()

    bonds = read_bondStats(args.infile[0])

    # --------------------------------------------------------------------------
    # plot all bonds
    import matplotlib.pyplot as plt
    colors = ['black', 'red', 'green', 'blue', 'yellow', 'cyan']

    plt.figure(figsize=(8,6))
    for s in bonds:

        plt.semilogy(s.p, s.y,linewidth=1, linestyle='--', label=s.title())

        #nx = len([p for p in s.p if p < 0])
        #plt.semilogy(s.p[:nx], s.y[:nx],linewidth=2)#,color=colors[s.ids[1]-1])
        #plt.semilogy(s.p[nx:], s.y[nx:],linewidth=2)#,color=colors[s.ids[2]-1])

    plt.title(args.infile[0])
    plt.legend()
    plt.ylabel('charge')
    plt.xlabel('distance from saddle along bond path')
    #plt.show()

    outfile = args.infile[0]+'.pdf'
    plt.savefig(outfile)
    print 'Saved plot:', outfile
