import argparse

class BondInfo:

    def __init__(self):
        self.ids = []
        self.area = 0
        self.chg = 0
        self.p = []
        self.y = []

    def disp(self):
        print 'bond', self.ids[0], 'connects to', self.ids[1:], ':',
        print 'area =', self.area, 'chg =', self.chg,

        l = len(self.p)
        print '#vals =', l,
        if l > 0:
            print ': [', min(self.p), max(self.p), ']',
            print '[', min(self.y), max(self.y), ']'
        else:
            print ''

def read_bondStats(infilename):

    text_file = open(infilename, 'r')
    lines = text_file.read().split('\n')
    text_file.close()

    bonds = []
    s = None
    for l in lines:

        if len(l) == 0:
            continue

        elif 'saddle' == l[:6]:
            if s != None:
                s.disp()
                bonds.append(s)

            s = BondInfo()
            vals = [v.split() for v in l.split(',')]
            s.ids.append(int(vals[0][1]))
            s.ids.append(int(vals[1][1]))
            s.ids.append(int(vals[1][2]))
            s.area = float(vals[2][1])
            s.chg = float(vals[3][1])

        else:
            l = [float(k) for k in l.split(',')]
            s.p.append(l[0])
            s.y.append(l[1])

    print 'read', len(bonds), 'bonds'
    return bonds


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Plot bond statistics output from TopoMS')

    parser.add_argument('--infile', metavar='(infile)', required=True, nargs=1, help='Input stats file')
    args = parser.parse_args()

    bonds = read_bondStats(args.infile[0])


    # plot all bonds
    import matplotlib.pyplot as plt
    colors = ['black', 'red', 'green', 'blue', 'yellow', 'cyan']

    for s in bonds:

        plt.semilogy(s.p, s.y,linewidth=1.5)
        #nx = len([p for p in s.p if p < 0])

        #plt.semilogy(s.p[:nx], s.y[:nx],linewidth=2)#,color=colors[s.ids[1]-1])
        #plt.semilogy(s.p[nx:], s.y[nx:],linewidth=2)#,color=colors[s.ids[2]-1])

    plt.title(args.infile[0])
    plt.ylabel('charge')
    plt.xlabel('distance from saddle along bond path')
    #plt.show()
    plt.savefig(args.infile[0]+'.png')
