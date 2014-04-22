

def write_result(R, fname):
    fout=open(fname, 'w')
    label = '\t'.join(map(str, R[R.keys()[0]].keys()))
    fout.write('\t'.join(['0.00', label]) + '\n')
    for p in R:
        value = '\t'.join(map(str, R[p].values()))
        fout.write('\t'.join([str(p), value]) + '\n')
    fout.close()

def average(T):
    TT = {}
    for p in T:
        TT[p] = {}
        for c in T[p]:
            TT[p][c] = sum(T[p][c])/float(len(T[p][c]))
    return TT

def variance(T):
    TT = {}
    for p in T:
        TT[p] = {}
        for c in T[p]:
            mean = sum(T[p][c])/float(len(T[p][c]))
            TT[p][c] = sum([(x-mean)**2 for x in T[p][c]])

    return TT

T = {'s':{}, 'm':{}, 'd':{}, 'it':{}}

for line in open('result.dat'):
    n, c, p, s, m, d, it = map(float, line.rstrip().split('\t'))
    T['s'].setdefault(p, {})
    T['s'][p].setdefault(c, [])
    T['s'][p][c].append(s)
    T['m'].setdefault(p, {})
    T['m'][p].setdefault(c, [])
    T['m'][p][c].append(m)
    T['d'].setdefault(p, {})
    T['d'][p].setdefault(c, [])
    T['d'][p][c].append(d)
    T['it'].setdefault(p, {})
    T['it'][p].setdefault(c, [])
    T['it'][p][c].append(it)

S = average(T['s'])
SV = variance(T['s'])
M = average(T['m'])
D = average(T['d'])
IT = average(T['it'])

write_result(S, 'gcsize.dat')
write_result(SV, 'gcvari.dat')
write_result(M, 'monomer.dat')
write_result(D, 'dimer.dat')
write_result(IT, 'iter.dat')
