import os
import sys
import re
import numpy as np
import math
from itertools import repeat
from paraparser.indexparser import parse_index, rev_parse_index
from fileparser.parse49 import parse_49
from filegenerator.genfchk import quicksave

VERSION = 'V4.0'
VERSIONTEXT = 'Original published version honorably delivered by Acid&Francium (with a blast).\nRevised on 16 Jan 2020 by Acid.\n'
REFERENCE_LIST = '''
'''

def myeigh(mat, rev=False):
    evals, evecs = np.linalg.eigh(mat)
    evecs = evecs.T
    if rev:
        order = evals.argsort()[::-1]
    else:
        order = evals.argsort()
    evals = evals[order]
    evecs = evecs[order]
    return evals, evecs

def PIO(ffchk, fn49='', fragmentation=None, silent=False):
    path, fn = os.path.split(ffchk)
    title, ext = os.path.splitext(fn)

    if fn49 == '':
        for fn in [os.path.join(path, title + '.49'), 
                os.path.join(path, title.upper() + '.49'), 
                os.path.join(path, title.upper() + 'FILE'[len(title)-4:] + '.49')]:
            if os.path.isfile(fn):
                fn49 = fn
                break
        else:
            raise UserWarning('.49 file not found.')

    raw = parse_49(fn49)
    naoao = raw['NAOAO']
    naolabels = raw['NAOAOlabel']
    if 'DNAO' in raw:
        dim_spin = 1
        dmnao = np.array([raw['DNAO']])
        if 'FNAO' in raw:
            fmnao = np.array([raw['FNAO']])
        else:
            fmnao is None
    else:
        assert 'DNAO1' in raw
        assert 'DNAO2' in raw
        dim_spin = 2
        dmnao = np.array([raw['DNAO1'], raw['DNAO2']])
        if 'FNAO1' in raw and 'FNAO2' in raw:
            fmnao = np.array([raw['FNAO1'], raw['FNAO2']])
        else:
            fmnao is None
    dim_nao, dim_bs = naoao.shape
    assert dmnao.shape == (dim_spin, dim_nao, dim_nao)
    if not silent:
        print 'Density matrix read in.'

    if fragmentation is None:
        prompt = '$ Please input the atom ID of two fragments: (e.g. 1-5,8,13 6-7,9-12)\n$ '
        fragmentation = raw_input(prompt)

    fragments = fragmentation.strip().lower().split()
    assert len(fragments) == 2
    atoms1 = parse_index(fragments[0])
    atoms2 = parse_index(fragments[1])
    oids1 = np.array([orb for orb in range(dim_nao) if naolabels[orb].center in atoms1])
    oids2 = np.array([orb for orb in range(dim_nao) if naolabels[orb].center in atoms2])
    assert len(set(oids1)) == len(oids1)
    assert len(set(oids2)) == len(oids2)
    assert not set(oids1).intersection(set(oids2))
    if set(oids1).union(set(oids2)) != set(range(dim_nao)):
        print 'Warning: Incomplete fragmentation detected!'
    if set([naolabels[orb].center for orb in oids1]) != set(atoms1):
        print 'Warning: Ghost atoms detected in fragment A!'
    if set([naolabels[orb].center for orb in oids2]) != set(atoms2):
        print 'Warning: Ghost atoms detected in fragment B!'
 
    data = {}
    data['srcfile'] = fn49
    data['title'] = title
    data['dim'] = (dim_spin, dim_nao, dim_bs)
    data['fragments'] = (sorted(set([naolabels[orb].center for orb in oids1])), 
        sorted(set([naolabels[orb].center for orb in oids2])))
    data['oids'] = (oids1, oids2)
    data['NAOAO'] = naoao
    data['DNAO'] = dmnao
    if fmnao is not None:
        data['FNAO'] = fmnao
    data['NAOAOlabel'] = naolabels
    
    data.update(genPIO(dmnao, oids1, oids2, fmnao=fmnao, silent=silent))

    npio = data['#PIO']
    pbi = data['PBI']
    pbio = np.zeros((dim_spin, dim_nao))
    for spin in xrange(dim_spin):
        pbio[spin,oids1[:npio]] = pbi[spin,:npio]
        pbio[spin,oids2[:npio]] = pbi[spin,:npio]
    data['PBIO'] = pbio

    save(data, ffchk, silent=silent)
    return data

def genPIO(dmnao, oids1, oids2, fmnao=None, tol=1e-6, silent=False):
    dim_spin = dmnao.shape[0]
    dim = dmnao.shape[-1]
    assert dmnao.shape == (dim_spin, dim, dim)
    oids1, oids2 = np.array(oids1), np.array(oids2)

    npio = min(len(oids1), len(oids2))
    pionao = np.zeros((dim_spin, dim, dim))
    pbi = np.zeros((dim_spin, npio))
    pimonao = np.zeros((dim_spin, dim, dim))
    for spin, dmnao_spin in enumerate(dmnao):
        vals1, vecs1 = myeigh(dmnao_spin[oids1[:,None],oids1])
        vals2, vecs2 = myeigh(dmnao_spin[oids2[:,None],oids2])
        off_block = vecs1.dot(dmnao_spin[oids1[:,None],oids2]).dot(vecs2.T)
        off_block[abs(off_block) < tol] = 0
        U, D, V = np.linalg.svd(off_block, full_matrices=True)
        U = vecs1.T.dot(U)
        V = V.dot(vecs2)
        D2 = D**2
        pionao[spin,oids1[:,None],oids1] = U.T
        pionao[spin,oids2[:,None],oids2] = V
        pbi[spin,:npio] = D2 * dim_spin

        dmpio = pionao[spin].dot(dmnao[spin]).dot(pionao[spin].T)
        nep = int(math.ceil(np.round(dmnao[spin].trace()) * dim_spin / 2))
        for nullspace in (np.array(oids1)[min(nep,npio):], np.array(oids2)[min(nep,npio):]):
            evals, evecs = myeigh(dmpio[nullspace[:,None],nullspace], rev=True)
            pionao[spin,nullspace] = evecs.dot(pionao[spin,nullspace])
            try:
                fm = pionao[spin].dot(fmnao[spin]).dot(pionao[spin].T)
                null2 = nullspace[abs(evals-2)<tol]
                null0 = nullspace[abs(evals-0)<tol]
                evals, evecs = myeigh(fm[null2[:,None],null2])
                pionao[spin,null2] = evecs.dot(pionao[spin,null2])
                evals, evecs = myeigh(fm[null0[:,None],null0])
                pionao[spin,null0] = evecs.dot(pionao[spin,null0])
            except NameError:
                pass

        pimonao[spin] = pionao[spin]
        for oid1, oid2 in zip(oids1, oids2)[:nep]:
            D = dmpio[[oid1, oid2]][:,[oid1,oid2]]
            evals, evecs = myeigh(D,rev=True)
            # assert np.allclose(evals, [2/dim_spin,0], atol=tol)
            evecs *= np.sign(evecs.sum(axis=1))[:,None]
            pimonao[spin,[oid1,oid2]] = evecs.dot(pionao[spin,[oid1,oid2]])
 
    if not silent:
        print 'PIO analysis completed.'
    data = {'#PIO':npio, 'PIONAO': pionao, 'PIMONAO': pimonao, 'PBI': pbi}
    return data

def save(data, ffchk, silent=False):
    savetxt(data, os.path.splitext(ffchk)[0]+'_pio.txt', silent=silent)
    savefchk(data, ffchk, silent=silent)
    # saveawc(data, os.path.splitext(ffchk)[0]+'_pio.awc', silent=silent) # Internal use only for now

def savefchk(data, ffchk, silent=False):
    naoao = data['NAOAO']
    pionao = data['PIONAO']
    pimonao = data['PIMONAO']
    pbio = data['PBIO']

    pioao = pionao.dot(naoao)
    pimoao = pimonao.dot(naoao)

    quicksave(ffchk, pioao, pbio, suffix='_pio', overwrite=True)
    quicksave(ffchk, pimoao, pbio, suffix='_pimo', overwrite=True)
    if not silent:
        print 'PIO/PIMO saved to fchk files.'

def savetxt(data, ofn, silent=False):
    atoms1, atoms2 = data['fragments']
    dim_spin, dim_pio, dim_bs = data['dim']
    oids1, oids2 = data['oids']
    naoao = data['NAOAO']
    pbi = data['PBI']
    pionao = data['PIONAO']
    pimonao = data['PIMONAO']
    dmnao = data['DNAO']
    fmnao = data.get('FNAO', np.zeros_like(dmnao))
    npio = data['#PIO']

    dmpio = np.array([pionao[_].dot(dmnao[_]).dot(pionao[_].T) for _ in xrange(dim_spin)])
    fmpio = np.array([pionao[_].dot(fmnao[_]).dot(pionao[_].T) for _ in xrange(dim_spin)])
    dmpimo = np.array([pimonao[_].dot(dmnao[_]).dot(pimonao[_].T) for _ in xrange(dim_spin)])
    fmpimo = np.array([pimonao[_].dot(fmnao[_]).dot(pimonao[_].T) for _ in xrange(dim_spin)])
    total = pbi.sum()

    mod = '%5d%1s %10.5f %10.5f %5d%1s %10.5f %10.5f %10.5f %10.5f %8.2f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f'
    modt = re.sub('(\.\d+)?[fd]', 's', mod)
    with open(ofn, 'w') as f:
        f.write('PIO %s\n' % VERSION)
        f.write('%s\n' % VERSIONTEXT)
        f.write('Fragment A: %s (Orbitals: %s)\n' % (rev_parse_index(atoms1), rev_parse_index(oids1+1)))
        f.write('Fragment B: %s (Orbitals: %s)\n' % (rev_parse_index(atoms2), rev_parse_index(oids2+1)))
        if dim_spin == 1:
            f.write('Total interaction: %f\n' % total)
        else:
            f.write('Total interaction: %f\n' % total)
            f.write('Alpha interaction: %f\n' % pbi[0].sum())
            f.write('Beta  interaction: %f\n' % pbi[1].sum())
        f.write('\n')

        for spin in xrange(dim_spin):
            nep = int(math.ceil(np.round(dmnao[spin].trace()) * dim_spin / 2))
            pbispin = pbi[spin]
            pop1 = np.diag(dmpio[spin])[oids1]
            pop2 = np.diag(dmpio[spin])[oids2]
            energy1 = np.diag(fmpio[spin])[oids1]
            energy2 = np.diag(fmpio[spin])[oids2]
            popb = np.diag(dmpimo[spin])[oids1]
            popa = np.diag(dmpimo[spin])[oids2]
            energyb = np.diag(fmpimo[spin])[oids1]
            energya = np.diag(fmpimo[spin])[oids2]
            contrib = pbispin / total * 100
            cum = np.cumsum(contrib)
            ie = [fmpio[spin][i][j] for i, j in zip(oids1, oids2)]
            rie = [fmpimo[spin][i][j] for i, j in zip(oids1, oids2)]

            if dim_spin == 2:
                spincode = 'ALPHA' if spin == 0 else 'BETA'
                f.write('---%s---\n' % spincode)
            f.write('Fragment A'.center(26))
            f.write('Fragment B'.center(26))
            f.write('\n')
            f.write(modt % ('Or', 'b', 'Pop   ', 'E    ', 'Or', 'b', 'Pop   ', 'E    ', 
                'Ixn   ', 'Contrib%', 'Cum%', 'IE   ', 
                'EB   ', 'EA   ', 'PopB  ', 'PopA  ', 'RIE   '))
            f.write('\n')

            if dim_spin == 2:
                orbspin = 'a' if spin == 0 else 'b'
            else:
                orbspin = ''
            text = zip(oids1 + 1, repeat(orbspin), pop1, energy1, 
                    oids2 + 1, repeat(orbspin), pop2, energy2, 
                    pbispin[:min(nep,npio)], contrib, cum, ie, 
                    energyb, energya, popb, popa, rie)
            while text:
                f.write(mod % text.pop(0))
                f.write('\n')
            f.write('\n\n')

        f.close()
    if not silent:
        print 'PIO result saved to file.'

def saveawc(data, ofn, silent=False):
    dim_spin, dim_pio, dim_bas = data['dim']
    naolabels = data['NAOAOlabel']
    pionao = data['PIONAO']
    pbio = data['PBIO']

    assert pbio.shape == (dim_spin, dim_pio)
    assert len(naolabels) == dim_pio
    assert pionao.shape == (dim_spin, dim_pio, dim_pio)

    atoms = sorted(set([nao.center for nao in naolabels]))
    awc = np.zeros((dim_spin, dim_pio, len(atoms)))
    for spin in xrange(dim_spin):
        for pio in xrange(dim_pio):
            for nao in xrange(dim_pio):
                awc[spin, pio, naolabels[nao].center-1] += pionao[spin,pio,nao] ** 2

    with open(ofn, 'w') as f:
        for spin in xrange(dim_spin):
            if dim_spin == 2:
                spincode = 'ALPHA' if spin == 0 else 'BETA'
                f.write('---%s---\n' % spincode)
            for pio in xrange(dim_pio):
                if pbio[spin,pio] >= 0.1 / dim_spin:
                    f.write('PIO %-4d\n' % (pio + 1))
                    for at in xrange(len(atoms)):
                        if awc[spin, pio, at] > 0.1:
                            f.write('        %-4d         %8.4f\n' % (at + 1, awc[spin, pio, at]))
                        for nao in xrange(dim_pio):
                            if naolabels[nao].center == at + 1 and pionao[spin,pio,nao] ** 2 >= 0.1:
                                f.write('            %-8s %8.4f\n' % (naolabels[nao].label(), pionao[spin,pio,nao] ** 2))
            f.write('\n\n')
        f.close()
    if not silent:
        print 'Atom/NAO contributions saved to file.'

def main():
    ffchk = sys.argv[1]
    fn49 = sys.argv[2] if sys.argv[2:] else ''
    fragmentation = sys.argv[3:5] if sys.argv[3:] else None
    PIO(ffchk, fn49, fragmentation)

if __name__ == '__main__':
    main()


