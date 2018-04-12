import os
import re
import numpy as np
from fileparser.parsefchk import parse_fchk

LAST_UPDATE = '201706021601'

def gen_fchk(mfn, data, ofn='', title=None, suffix='_rlo', overwrite=False):
    # text: info, energy, coeff, rest
    # data: energy, coeff, dim
    
    if not ofn:
        t, ext = os.path.splitext(mfn)
        ofn = t + suffix + ext
    if not overwrite and os.path.exists(ofn):
        overwriteflag = raw_input('Overwrite file %s? Y/N ' % ofn)
        if overwriteflag.lower() != 'y':
            print ofn, 'skipped.'
            return False
    title = title or os.path.splitext(ofn)[0]
    
    fchk = parse_fchk(mfn)
    text = fchk['text']
    dim_mo, dim_bs = fchk['data']['dim']
    
    dim_nmo, dim_nbs = data.get('dim', (dim_mo, dim_bs))
    assert dim_bs == dim_nbs
    if dim_nmo == dim_mo:
        for e in data['energy']:
            assert e.shape == (dim_nmo,)
        for c in data['coeff']:
            assert c.shape == (dim_nmo, dim_bs)
    elif dim_nmo < dim_mo:
        for spin, e in enumerate(data['energy']):
            data['energy'][spin] = np.append(e, np.zeros(dim_mo - dim_nmo), axis = 0)
            assert data['energy'][spin].shape == (dim_mo,)
        for spin, c in enumerate(data['coeff']):
            data['coeff'][spin] = np.append(c, np.zeros((dim_mo - dim_nmo, dim_bs)), axis = 0)
            assert data['coeff'][spin].shape == (dim_mo, dim_bs)
        dim_nmo = dim_mo
    else:
        # raise AssertionError('Number of new orbitals exceeds number of original orbitals')
        for lid, line in enumerate(text['info']):
            text['info'][lid] = line.replace(str(dim_mo), str(dim_nmo))
        text['energy'][0] = text['energy'][0].replace(str(dim_mo), str(dim_nmo))
        text['coeff'][0] = text['coeff'][0].replace(str(dim_mo*dim_bs), str(dim_nmo*dim_bs))
    if len(data['energy']) != len(text['energy']):
        print 'Warning: Spin symmetry not match!'
    amoeline = text['energy'][0]
    bmoeline = amoeline.replace('Alpha Orbital Energies', 'Beta Orbital Energies ')
    text['energy'] = [amoeline, bmoeline]
    amocline = text['coeff'][0]
    bmocline = amocline.replace('Alpha MO coefficients', 'Beta MO coefficients ')
    text['coeff'] = [amocline, bmocline]
    
    
    # print '%s generated.' % ofn
    f = open(ofn, 'w')
    write = f.write
    writes = lambda s: write(s + '\n')
    writel = lambda l: write('\n'.join(l) + '\n')
    
    writes(title)
    writel(text['info'][1:])
    for spin, energy in enumerate(data['energy']):
        writes(text['energy'][spin])
        for i, e in enumerate(['%16.8E' % e for e in energy]):
            write(e)
            if (i + 1) % 5 == 0:
                write('\n')
        if (i + 1) % 5 != 0:
            write('\n')
    for spin, coeff in enumerate(data['coeff']):
        writes(text['coeff'][spin])
        # print coeff
        for i, e in enumerate(['%16.8E' % e for e in coeff.reshape(1,-1)[0]]):
            write(e)
            if (i + 1) % 5 == 0:
                write('\n')
        if (i + 1) % 5 != 0:
            write('\n')
    writel(text['rest'])
    
    f.close()
    return ofn

def quicksave(mfn, xoao, fmxo, suffix='_rlo', overwrite=False):
    if xoao.ndim == 2:
        if fmxo.ndim == 2:
            energies = np.diag(fmxo)
        elif fmxo.ndim == 1:
            energies = fmxo
        else:
            raise UserWarning('Invalid shape of fock matrix')
        assert energies.shape[0] == xoao.shape[0]
        data = {'energy': [energies], 
                'coeff': [xoao], 
                'dim': xoao.shape}
    elif xoao.ndim == 3:
        if fmxo.ndim == 3:
            energies = np.diagonal(fmxo, axis1=1, axis2=2)
        elif fmxo.ndim == 2:
            energies = fmxo
        else:
            raise UserWarning('Invalid shape of fock matrix')
        assert energies.shape[1] == xoao.shape[1]
        assert energies.shape[0] == 2
        assert xoao.shape[0] == 2
        data = {'energy': energies, 
                'coeff': xoao, 
                'dim': xoao.shape[1:]}
    else:
        raise UserWarning('Unrecognized data shape: %s' % xoao.shape)
    
    return gen_fchk(mfn, data, suffix=suffix, overwrite=overwrite)

if __name__ == '__main__':
    import sys
    mfn = sys.argv[1]
    gen_fchk(mfn, parse_fchk(mfn)['data'])

