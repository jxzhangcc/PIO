print 'Parse FChk'
import os
import sys
import re
import numpy as np

def parse_fchk(filename):
    if os.path.splitext(filename)[1].lower() not in ('.fchk', '.fch'):
        raise Exception('%s is not a fchk file.' % filename)
    
    f = open(filename)
    newline = lambda: f.readline().strip('\r\n')
    text = {'info': [], 
            'energy': [], 
            'coeff': [], 
            'rest': []}
    data = {'energy': [], 
            'coeff': [], 
            'dim': (0, 0)}
    data_pat = re.compile(r'^( [ -]\d\.\d{8}E[+-]\d{2}){,5}$')
    
    line = newline()
    while 'Orbital Energies' not in line:
        text['info'].append(line)
        line = newline()
    
    while 'Orbital Energies' in line:
        text['energy'].append(line)
        sec = []
        line = newline()
        while data_pat.match(line):
            sec.append(line)
            line = newline()
        data['energy'].append(np.fromstring(''.join(sec), sep = ' '))
    
    assert 'MO coefficients' in line
    while 'MO coefficients' in line:
        text['coeff'].append(line)
        sec = []
        line = newline()
        while data_pat.match(line):
            sec.append(line)
            line = newline()
        data['coeff'].append(np.fromstring(''.join(sec), sep = ' '))
    
    while line:
        text['rest'].append(line)
        line = newline()
    
    f.close()
    
    for line in text['info']:
        m1 = re.match(r'^Number of basis functions +I +(\d+)$', line)
        m2 = re.match(r'^Number of independent functions +I +(\d+)$', line)
        if m1:
            dim_bs = int(m1.group(1))
        if m2:
            dim_mo = int(m2.group(1))
    data['dim'] = (dim_mo, dim_bs)
    
    assert len(text['energy']) == len(data['energy'])
    assert len(text['coeff']) == len(data['coeff'])
    for e in data['energy']:
        assert e.shape == (dim_mo,)
    for c in data['coeff']:
        c.shape = (-1, dim_bs)
        assert c.shape == (dim_mo, dim_bs)
    
    if dim_mo < dim_bs:
        print 'Warning! Truncated dimension of orbitals detected.'
    
    return {'data': data, 'text': text}

if __name__ == '__main__':
    fchk = parse_fchk(sys.argv[1])
    for e in fchk['data']['energy']:
        print e.shape
    for c in fchk['data']['coeff']:
        print c.shape
    print fchk['data']['dim']
    for k, v in fchk['text'].items():
        print k
        if raw_input() == 'y':
            print v

