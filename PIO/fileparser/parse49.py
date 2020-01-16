import os
import sys
import re
from math import ceil, sqrt
import numpy as np
from chemistry.AO import AtomicOrbital
from chemistry.constants import ATOMLIST

ALLOWED_XOS = "AO PAO PNAO NAO PNHO NHO PNBO NBO PNLMO NLMO MO NO"
ALLOWED_MS = "S DM F DI K V"
ALLOWED_MATRICES = \
'''
    AOPAO AOPNAO AOPNHO AOPNBO AOPNLMO AOMO AONO AONAO AONHO AONBO AONLMO
    PAOPNAO
    NAONHO NAONBO NAONLMO NAOMO NAONO
    NHONBO NHONLMO NHOMO NHONO
    NBONLMO NBOMO NBONO
    # NLMOAO NLMOPAO NLMOPNAO NLMONAO NLMOPNHO NLMONHO NLMOPNBO NLMONBO NLMOPNLMO 
    NLMOMO NLMONO
    NOMO
    SAO SPAO SPNAO SPNHO SPNBO SPNLMO
    DMAO DMPNAO DMNAO DMNHO DMNBO DMNLMO DMMO DMNO
    FAO FNAO FPNHO FNHO FPNBO FNBO FNLMO FMO FNO
    DIAO DINAO DINHO DINBO DINLMO DIMO DINO
    KAO KNAO KNHO KNBO KNLMO KMO KNO
    VAO VNAO VNHO VNBO VNLMO VMO VNO
'''

class EmptyFile(UserWarning):
    pass

class WrongType(UserWarning):
    pass

class NBOMatrix(object):
    def __init__(self, title, type, data, label_text, flag_spin=0):
        self.title = title
        self.type = type
        self.data = data
        self.label_text = label_text
        assert flag_spin in (0, 1, 2)
            # restricted, unrestricted alpha, unrestricted beta
        self.flag_spin = flag_spin
    
    def set_labels(self, labels):
        self.labels = labels

class NaturalAtomicOrbital(AtomicOrbital):
    def __init__(self, center, label, Lewis_flag):
        AtomicOrbital.__init__(self, center, label)
        self.Lewis_flag = bool(Lewis_flag)



abbrs = {'overlap':'S', 'density':'D', 'Fock':'F', '1-e potential energy':'V', 'kinetic energy':'K'}
def type2key(type):
    m = re.match(r'^([A-Z]+) (overlap|density|Fock|1-e potential energy|kinetic energy) matrix$', type)
    if m:
        return abbrs[m.group(2)] + m.group(1)
    m = re.match(r'^([A-Z]+)s in the ([A-Z]+) basis$', type)
    if m:
        return ''.join(m.groups())
    m = re.match(r'^([A-Z]+) (x|y|z) dipole integrals$', type)
    if m:
        return 'DI' + m.group(2).upper() + m.group(1)
    return type

def type2class(type):
    if re.match(r'^([A-Z]+)s in the ([A-Z]+) basis$', type):
        return 0
    if re.match(r'^([A-Z]+) (overlap|density|Fock|1-e potential energy|kinetic energy) matrix$', type):
        return 1
    if re.match(r'^([A-Z]+) (x|y|z) dipole integrals$', type):
        return 1
    return -1

def read_49(filename):
    with open(filename) as f:
        lines = f.read().splitlines()
    if not lines:
        raise EmptyFile(filename)
    
    cuts = []
    for k, line in enumerate(lines):
        if re.match(r'^ -+$', line):
            cuts.append(k-2)
    cuts.append(None)
    sections = [lines[start:end] for start, end in zip(cuts[:-1], cuts[1:])]
    
    d = {}
    for sec in sections:
        title = sec[0].strip()
        type = sec[1].strip().rstrip(':')
        
        cuts = []
        for k, line in enumerate(sec):
            if re.match(r'^(alpha|beta ) spin$', line.strip().lower()):
                cuts.append((k, line))
        
        if len(cuts) == 0:
            data_text = []
            for k, line in enumerate(sec):
                if k < 3:
                    continue
                if re.match(r'^( +-?\d+\.\d+){,5}$', line):
                    data_text.append(line)
                else:
                    break
            else:
                k += 1
            label_text = '\n'.join(sec[k:])
            data = np.fromstring(''.join(data_text), sep = ' ')
            d[type2key(type)] = NBOMatrix(title, type, data, label_text, 0)
        elif len(cuts) == 2:
            assert re.match(r'^alpha spin$', cuts[0][1].strip().lower())
            assert re.match(r'^beta  spin$', cuts[1][1].strip().lower())
            cuta, cutb = cuts[0][0], cuts[1][0]
            spinsecs = (sec[cuta+1:cutb], sec[cutb+1:])
            for k, spinsec in enumerate(spinsecs):
                data_text = []
                for l, line in enumerate(spinsec):
                    if re.match(r'^( +-?\d+\.\d+){,5}$', line):
                        data_text.append(line)
                    else:
                        break
                else:
                    l += 1
                label_text = '\n'.join(spinsec[l:])
                data = np.fromstring(''.join(data_text), sep = ' ')
                d[type2key(type)+str(k+1)] = NBOMatrix(title, type, data, label_text, k+1)
        else:
            raise UserWarning('More than two spins found in each section.')
    return d

def test_read49():
    class WrongShape(UserWarning):
        pass
    path = r'.\PIO-test-set\nbo3\xx'
    # path = r'.\PIO-test-set\nbo6\xx'
    dirs = os.listdir(path)
    dim_bs = 62
    dim = dim_bs
    for dir in dirs:
        print dir
        for f in os.listdir(os.path.join(path,dir)):
            # if f in ('CN-U.70', 'CN-RO.70'):
                # continue
            if f.endswith('.nbo') or f.endswith('.47'):
                continue
            try:
                d = read_49(os.path.join(path,dir,f))
            except EmptyFile as e:
                print f
                continue
            try:
                for k, v in d.iteritems():
                    # print k, v.data.shape
                    if type2class(v.type) == 0:
                        if type2key(v.type) in ('NBOAO', 'NLMOAO'):
                            assert v.data.shape == (dim * dim_bs + dim,)
                        else:
                            assert v.data.shape == (dim * dim_bs,)
                    elif type2class(v.type) == 1:
                        assert v.data.shape == (dim_bs * (dim_bs + 1) / 2,)
                    else:
                        raise WrongShape
            except AssertionError as e:
                print f
                raise AssertionError
            except WrongShape as e:
                print f
                raise WrongShape



def nbofn(f):
    dir, fn = os.path.split(f)
    title = os.path.splitext(fn)[0].upper()
    if len(title) < 4:
        title += 'FILE'[len(title)-4:]
    return os.path.join(dir, title)

def find_n_from_n_sq(M, bias=0, warning=True):
    # M = n*(n-d+b), given b, find n for smallest d
    i = int(ceil((sqrt(4*M+bias*bias)-bias)/2))
    while M % i:
        i += 1
    if warning and i - int(ceil((sqrt(4*M+bias*bias)-bias)/2)) >= 1:
        print 'Warning: Large truncation detected with ' \
            'M=%d, n=%d, truncation=%d' % \
            (M, i, i - int(ceil((sqrt(4*M+bias*bias)-bias)/2)))
    return i
    
def find_n_from_n_tri(M):
    i = int((sqrt(8*M+1) - 1) / 2)
    assert M == i * (i + 1) / 2
    return i

def tri2square(tri):
    assert tri.ndim == 1
    dim = int((tri.shape[0]*2) ** 0.5)
    assert tri.shape == (dim*(dim+1)/2,)
    square = np.zeros((dim, dim))
    id = np.tril_indices(dim)
    square[id] = tri
    square[id[::-1]] = tri
    return square

def parse_XOAO(nbomatrix):
    assert re.match(r'^([A-Z]+)s in the ([A-Z]+) basis$', nbomatrix.type)
    factor = 2 if nbomatrix.flag_spin else 1
    bias = 1 if type2key(nbomatrix.type) in ('NBOAO', 'NLMOAO') else 0
    dim_bs = find_n_from_n_sq(nbomatrix.data.size/factor, bias)
    
    data = nbomatrix.data.reshape()
    
    data = nbomatrix.data
    labels = None
    # if type2key(nbomatrix.type) == 'NAOAO':
        # labels = analyze_NAO_labels(nbomatrix.label_text)
    # elif type2key(nbomatrix.type) == 'NAOAO':
        # labels = analyze_NBO_labels(nbomatrix.label_text)
    return data, labels
    # if re.match(r'^[A-Z]+ (overlap|density|Fock) matrix$', type):
            # m = re.match(r'^([A-Z]+) (overlap|density|Fock) matrix$', type)
            # if not flag_spin:
                # if m.group(2) == 'overlap':
                    # d['S'+m.group(1)] = tri2square(data)
                # if m.group(2) == 'density':
                    # d['D'+m.group(1)] = tri2square(data)
                # if m.group(2) == 'Fock':
                    # d['F'+m.group(1)] = tri2square(data)
            # else:
                # if m.group(2) == 'overlap':
                    # d['S'+m.group(1)] = np.array([tri2square(_) for _ in data.reshape(2,-1)])
                # if m.group(2) == 'density':
                    # d['D'+m.group(1)] = np.array([tri2square(_) for _ in data.reshape(2,-1)])
                # if m.group(2) == 'Fock':
                    # d['F'+m.group(1)] = np.array([tri2square(_) for _ in data.reshape(2,-1)])
        # elif re.match(r'^[A-Z]+s in the [A-Z]+ basis$', type):
            # assert not flag_spin
            # m = re.match(r'^([A-Z]+)s in the ([A-Z]+) basis$', type)
            # dimxo = find_n_0(data.size)
            # d[''.join(m.groups())] = data.reshape(-1, dimxo)
            # if m.groups() == ('NAO', 'AO'):
                # d['NAOlabel'] = analyze_NAO_labels(label_text, dimxo)
            # else:
                # d[''.join(m.groups())+'_labeltext'] = label_text
        # else:
            # print type, data.shape
    pass

def analyze_NAO_labels(label_text, dim):
    n_sections = 3
    l = np.fromstring(' '.join(label_text.strip().split()), sep=' ', dtype=int)
    l.shape = (3, -1)
    labels = [NaturalAtomicOrbital(*orb) for orb in l.T]
    counter = dict()
    for nao in labels:
        label = str(nao)
        if label in counter:
            counter[label] += 1
        else:
            counter[label] = 1
        nao.set_n(counter[label])
    return labels

def analyze_NBO_labels(label_text, dim):
# #Atom#		   1  2  3  4  5  6  1  1  2  2  3  3  4  1  2  3  4  5  6  1  2  3  4  4  5  5  5  6  6  6  1  2  3  4  5  6  1  1  1  1  1  1  1  1  1  1  2  2  2  2  2  2  2  2  2  2  3  3  3  3  3  3  3  3  3  3  4  4  4  4  4  4  4  4  4  4  5  5  5  5  5  5  5  5  5  5  6  6  6  6  6  6  6  6  6  6
# Label#	   1  1  1  1  1  1 11 11 11 11 11 11 11  2  2  2  2  2  2 21 21 21 21 21 21 21 21 21 21 21 22 22 22 22 22 22 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21
# Label	  CR CR CR CR CR CR LP LP LP LP LP LP LP BD BD BD BD BD BD LV LV LV LV LV LV LV LV LV LV LV BD BD BD BD BD BD RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY RY
# Anti?	                                                                                           *  *  *  *  *  *                                                                                                                                                                                      
# #NBO_At	   1  1  1  1  1  1  1  2  1  2  1  2  1  1  1  1  1  1  1  1  1  1  1  2  1  2  3  1  2  3  1  1  1  1  1  1  1  2  3  4  5  6  7  8  9  1  1  2  3  4  5  6  7  8  9  1  1  2  3  4  5  6  7  8  9  1  1  2  3  4  5  6  7  8  9  1  1  2  3  4  5  6  7  8  9  1  1  2  3  4  5  6  7  8  9  1
# #Atom1	   1  3  5  7  9 11  1  1  3  3  5  5  7  1  3  5  7  9 11  1  3  5  7  7  9  9  9 11 11 11  1  3  5  7  9 11  1  1  1  1  1  1  1  1  1  2  3  3  3  3  3  3  3  3  3  4  5  5  5  5  5  5  5  5  5  6  7  7  7  7  7  7  7  7  7  8  9  9  9  9  9  9  9  9  9 10 11 11 11 11 11 11 11 11 11 12
# #Atom2	   0  0  0  0  0  0  0  0  0  0  0  0  0  2  4  6  8 10 12  0  0  0  0  0  0  0  0  0  0  0  2  4  6  8 10 12  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
    dim = int(dim)
    # return len(''.join([line[1:] for line in label_text.strip('\n').split('\n')]))
    text = ''.join([line[1:] for line in label_text.strip('\n').split('\n')])
    n_sections = 14
    label_text = []
    for section in range(n_sections):
        label_text.append([])
        # print text[section * dim * 3 : (section + 1) * dim * 3 ]
        for nbo in range(dim):
            label_text[-1].append(text[section * dim * 3 + nbo * 3:
                section * dim * 3 + nbo * 3 + 3].strip())
    # for section in label_text:
        # print '\t'.join(section)
    atoms = [int(_) for _ in text[n_sections * dim * 3:].strip().split()]
    atoms.append(0)
    otype, antiflag, n_nbo, at1, at2 = label_text[2:7]
    tostring = lambda num: ATOMLIST.an2symbol(atoms[int(num)-1]) + (num if int(num) else '')
    labels = [tostring(at1[nbo]) + tostring(at2[nbo]) + 
        otype[nbo] + antiflag[nbo] + n_nbo[nbo] for nbo in range(dim)]
    # labels = [(atoms[int(at1[nbo])-1] , atoms[int(at2[nbo])-1]) for nbo in range(dim)]
    # labels = [(at1[nbo], at2[nbo]) for nbo in range(dim)]
    # for label in labels:
        # print label
    return (labels, atoms)

def analyze_NBONAO_labels(label_text, dim):
    pass

def analyze_PNAOPAO_labels(label_text, dim):
    pass

def parse_49(filename):
    d = read_49(filename)
    for k, v in d.iteritems():
        if type2class(v.type) == 0:
            if type2key(v.type) in ('NBOAO', 'NLMOAO'):
                dim_bs = find_n_from_n_sq(v.data.size, 1)
                dim = v.data.size / (dim_bs + 1)
                assert v.data.shape == (dim * dim_bs + dim,)
                v.data = v.data[:dim*dim_bs].reshape(dim, dim_bs)
            else:
                dim_bs = find_n_from_n_sq(v.data.size, 0)
                dim = v.data.size / dim_bs
                assert v.data.shape == (dim * dim_bs,)
                v.data.shape = (dim, dim_bs)
            if type2key(v.type) == 'NAOAO':
                v.set_labels(analyze_NAO_labels(v.label_text, dim))
            elif type2key(v.type) in ('NBOAO', 'NLMOAO'):
                v.set_labels(analyze_NBO_labels(v.label_text, dim))
            elif type2key(v.type) == 'NBONAO':
                v.set_labels(analyze_NBONAO_labels(v.label_text, dim))
            elif type2key(v.type) == 'PNAOPAO':
                v.set_labels(analyze_PNAOPAO_labels(v.label_text, dim))
            else:
                assert not v.label_text
        elif type2class(v.type) == 1:
            dim_bs = find_n_from_n_tri(v.data.size)
            assert v.data.shape == (dim_bs * (dim_bs + 1) / 2,)
            assert not v.label_text
            v.data = tri2square(v.data)
        else:
            raise WrongType(filename, k)
    
    ans = {}
    for k, v in d.iteritems():
        ans[k] = v.data
        try:
            ans[k+'label'] = v.labels
        except AttributeError:
            pass
    return ans

def test_parse49():
    class WrongShape(UserWarning):
        pass
    # path = r'.\PIO-test-set\nbo3\xx'
    path = r'.\PIO-test-set\nbo6\xx'
    dirs = os.listdir(path)
    dim_bs = 62
    dim = dim_bs
    for dir in dirs:
        print dir
        for f in os.listdir(os.path.join(path,dir)):
            # if f in ('CN-U.70', 'CN-RO.70'):
                # continue
            if f.endswith('.nbo') or f.endswith('.47'):
                continue
            try:
                d = parse_49(os.path.join(path,dir,f))
            except EmptyFile as e:
                print f
                continue

def main():
    parse_49(sys.argv[1])

if __name__ == '__main__':
    # test_read49()
    # test_parse49()
    main()
    pass

