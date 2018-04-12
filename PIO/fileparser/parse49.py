import os
import sys
import re
import numpy as np
from chemistry.AO import AtomicOrbital
# from chemistry import ATOMLIST
# from fileparser.parserException import EmptyFile, BadFile, FileExtNotMatch, NonDataFile, NoFile

def nbofn(f):
    dir, fn = os.path.split(f)
    title = os.path.splitext(fn)[0].upper()
    if len(title) < 4:
        title += 'FILE'[len(title)-4:]
    return os.path.join(dir, title)

def smallest_factor_larger_than_square_root(n):
    i = int((n-1) ** 0.5) + 1
    while n % i:
        i += 1
    return i

def tri2square(tri, dim = None):
    assert tri.ndim == 1
    if not dim:
        dim = int((tri.shape[0]*2) ** 0.5)
    assert tri.shape == (dim*(dim+1)/2,)
    square = np.zeros((dim, dim))
    id = np.tril_indices(dim)
    square[id] = tri
    square[id[::-1]] = tri
    return square

def parse_49(filename):
    with open(filename) as f:
        text = f.read().splitlines()
    print text[0]
    cuts = []
    for k, line in enumerate(text):
        if re.match(r'^ -+$', line):
            cuts.append(k-2)
    cuts.append(None)
    sections = [text[start:end]for start, end in zip(cuts[:-1], cuts[1:])]
    d = dict()
    for sec in sections:
        title = sec[0].strip()
        d.setdefault('title', title)
        type = sec[1].strip().rstrip(':')
        
        # print 'Read in: %s of %s' % (type, title)
        data_text = []
        flag_spin = False
        for k, line in enumerate(sec):
            if k < 3:
                continue
            if re.match(r'^(alpha|beta ) spin$', line.strip()):
                flag_spin = True
                continue
            if re.match(r'^( +-?\d+\.\d+){,5}$', line):
                data_text.append(line)
            else:
                break
        else:
            k += 1
        label_text = '\n'.join(sec[k:])
        data = np.fromstring(''.join(data_text), sep = ' ')
        if re.match(r'^[A-Z]+ (overlap|density|Fock) matrix$', type):
            m = re.match(r'^([A-Z]+) (overlap|density|Fock) matrix$', type)
            if not flag_spin:
                if m.group(2) == 'overlap':
                    d['S'+m.group(1)] = tri2square(data)
                if m.group(2) == 'density':
                    d['D'+m.group(1)] = tri2square(data)
                if m.group(2) == 'Fock':
                    d['F'+m.group(1)] = tri2square(data)
            else:
                if m.group(2) == 'overlap':
                    d['S'+m.group(1)] = np.array([tri2square(_) for _ in data.reshape(2,-1)])
                if m.group(2) == 'density':
                    d['D'+m.group(1)] = np.array([tri2square(_) for _ in data.reshape(2,-1)])
                if m.group(2) == 'Fock':
                    d['F'+m.group(1)] = np.array([tri2square(_) for _ in data.reshape(2,-1)])
        elif re.match(r'^[A-Z]+s in the [A-Z]+ basis$', type):
            assert not flag_spin
            m = re.match(r'^([A-Z]+)s in the ([A-Z]+) basis$', type)
            dimxo = smallest_factor_larger_than_square_root(data.size)
            d[''.join(m.groups())] = data.reshape(-1, dimxo)
            if m.groups() == ('NAO', 'AO'):
                d['NAOlabel'] = analyze_NAO_labels(label_text, dimxo)
            else:
                d[''.join(m.groups())+'_labeltext'] = label_text
        else:
            print type, data.shape
    return d

class NaturalAtomicOrbital(AtomicOrbital):
    def __init__(self, center, label, Lewis_flag):
        AtomicOrbital.__init__(self, center, label)
        self.Lewis_flag = bool(Lewis_flag)

def analyze_NAO_labels(label_text, dim):
    n_sections = 3
    l = np.fromstring(' '.join(label_text.strip().split()), sep=' ', dtype=int)
    l.shape = (3, -1)
    return [NaturalAtomicOrbital(*orb) for orb in l.T]

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

if __name__ == '__main__':
    parse_49(sys.argv[1])

