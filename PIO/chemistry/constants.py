Bohr_Radius = 0.52917721067

class _AtomList():
    '''Class _AtomList Contains basic information of atoms,
    including atomic numbers and symbols.'''
    def __init__(self):
        self.atoms = ('X',
            'H',                                    'He',
            'Li','Be',     'B', 'C', 'N', 'O', 'F', 'Ne',
            'Na','Mg',     'Al','Si','P', 'S', 'Cl','Ar',
            'K', 'Ca',
                 'Sc','Ti','V', 'Cr','Mn','Fe','Co','Ni','Cu','Zn',
                           'Ga','Ge','As','Se','Br','Kr',
            'Rb','Sr',
                 'Y', 'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
                           'In','Sn','Sb','Te','I', 'Xe',
            'Cs','Ba',
                 'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd',
                      'Tb','Dy','Ho','Er','Tm','Yb','Lu',
                      'Hf','Ta','W', 'Re','Os','Ir','Pt','Au','Hg',
                           'Tl','Pb','Bi','Po','At','Rn',
            'Fr','Ra',
                 'Ac','Th','Pa','U', 'Np','Pu','Am','Cm'
            )
        self.atomlist = [atom.lower() for atom in self.atoms]
    
    def symbol2an(self, atomic_symbol):
        '''input: atomic symbol (str)
        output: atomic number (int)'''
        try:
            return self.atomlist.index(atomic_symbol.strip().lower())
        except ValueError:
            raise Exception('Atom symbol unrecognized.')
    
    def an2symbol(self, an):
        '''an2symbol(an)
        input: atomic number (int)
        output: atomic symbol (str)'''
        try:
            return self.atoms[an]
        except IndexError:
            raise Exception('Atomic number out of range.')
    
ATOMLIST = _AtomList()

