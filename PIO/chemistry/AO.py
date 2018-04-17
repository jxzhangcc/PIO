class AtomicOrbital(object):
    Llabels = ['s', 'p', 'd', 'f']
    def __init__(self, center, label):
        self.center = center
        label = int(label)
        self.n = 0
        self.l = label / 100
        self.pure = bool((label % 100) / 50) or self.l <= 1
        self.m = label % 50
        # AO label = 100*l+k+m (See NBO mannual B-74)
        # l = angular quantum number
        # k = 0 for cartesian and 50 for pure
        # m = magnetic quantum number
        # AOLabels = {0: ['s'], 
                    # 10: ['px', 'py', 'pz'], 
                    # 20: ['dxx', 'dxy', 'dxz', 'dyy', 'dyz', 'dzz'], 
                    # 25: ['dxy', 'dxz', 'dyz', 'dx2y2', 'dz2'], 
                    # 35: ['f(0)', 'f(c1)', 'f(s1)', 'f(c2)', 
                        # 'f(s2)', 'f(c3)', 'f(s3)']}
        # int2label = lambda i: AOLabels[i / 10][i % 10 - 1]
        # labels = [int2label(int(_)) for _ in l.strip().split()]
    
    def set_n(self, n):
        self.n = n
    
    @classmethod
    def Llabel(self, l):
        assert l >= 0
        if l < len(AtomicOrbital.Llabels):
            return AtomicOrbital.Llabels[l]
        return chr(ord(AtomicOrbital.Llabels[-1]) + l + 1
            - len(AtomicOrbital.Llabels))
    
    def label(self):
        return (str(self.n) if self.n else '') + \
            self.Llabel(self.l) + str(self.m)
    
    def __str__(self):
        return str(self.center) + '-' + self.label()
    
    def shift(self, delta):
        self.center += delta
        return self
    
    def __eq__(self, ao2):
        return (self.center == ao2.center and self.n == ao2.n 
            and self.l == ao2.l and self.m == ao2.m and self.pure == ao2.pure)

def main():
    pass

if __name__ == '__main__':
    main()

