import mockgallib._mockgallib as c

class Const:
    def __init__(self):
        self.G = c._const_G()
        self.rho_crit_0 = c._const_rhocrit0()

    def __repr__(self):
        return ("Constants in internal unit:\n" +
                "  G          = %.4e\n" +
                "  rho_crit_0 = %.4e\n") % (self.G, self.rho_crit_0) 
