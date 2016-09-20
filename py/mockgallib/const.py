import mockgallib._mockgallib as c

G = None
rho_crit_0 = None
delta_c = None

def init():
    global G, rho_crit_0, delta_c
    G = c._const_G()
    rho_crit_0 = c._const_rhocrit0()
    delta_c = c._const_deltac();

