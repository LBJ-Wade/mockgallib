HaloMass
========

<span class="kcl">class</span> **HaloMassFoF**

halo_mass = **HaloMassFoF**(<span class="kt">char</span> *filename[]*,
<span class="kt">double</span> *M0* = 1.0e14)

> Reads a relation between number of FoF particle *nfof* and halo mass. Above *M0*, linear extrapolation is used instead of the tabulated data.

halo_mass&rarr;**mass**(<span class="kt">int</span> *nfof*)

> Returns interpolated halo mass
