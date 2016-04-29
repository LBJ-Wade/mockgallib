Remap
=====

Remaps a periodic cube to a cuboid.

<span class="nc">class</span> **Remap**(*u*, *boxsize*)

> - *u* is a 3x3 integer matrix that characterise the remapping.
> - *boxsize* is the boxsize of the cube before remapping

<span class="kt">bool</span> remap&rarr;corrdinate(<span class="kt">Halo</span>* h)

> - remap h&rarr;x to the cuboid corrdinate
> - returns true when the remapping is sucessful
> - returns false and h&rarr;x = 0 if remapping not found
