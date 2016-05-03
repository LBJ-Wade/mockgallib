mockgallib
==========

A Library for generating mock galaxy catalogues using C++/Python (cosmology)

[https://github.com/junkoda/mockgallib](https://github.com/junkoda/mockgallib)

[http://junkoda.github.io/codes/mockgallib/](http://junkoda.github.io/codes/mockgallib/)

## Compile and install

### Library `libs`

Run `make` in `libs` directory. GNU Scientific Library (GSL) and HDF5
are required. If you have these libraries in non-standard directories,
add `DIR_PATH` in `libs/Makefile`.

```bash
% make
```

This creates a static library `libs/libmockgal.a` and a shared library
`libs/libmockgal.so`. You can copy the shared library to the standard
location (e.g., /usr/local/lib ), or set an environment variable:

```bash
export LD_LIBRARY_PATH=/path/to/mockgallib/libs:$LD_LIBRARY_PATH
```


### Python module `python`

Run setup.py in the `python` directory:

Python3 and NumPy are required. If you do not have python3 and do not
have privileges to install Python3, I recommend using the
[anaconda/miniconda](https://www.continuum.io/downloads) package.


```bash
% python3 setup.py build_ext --inplace
```

You can now load `mockgallib` package in the `python/` directory

```bash
% python3
>>> import mockgallib
```

You can also install the module so that you can import from any location:

```bash
% python3 setup.py install
```

