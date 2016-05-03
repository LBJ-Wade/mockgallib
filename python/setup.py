from distutils.core import setup, Extension
import numpy as np

setup(name='mockgallib',
      version='0.0.1',
      description='Mock galaxy library',
      author='Jun Koda',
      url='https://github.com/junkoda/mockgallib',
      py_modules=['mockgallib.power', 'mockgallib.const', 'mockgallib.sigma',
                  'mockgallib.mf', 'mockgallib.hod',
                  'mockgallib.lightcones', 'mockgallib.catalogues',
                  'mockgallib.sky', 'mockgallib.distance',
                  'mockgallib.remap'
      ],
      ext_modules=[
          Extension('mockgallib._mockgallib',
                    ['py_package.cpp', 'py_msg.cpp', 'py_const.cpp',
                     'py_cosmology.cpp',
                     'py_power.cpp',  'py_growth.cpp',
                     'py_sigma.cpp', 'py_mf.cpp', 'py_hod.cpp',
                     'py_nbar.cpp', 'py_nbar_fitting.cpp',
                     'py_lightcones.cpp', 'py_catalogues.cpp',
                     'py_sky.cpp', 'py_distance.cpp',
                     'py_remap.cpp'
                    ],
                    include_dirs = ['../libs', np.get_include()],
                    library_dirs =  ['../libs'],
                    libraries = ['mockgal'],
          )
      ],
      packages=['mockgallib'],
)


