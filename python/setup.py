from distutils.core import setup, Extension
import numpy as np

print("np.get_include()")
print(np.get_include())

setup(name='mockgallib',
      version='0.0.1',
      author='Jun Koda',
      py_modules=['mockgallib.power', 'mockgallib.const', 'mockgallib.sigma',
                  'mockgallib.mass_function', 'mockgallib.hod',
                  'mockgallib.lightcones', 'mockgallib.catalogues',
                  'mockgallib.sky'
      ],
      ext_modules=[
          Extension('mockgallib._mockgallib',
                    ['py_package.cpp', 'py_msg.cpp', 'py_const.cpp',
                     'py_cosmology.cpp',
                     'py_power.cpp',  'py_growth.cpp',
                     'py_sigma.cpp', 'py_mf.cpp', 'py_hod.cpp',
                     'py_nbar.cpp', 'py_nbar_fitting.cpp',
                     'py_lightcones.cpp', 'py_catalogues.cpp',
                     'py_sky.cpp'
                    ],
                    include_dirs = ['../libs', np.get_include()],
                    library_dirs =  ['../libs'],
                    libraries = ['mockgal'],
          )
      ],
      packages=['mockgallib'],
)


