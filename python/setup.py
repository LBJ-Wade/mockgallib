from distutils.core import setup, Extension
import numpy as np

print("np.get_include()")
print(np.get_include())

setup(name='mockgallib',
      version='0.0.1',
      author='Jun Koda',
      py_modules=['mockgallib.power'],
      ext_modules=[
          Extension('mockgallib._mockgallib',
                    ['py_package.cpp', 'py_power.cpp', 'py_msg.cpp'],
#                    ['py_package.cpp'],
                    include_dirs = ['../libs', np.get_include()],
                    library_dirs =  ['../libs'],
                    libraries = ['mockgal'],
          )
      ],
      packages=['mockgallib'],
)


