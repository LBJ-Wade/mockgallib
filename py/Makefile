# make
#   to compile this python module
# make install
#   to install it to a python module directory

all: build install


SRC  := py_package.cpp py_power.cpp py_msg.cpp

build: $(SRC)
	python3 setup.py build_ext --inplace

# You might need extra options for libraries in non-standard location:
# python3 setup.py build_ext --swig-opts="-I... -L..."

install:
	python3 setup.py install


.PHONY: clean

clean:
	rm -rf build mockgallib/*.so

