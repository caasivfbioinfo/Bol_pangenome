###This bash could be used to visualize the TAD structure of genome.

Prerequisite: HDF5

If you are on Linux, download the source code of the latest version from the HDF5 website and unpack it.
``````
# make a new directory
mkdir hdf5-build
cd hdf5-build
# replace xx with current version number
wget https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.xx.tar.gz
# unpack
tar xzf hdf5-1.8.xx.tar.gz
cd hdf5-1.8.xx/
# use --prefix to set the folder in which HDF5 should be installed
# alternatively, you can omit --prefix=... here and run
# sudo make install to install globally (requires admin rights)
./configure --prefix=/path/to/hdf5/dir
make
make install
``````

FAN-C is a Python (3.6+) toolkit for the analysis and visualisation of Hi-C data
``````
conda activate python3
``````
The simplest way to install FAN-C is via pip:
``````
pip install --user fanc
``````