SyRI needs Python3.5 and depends on the following python libraries:
```
conda install  python=3.5
conda install cython numpy=1.16. scipy pandas=0.23.4 biopython psutil matplotlib=3.0.0
conda install -c conda-forge python-igraph
conda install -c bioconda pysam
conda install -c bioconda longestrunsubsequence
```

Download Syri:
```
git clone https://github.com/schneebergerlab/syri.git
```

Open the directory where setup.py is (we will call it as cwd) and in the terminal run:
```
python3 setup.py install		            # Install syri
chmod +x syri/bin/syri syri/bin/chroder	syri/bin/plotsr	# Make files executable
```