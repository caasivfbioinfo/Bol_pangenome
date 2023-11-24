##How to install and use the mSynOrths 
###Install:
``````
conda create -n msynorths python=3.7
git clone https://gitee.com/zhanglingkui/msynorths.git
pip install networkx==2.2
pip install numpy==1.21.6
pip install matplotlib==3.0.2
pip install scikit-learn==1.0.2
pip install scipy==1.1.0
``````
###Usage

```
mSynOrths.py [-h] [-f F] [-c] [-a] [-o O] [-t T] [-s S] [-n N] [-e E]
                    [-g G]

mSynOrths_0.1

optional arguments:
  -h, --help        show this help message and exit
  -f F              input genome data location file
  -c, --continueDo  If the program breaks,you can use this parameter to
                    continue
  -a, --addGenome   Add the genome to the existing results
  -o O              output folder
  -t T              input num of threads
  -s S              blast or diamond
  -n N              The number of flanking homologous gene pairs supporting
                    syntenic gene pairs.
  -e E              blast Evalue
  -g G              The maximum number of genes in the syntenic fragment gap
```