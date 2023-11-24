##一、评估基因组大小及杂合度

```
jellyfish count -C -t 24 -m 21 -s 4G -o kmer21.out 1.clean.fq 2.clean.fq
jellyfish histo -t 24 kmer21.out -o kmer21.histo
genomescope.R kmer21.histo 21 150
```

##二、NextDenovo进行三代数据组装
NextDenovo需要用Python2.7
####1.用seq_stat确定seed_cutoff值，NextDenovo配置文件需要该值
```
~/NextDenovo/bin/seq_stat -g 568Mb -d 45 input.fofn >seq_stat.log
```
####2.配置NextDenovo运行文件
复制配置文件run.cfg：cp ~/NextDenovo/doc/run.cfg ./run_NextDenovo.cfg
修改配置文件run_NextDenovo.cfg：
input_NextDenovo.fofn放三代数据存储的路径，其他参数视情况修改。
```
[General]
job_type = local
job_prefix = nextDenovo
task = all # 'all', 'correct', 'assemble'
rewrite = yes # yes/no
deltmp = yes
rerun = 3
parallel_jobs = 8
input_type = raw
input_fofn = input_NextDenovo.fofn
workdir = 01_rundir

[correct_option]
read_cutoff = 1k
seed_cutoff = 37159
blocksize = 1g
pa_correction = 2
seed_cutfiles = 30
sort_options = -m 1g -t 2 -k 50
minimap2_options_raw = -x ava-ont -t 5
correction_options = -p 15

[assemble_option]
random_round = 20
minimap2_options_cns = -x ava-ont -t 5 -k17 -w17
nextgraph_options = -a 1
```
####3.运行NextDenovo
```
nextDenovo run_NextDenovo.cfg
```
###三、利用二代，三代数据纠错
####1.准备配置文件
（1）创建文件，记录二代序列的位置：realpath  1.fastq  2.fastq  > sgs.fofn （需解压为fq）
（2）复制配置文件run.cfg：cp NextPolish/doc/run.cfg ./run_nextPolish.cfg
（3）修改配置文件：
``````
[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = 3
multithread_jobs = 10
genome = ../../racon/racon_r3.fasta
genome_size = auto
workdir = ../racon_nextpolish/01_rundir
polish_options = -p {multithread_jobs}

[sgs_option] #2代数据
sgs_fofn = ./sgs.fofn
sgs_options = -max_depth 100

[lgs_option] #3代数据
lgs_fofn = ./lgs.fofn
lgs_options = -min_read_len 10k -max_read_len 150k -max_depth 60
lgs_minimap2_options = -x map-ont
``````
（sgs_option，lgs_option若同时设置，可实现对2代+3代同时纠错）

####2.运行nextPolish
```
nextPolish run_nextPolish.cfg
```
