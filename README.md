# Buxus-Tetracentron-Genomes
## Scripts and command-line arguments
### Genome Assembly

```
fc_run.py pb-asm.cfg
```

### pb-asm.cfg

###### Input
```
input_fofn=input.fofn
input_type=raw
pa_DBdust_option=
pa_fasta_filter_option=streamed-median
target=assembly
skip_checks=False
LA4Falcon_preload=false
```
###### Data Partitioning
```
pa_DBsplit_option=-x500 -s400
ovlp_DBsplit_option=-s400
```
###### Repeat Masking
```
pa_HPCTANmask_option=
#no-op repmask param set
pa_REPmask_code=0,300;0,300;0,300
```
#### Pre-assembly
```
genome_size = 900000000
length_cutoff = 5000
pa_HPCdaligner_option=-v -B128 -M24
pa_daligner_option= -k18 -e0.80 -l1000 -h480 -w8 -s100
falcon_sense_option=--output-multi --min-idt 0.70 --min-cov 4 --max-n-read 200
falcon_sense_greedy=False
```
#### Pread overlapping
```
ovlp_HPCdaligner_option=-v -B128 -M24
ovlp_daligner_option= -k24 -e.96 -l1800 -h480 -s100
```
#### Final Assembly
```
length_cutoff_pr=1000
overlap_filtering_setting=--max-diff 100 --max-cov 100 --min-cov 2
fc_ovlp_to_graph_option=

[job.defaults]
job_type=slurm
pwatcher_type=
JOB_QUEUE=
MB=6000
NPROC=8
njobs=24
submit = sbatch --time=12:00:00 -J ${JOB_NAME} --mem-per-cpu=${MB}M --cpus-per-task=${NPROC} ${JOB_SCRIPT}

[job.step.da]
NPROC=8
MB=6000
njobs=24

[job.step.la]
NPROC=8
MB=6276
njobs=24

[job.step.cns]
NPROC=8
MB=6553
njobs=24

[job.step.pda]
NPROC=4
MB=6276
njobs=24

[job.step.pla]
NPROC=4
MB=6276
njobs=24

[job.step.asm]
NPROC=16
MB=6000
njobs=1
```

### phase and polish
```
fc_unzip.py pb-unzip.cfg 
```

#### pb-unzip.cfg
```
[General]
max_n_open_files = 300

[Unzip]
input_fofn=input.fofn
input_bam_fofn=input_bam.fofn
polish_include_zmw_all_subreads = true

[job.defaults]
job_type=slurm
pwatcher_type=fs_based
JOB_QUEUE=serial
MB=6000
NPROC=8
njobs=24
TIME=24:00:00
submit = sbatch --time=${TIME} -J ${JOB_NAME} --mem-per-cpu=${MB}M --cpus-per-task=${NPROC} ${JOB_SCRIPT}

[job.step.unzip.track_reads]
njobs=1
NPROC=8
MB=6000

[job.step.unzip.blasr_aln]
njobs=24
NPROC=8
MB=6000

[job.step.unzip.phasing]
njobs=24
NPROC=8
MB=6000

[job.step.unzip.hasm]
njobs=1
NPROC=8
MB=10000
TIME=48:00:00

[job.step.unzip.quiver]
njobs=24
NPROC=8
MB=10000
TIME=00:30:00
```
## Genome Annotation

### annotate using evidence from transcriptome assemblies, reference proteins and a custom repeat library
```
maker --genome Genome.fsa evidence_maker_opts.ctl -cpus 24
```
#### evidence_maker_opts.ctl settings
```
est=Genome.RNASeq.assembly
protein=referenceproteomes.fasta
rmlib=Genome.Repeats.lib
est2genome=1
protein2genome=1 
```
### collect evidence-based annotations
```
awk '{ if ($2 == "est2genome") print $0 }' Genome.evidence.maker.gff > Genome.est2genome.gff 

awk '{ if ($2 == "protein2genome") print $0 }' Genome.evidence.maker.gff > Genome.protein2genome.gff 

awk '{ if ($2 == "repeatmasker") print $0 }' Genome.evidence.maker.gff > Genome.repeats.gff

grep -v -e "Satellite" -e "Simple_repeat" -e "-rich" Genome.repeats.gff > Genome.complex.repeats.gff
```

### Train gene predictors:

#### Train SNAP
```
maker2zff : -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl Genome.snap1 . > Genome.snap1.hmm
```
#### Train AUGUSTUS
```
awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' Genome.evidence.maker.gff | awk -v OFS="\t" '{if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000}' | bedtools getfasta -fi Genome.fasta -bed - -fo Genome.evidence.transcripts1000.fasta

busco --in Genome.evidence.transcripts1000.fasta --out rnd1_Augustus --lineage embryophyta_odb9 --mode genome --cpu 6 --long -z --augustus_parameters='--progress=true'

mv run_rnd1_Augustus/augustus_output/retraining_parameters/* Genome.augustus.1
```
### Run maker with trained predictors
```
maker --genome Genome.fsa predict_maker_opts.ctl -cpus 24
```
#### predict_maker_opts.ctl settings
```
est_gff= Genome.est2genome.gff 
protein_gff=Genome.protein2genome.gff
rm_gff=Genome.complex.repeats.gff
snaphmm=snap/Genome.snap1.hmm 
augustus_species=Genome_busco.1

est2genome=0 
protein2genome=0 
```
```
Repeat two more rounds of training and gene prediction 
```

## Transcriptome assembly and gene prediction
```
Trinity --seqType fq --CPU 10 --max_memory 80G \
--left \
./species_sample_R1.fastq.gz \
--right \
./Species_sample_R2.fastq.gz \
--trimmomatic \
--output trinity_out_Species \
--full_cleanup


TransDecoder.LongOrfs -t trinity_out_Species.Trinity.fasta
TransDecoder.Predict -t trinity_out_Species.Trinity.fasta
```

## Phylogenetic analyses

### align protein sequences with MAFFT:
```
mafft --localpair --maxiterate 1000 --thread 1 $file.pep > $file.pep.mafft
```
### Convert protein alignments to cdna alignments:
```
pal2nal.v14/pal2nal.pl $file.pep.mafft $file.cdna -output fasta > $file.pep.mafft.p2n
```
### trim alignments with trimAl
```
for f in alignments/*.pep.mafft.p2n

do
        trimal -in $f -out $f.rmspurious -resoverlap 0.50 -seqoverlap 70
        trimal -in $f.rmspurious -out $f.trim1 -automated1
        trimal -in $f.trim1 -out $f.trim2 -resoverlap 0.50 -seqoverlap 70
        sed -e 's/^\(>[^[:space:]]*\).*/\1/' $f.trim2 > $f.trimal
        rm $f.rmspurious
        rm $f.trim1
        rm $f.trim2
done
```
### Strip gene IDs from fasta header and rename files for FASCONCAT
```
for file in alignments/*.pep.mafft.p2n.trimal; do awk -F"|" '{print $1}' $file > $file.fas; done
```
### generate supermatrix and corresponding partition files
```
perl ~/bin/FASconCAT-G_v1.02.pl -l -l -s 
```
### RAxML with 1000 bootstrap replicates for individual alignments
```
for file in alignments/*.trimal.fas
do
raxmlHPC-PTHREADS-SSE3 -s $file -f a -m GTRGAMMA -p $RANDOM -# 1000 -T 10 -x $RANDOM -n $file
done
```
### ASTRAL to generate coalescence-based tree from individual gene trees
```
cat alignments/RAxML_bipartitions.* > cat.RAxML_bipartitions.genetrees #concatenate gene trees
astral -i cat.RAxML_bipartitions.genetrees -t 1 -o gene.Astral.tre
```
### RAxML with partitioned supermatrix and 1000 bootstrap replicates 
```
raxmlHPC-PTHREADS-SSE3 -s alignment.supermatrix.fas -q alignment.supermatrix_partition.txt -f a -m GTRGAMMA -p $RANDOM -# 1000 -x 1234 -T 8 -n alignment.partition.supermatrix
```
###  MRBAYES with partitioned supermatrix 
```
mpiexec -np 1 mb alignment.supermatrix_partition.nex
```
### IQ-Tree with 1000 bootstrap replicates for individual orthogroups alignments
```
for file in Orthogroup_alignments/*pep.mafft.p2n.trimal
do
iqtree2 -s $file -B 1000 -m TEST -mset raxml -T AUTO
done
```
### Run ASTRAL-Pro on Orthogroup trees
```
strip "|sequenceID" from newick trees leaving species names 

sed 's/|[^:]*:/:/g' OrthogroupGenetrees.nwk > OrthogroupGenetrees4Apro.nwk

java -Xmx3000M -D"java.library.path=~/bin/A-pro/lib" -jar ~/bin/A-pro/astral.1.1.3.jar -i OrthogroupGenetrees4Apro.nwk -o ASTRAL-Pro.GeneFamilies.tree.nwk 2>out.log
```
## Subgenome Synteny 

### Pairwise synteny searches
```
python -m jcvi.compara.catalog ortholog grape buxus

python -m jcvi.compara.catalog ortholog grape tetra 

python -m jcvi.compara.catalog ortholog grape aquil

python -m jcvi.compara.catalog ortholog grape nelumbo

python -m jcvi.compara.catalog ortholog grape grape

python -m jcvi.compara.catalog ortholog grape amaranth
```
### get syntenic blocks with specified syntenic depths
```
python -m jcvi.compara.synteny mcscan grape.bed grape.aquil.lifted.anchors --iter=2 -o grape.aquil.i2.blocks

python -m jcvi.compara.synteny mcscan grape.bed grape.buxus.lifted.anchors --iter=2 -o grape.buxus.i2.blocks

python -m jcvi.compara.synteny mcscan grape.bed grape.grape.lifted.anchors --iter=3 -o grape.grape.i3.blocks

python -m jcvi.compara.synteny mcscan grape.bed grape.nelumbo.lifted.anchors --iter=2 -o grape.nelumbo.i2.blocks

python -m jcvi.compara.synteny mcscan grape.bed grape.tetra.lifted.anchors --iter=4 -o grape.tetra.i4.blocks
```
### merge pairwise syntenic blocks into a multi-species “blocks” file
```
python -m jcvi.formats.base join grape.grape.i3.blocks grape.aquil.i2.blocks grape.buxus.i2.blocks grape.nelumbo.i2.blocks grape.tetra.i4.blocks --noheader | cut -f1,2,3,5,6,8,9,11,12,14,15,16,17,19,20,21 > five_species_merged.blocks
```
## CoGe analyses
```
Comparative genomic analyses conducted on the CoGe platform can be recreated at the following urls: 

Buxus mapped to Vitis: https://genomevolution.org/r/1fksc

Tetracentron mapped to Vitis: https://genomevolution.org/r/1fkss

Aquilegia mapped to Vitis: https://genomevolution.org/r/1apmi 

Nelumbo mapped to Vitis: https://genomevolution.org/r/1fkub 

Nelumbo mapped to Aquilegia: https://genomevolution.org/r/1aplx 

Vitis mapped to Buxus: https://genomevolution.org/r/1fkuu 

Vitis mapped to Tetracentron: https://genomevolution.org/r/1gm8e 

Vitis mapped to Aquilega: https://genomevolution.org/r/1fjcs 

Vitis mapped to Nelumbo: https://genomevolution.org/r/1fjcc 
```
## Ancestral Genomes
```
Code available at https://github.com/jin-repo/RACCROCHE
```
