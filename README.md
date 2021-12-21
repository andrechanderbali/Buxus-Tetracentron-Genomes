# Buxus-Tetracentron-Genomes
## Scripts and command-line arguments
### Genome Assembly

```
fc_run.py pb-asm.cfg
```

### phase and polish
```
fc_unzip.py pb-unzip.cfg 
```

#### pb-unzip.cfg
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
