# Buxus-Tetracentron-Genomes
## Scripts and command-line arguments
### Genome Assembly
##### assemble PacBio reads with parameters specified in pb-asm.cfg
```
fc_run.py pb-asm.cfg 
```
##### phase and polish assembly with parameters specified in pb-unzip.cfg
```
fc_unzip.py pb-unzip.cfg 
```

## Genome Annotation

#### Prepare custom repeat library following the MAKER-P tutorial 
```
http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/Repeat_Library_Construction-Advanced
```
#### Initial MAKER run: evidence-based annotation using transcriptomes, reference proteins, and repeat library
```
maker --genome Genome.fsa evidence_maker_opts.ctl #map ESTs, reference proteins, repeats and generate initial gene models
```
#### evidence_maker_opts.ctl settings
```
est=Genome.RNASeq.assembly
protein=referenceproteomes.fasta
rmlib=Genome.Repeats.lib
est2genome=1
protein2genome=1 
```
#### collect evidence-based annotations
```
awk '{ if ($2 == "est2genome") print $0 }' Genome.evidence.maker.gff > Genome.est2genome.gff 

awk '{ if ($2 == "protein2genome") print $0 }' Genome.evidence.maker.gff > Genome.protein2genome.gff 

awk '{ if ($2 == "repeatmasker") print $0 }' Genome.evidence.maker.gff > Genome.repeats.gff

grep -v -e "Satellite" -e "Simple_repeat" -e "-rich" Genome.repeats.gff > Genome.complex.repeats.gff
```

#### Train gene predictors:

##### Train SNAP
```
maker2zff : -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl Genome_snap1 . > Genome_snap1.hmm
```
##### Train AUGUSTUS
```
awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' Genome.evidence.maker.gff | awk -v OFS="\t" '{if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000}' | bedtools getfasta -fi Genome.fasta -bed - -fo Genome.evidence.transcripts1000.fasta

busco --in Genome.evidence.transcripts1000.fasta --out rnd1_Augustus --lineage embryophyta_odb9 --mode genome --cpu 6 --long -z --augustus_parameters='--progress=true'

mv run_rnd1_Augustus/augustus_output/retraining_parameters/* Genome.augustus.1
```
### Prediction MAKER run
```
maker --genome Genome.fsa predict_maker_opts.ctl #annotate masked assembly with trained predictors and evidence-based gffs
```
#### predict_maker_opts.ctl settings
```
est_gff= Genome.est2genome.gff 
protein_gff=Genome.protein2genome.gff
rm_gff=Genome.complex.repeats.gff
snaphmm=Genome_snap1.hmm 
augustus_species=Genome_busco.1

est2genome=0 
protein2genome=0 
```

#### Repeat two more rounds of training and gene prediction 
```
est_gff= Genome.est2genome.gff 
protein_gff=Genome.protein2genome.gff
rm_gff=Genome.complex.repeats.gff
snaphmm=Genome_snap2.hmm 
augustus_species=Genome_busco.2

est_gff= Genome.est2genome.gff 
protein_gff=Genome.protein2genome.gff
rm_gff=Genome.complex.repeats.gff
snaphmm=Genome_snap3.hmm 
augustus_species=Genome_busco.3

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
### Build data sets
#### Angiosperms353
```

#This script is modified from https://github.com/HeatherKates/GenomicDataforPhylo/blob/master/Blast_to_seqs.bash
#It uses blast to extract a set of target loci (queries) from a set of proteomes (in a dir called proteomes).
#The output is one file per proteome that contains the best blast hit per query.
 
#loop through a set of proteomes in a dir 
for i in path/to/proteomes/*.pep
do
#make blastdb (comment next line out if db already exists)
makeblastdb -dbtype prot -in $i -parse_seqids
#Set variable $Organism equal to the basename of the proteome file.
Organism=$(basename $i .pep)
#blast a list of target sequences (in Queries.fasta, change name if needed) against db
blastp -db $i -query path/to/Amborella.Angiosperms353.FAA -outfmt 6 -max_target_seqs 1 > $Organism.bls
#"get best blast hits per query"
perl get_best_blast_results.pl 1 $Organism.bls $Organism.bls.best
#convert blast table (-outfmt 6) to bed file format
grep -v '^#' ${Organism}.bls.best | perl -ane 'if($F[8]<=$F[9]){print join("\t",$F[1],$F[8]-1,$F[9],$F[0],"0","+"),"\n";}else{print join("\t",$F[1],$F[9]-1,$F[8],$F[0],"0","-"),"\n";}' | sort >> $Organism.bed
#use bedtools to extract the best blast hit sequence from the proteomes.
bedtools getfasta -s -name -fi path/to/$Organism.pep -bed $Organism.bed -fo $Organism.target_loci.fasta
#add the organism name to the headers
sed -i "s/>/>${Organism}_/g" $Organism.target_loci.fasta
done
```
#### BUSCOs
```
#find all full length single copy BUSCOs in folder of busco_runs 
for file in $(find ./busco_runs -name "full_table_*.tsv"); do grep -v "^#" $file | awk '{if ($2 == "Complete") print $3 >> $1}'; done
mv EOG0* complete_buscos/
cd complete_buscos

#get fastas of genes for each BUSCO
for file in $(find . -name "EOG0*"); do perl extractFromFasta.pl allproteomes.fasta list $file > $file.pep; done

```
#### Orthofinder
```
#run full OrthoFinder analysis on FASTA format proteomes in <dir>
orthofinder -f proteomes_for_Orthofinder -M msa -T fasttree -S blast

```

### Align protein sequences with MAFFT:
```
mafft --localpair --maxiterate 1000 --thread 1 $file.pep > $file.pep.mafft
```
### Convert protein alignments to cdna alignments:
```
pal2nal.pl $file.pep.mafft $file.cdna -output fasta > $file.pep.mafft.p2n
```
### Refine alignments with trimAl
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
#### Format alignment files for coalescence and concatenation phylogenetic data sets 
```
#Strip gene IDs from species names in fasta header and append .fas suffix to file names 

for file in alignments/*.pep.mafft.p2n.trimal; do awk -F"|" '{print $1}' $file > $file.fas; done
```
#### generate supermatrix and corresponding partition files
```
perl ~/bin/FASconCAT-G_v1.02.pl -l -l -s 
```
#### RAxML with 1000 bootstrap replicates for individual alignments
```
for file in alignments/*.trimal.fas
do
raxmlHPC-PTHREADS-SSE3 -s $file -f a -m GTRGAMMA -p $RANDOM -# 1000 -T 10 -x $RANDOM -n $file
done
```
#### ASTRAL coalescence-based tree from individual gene trees
```
cat alignments/RAxML_bipartitions.* > cat.RAxML_bipartitions.genetrees #concatenate gene trees
astral -i cat.RAxML_bipartitions.genetrees -t 1 -o gene.Astral.tre
```
#### RAxML - partitioned supermatrix and 1000 bootstrap replicates 
```
raxmlHPC-PTHREADS-SSE3 -s alignment.supermatrix.fas -q alignment.supermatrix_partition.txt -f a -m GTRGAMMA -p $RANDOM -# 1000 -x 1234 -T 8 -n alignment.partition.supermatrix
```
#### MRBAYES - partitioned supermatrix 
```
mpiexec -np 1 mb alignment.supermatrix_partition.nex
```
#### IQ-Tree - 1000 bootstrap replicates for individual orthogroups alignments
```
for file in Orthogroup_alignments/*pep.mafft.p2n.trimal
do
iqtree2 -s $file -B 1000 -m TEST -mset raxml -T AUTO
done
```
#### Run ASTRAL-Pro on Orthogroup trees
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
### merge pairwise syntenic blocks into a multi-species ???blocks??? file
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
