## Trim and Filter Mar2017 2bRAD samples
First, extract and change name to get rid of unecessary parts of file name
```sh
for file in *.fastq.gz;
    do gunzip $file;
done
for file in *fastq;
    do mv $file ${file#CP-KS-2B-5-};
    done
for file in *fastq;
   do mv $file ${file/_L001_R1_001.fastq/_R1.fastq};
   done
```
Following previous notebook.
```sh
./2bRAD_trim_launch.pl fastq sampleID=1 site=.{12}GCA.{6}TGC.{12} > trims_Adults
```
Using fastx_toolkit for quality filtering. Submitted script with fastx_toolkit call for all samples.
```sh
module load fastx_toolkit
ls *.tr0 | perl -pe 's/^(\S+)\.tr0$/cat $1\.tr0 \| fastq_quality_filter -q 20 -p 90 >$1\.trim/' >filt0
#Edit for Q33
cat filt0 | perl -pe 's/filter /filter -Q33 /' > filt
#makes the following command for each sample
cat SS5-8-L5.tr0 | fastq_quality_filter -q 20 -p 90 >SS5-8-L5.trim
#Change name so has a fastq file extension
for file in *trim;
   do mv $file ${file/.trim/_trim.fastq};
   done
```
Edit filt file to make it suitable for submission to cluster.

Use cutadapt to trim off any adapter sequences left on reads and filter reads shorter than 34 bp:
```sh
module load cutadapt/1.2.1
for file in *.trim;
    do cutadapt -a AGATCG -m 34 --too-short-output ${file/_q0.fastq/toosh.fastq} -o ${file/_q0.fastq/ca.fastq} $file;
done
```
Would usually use the extractreads.py file (described below), to get the number of reads after all filtering, but I accidently deleted the output of cutadapt. So I used a modified file to get the number of filtered reads from the output of bowtie2 mapping.
```python
import re
import pandas as pd
import sys
from os import listdir
from os.path import isfile, join

def extractReads(indir,insuf,outfile):
    onlyfiles = [f for f in listdir(indir) if isfile(join(indir,f)) and insuf in f]
    sample_dict = {}
    reads_re = re.compile('\d+')
    #Extract number of reads from mapping output files
    for f in onlyfiles:
        IN = open(f,"r")
        name = f.translate(None,insuf)
        for line in IN:
            if "reads" in line:
                sample_dict[name] = {'Filtered_Reads':0}
                t_reads = reads_re.search(line).group()
                sample_dict[name]['Filtered_Reads'] = t_reads
        IN.close()
    #Write to new file using pandas
    df = pd.DataFrame(sample_dict).T
    df.to_csv(outfile,sep='\t', na_rep="none")
#get infile and outfile name arguments from command line
def main(argv):
    indir = argv[1]
    insuf = argv[2]
    outfile = argv[3]
    extractReads(indir,insuf, outfile)

if __name__ == "__main__":
    status = main(sys.argv)
    sys.exit(status)
```
Sample|	Filtered_Reads|Kept
---|---|---
HC1-14-L5|	1612570
HC1-2B-L5|	2131048
HC2-11-L5|	1713414
HC2-15A-L5|	1928858
HC2-5-L5|	1796542
HC2-6-L5|	2465684
HC3-14-L5|	1531433
HC3-15-L5|	984383
HC3-16-L5|	1480649
HC3-18-L5|	11217|n
HC3-5-L5|	1612829
HC4-17-L5|	1714278
HC4-19-L5|	1547948
HC5-2-L5|	1467271
HC5-3-L5|	1589960
HC5-9-L5|	1223051
NF1-7-L5|	1507120
NF2-1-L5|	1224054
NF2-10-L5|	1567789
NF2-11-L5|	1422167
NF2-19-L5|	1479307
NF2-7-L5|	1536566
NF3-12-L5|	1134519
NF3-13-L5|	1054943
NF3-14-L5|	1318132
NF3-2-L5|	1399583
NF3-7-L5|	1265670
NF4-18-L5|	1528887
NF5-1B-L5|	1361847
SS1-10-L5|	1489274
SS1-11-L5|	800990
SS1-17-L5|	2054627
SS1-18-L5|	968838
SS1-20-L5|	2198282
SS1-6-L5|	1614222
SS1-7-L5|	1814984
SS1-9-L5|	1647128
SS2-10-L5|	1632597
SS2-13-L5|	1674229
SS2-14-L5|	2026125
SS2-18-L5|	2064709
SS2-5B-L5|	1884417
SS2-9-L5|	1885530
SS3-14-L5|	1561863
SS3-18-L5|	5325|n
SS3-4-L5|	1447342
SS4-5-L5|	1083243
SS5-10-L5|	1930523
SS5-3-L5|	1469142
SS5-8-L5|	1788717

## Mapping Mar2017 samples
Did not prepare reference for mapping as was done previously, in case it affected the way that MBDSeq and 2bRAD mapping could be compared.
```sh
module load bowtie2/2.1.0 
module load picard-tools/1.92 
module load samtools/1.3.1

export GENOME_FASTA=wq
export GENOME_DICT= 
bowtie2-build $GENOME_FASTA $GENOME_FASTA
samtools faidx $GENOME_FASTA
java -jar $PICARD/CreateSequenceDictionary.jar R=$GENOME_FASTA O=$GENOME_DICT
```
Script to automate making job submission scripts for each sample to submit to tarbell cluster.
```sh
export REF=/scratch/t.cri.ksilliman/Ostrea_lurida-Scaff-10k.fa

for ffile in *.fastq; do 
    echo "#!/bin/bash"  >> ${ffile/ca.fastq/.sh}
    echo "#PBS -N ${ffile/ca.fastq/m}" >> ${ffile/ca.fastq/.sh}
    echo "#PBS -S /bin/bash" >> ${ffile/ca.fastq/.sh}
    echo "#PBS -l nodes=1:ppn=1"  >> ${ffile/ca.fastq/.sh}
    echo "#PBS -o ${ffile/ca.fastq/m.out}"  >>  ${ffile/ca.fastq/.sh}
    echo "#PBS -e ${ffile/ca.fastq/m.err}" >> ${ffile/ca.fastq/.sh}
    echo "#PBS -l walltime=120:00:00" >> ${ffile/ca.fastq/.sh}
    echo "#PBS -V" >> ${ffile/ca.fastq/.sh}
    echo "cd /scratch/t.cri.ksilliman/CommonG/2bRAD_Mar2017/ca_Mar2017/" >> ${ffile/ca.fastq/.sh}
    echo "module load bowtie2/2.1.0" >> ${ffile/ca.fastq/.sh}
    echo "bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x $REF -U $ffile -S ${ffile/ca.fastq/.bt2.sam}" >> ${ffile/ca.fastq/.sh}
done
#Submit to cluster
for file in *L5.sh; do qsub $file; done
```

Making .bam files from .sam files. NOTE: samtools updated since the last notebook, so sort had to be change from -f to -o.
```sh
for ffile in *.bt2.sam; do 
    echo "#!/bin/bash"  >> ${ffile/bt2.sam/s2b.sh}
    echo "#PBS -N ${ffile/bt2.sam/s2b}" >> ${ffile/bt2.sam/s2b.sh}
    echo "#PBS -S /bin/bash" >> ${ffile/bt2.sam/s2b.sh}
    echo "#PBS -l mem=5gb,nodes=1:ppn=1"  >> ${ffile/bt2.sam/s2b.sh}
    echo "#PBS -o ${ffile/bt2.sam/s2b.out}"  >>  ${ffile/bt2.sam/s2b.sh}
    echo "#PBS -e ${ffile/bt2.sam/s2b.err}" >> ${ffile/bt2.sam/s2b.sh}
    echo "#PBS -l walltime=50:00:00" >> ${ffile/bt2.sam/s2b.sh}
    echo "cd /scratch/t.cri.ksilliman/CommonG/2bRAD_Mar2017/map_Mar2017/" >> ${ffile/bt2.sam/s2b.sh}
    echo "module load samtools/1.3.1" >> ${ffile/bt2.sam/s2b.sh} 
    echo "module load picard-tools/1.92" >> ${ffile/bt2.sam/s2b.sh} 
    echo "samtools view -bS $ffile > ${ffile/.bt2.sam/_unsorted.bam}" >> ${ffile/bt2.sam/s2b.sh}
    echo "samtools sort ${ffile/.bt2.sam/_unsorted.bam} -o ${ffile/.bt2.sam/_sorted.bam}" >> ${ffile/bt2.sam/s2b.sh}
    echo "java -Xmx5g -jar \$PICARD/AddOrReplaceReadGroups.jar INPUT=${ffile/.bt2.sam/_sorted.bam} OUTPUT=${ffile/bt2.sam/bam} RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${ffile/bt2.sam/m}" >> ${ffile/bt2.sam/s2b.sh} 
    echo "samtools index ${ffile/bt2.sam/bam}" >> ${ffile/bt2.sam/s2b.sh}
done

#Submit scripts
for file in *s2b.sh; do qsub $file; done
rm *sorted*
```
## GATK
Realigning around indels
```sh
ls *.bams > bams
cat bams | perl -pe 's/(\S+)\.bam/java -Xmx5g -jar \$GATK -T RealignerTargetCreator -R \$GENOME_REF -I $1\.bam -o $1\.intervals/' >intervals

#edit intervals so it can be submitted to cluster 
qsub intervals. sh

#Realigning
cat bams | perl -pe 's/(\S+)\.bam/java -Xmx5g -jar \$GATK -T IndelRealigner -R \$GENOME_REF -targetIntervals $1\.intervals -I $1\.bam -o $1\.real.bam -LOD 0\.4/' >realign

qsub realign.sh
```
# Assembly and Genotyping of "best 90" Samples
Combined samples with the most post-filtering reads from Dec 2015 and Mar 2017 sequencing runs, as well as all successfully sequenced replicates and samples used for MBD-Seq. Have 106 samples total (16 replicates).

Run UnifiedGenotyper for a first pass in order to get base quality recalibration. Needed 20gb of RAM for job to run.
```sh
#PBS -l nodes=1:ppn=4
#PBS -l mem=20gb
#PBS -o round1b90.out
#PBS -e round1b90.err
#PBS -l walltime=50:00:00
cd $PATH/best90_map/
module load gatk/3.6

export GENOME_REF=/scratch/t.cri.ksilliman/Ostrea_lurida-Scaff-10k.fa

java -jar $GATK -T UnifiedGenotyper -R $GENOME_REF -nt 4 -nct 1 \
-I HC1-14-L5.real.bam \ 
ect...
-o round1best90.vcf

python ../z0on_2brad/GetHighQualVcfs.py -i round1.vcf --percentile 75 -o .
```
Recalibrate base quality scores
```sh
ls *real.bams >bams
cat bams | perl -pe 's/(\S+)\.real\.bam/java -Xmx20g -jar \$GATK -T BaseRecalibrator -R \$GENOME_REF -knownSites $1.m_HQ\.vcf -I $1\.real\.bam -o $1\.real\.recalibration_report.grp/' >bqsr
#Edit to submit to cluster
```
Rewrite bams according to recalibration reports
```sh
cat bams | perl -pe 's/(\S+)\.bam/java -Xmx10g -jar \$GATK -T PrintReads -R \$GENOME_REF -I $1\.bam -BQSR $1\.recalibration_report.grp -o $1\.recal\.bam /' >bqsr2
#Example:
java -Xmx20g -jar $GATK -T BaseRecalibrator -R $GENOME_REF -knownSites SS5-8-L5.m_HQ.vcf -I SS5-8-L5.real.bam -o SS5-8-L5.real.recalibration_report.grp
```
2nd iteration of UnifiedGenotyper on quality-recalibrated files:
```sh
 ls *.recal.bam > bams
 cat bams | perl -pe 's/(\S+\.bam)/-I $1 \\/' >> unig2.sh
#edit for submissing to cluster.
#PBS -N round2b90
#PBS -l nodes=1:ppn=8
#PBS -l mem=20gb
#PBS -o round2b90.out
#PBS -e round2b90.err
#PBS -l walltime=100:00:00

module load gatk/3.6
export GENOME_REF=/scratch/t.cri.ksilliman/Ostrea_lurida-Scaff-10k.fa

java -jar $GATK -T UnifiedGenotyper -R $GENOME_REF -nt 8 -nct 1 \
-I HC1-11.real.recal.bam \
#etc.
-o round2b90.vcf
```
Edit .vcf to get rid of .m in file names. Remove individuals with fewer than 500k reads after all filtering:
```sh
cat round2.vcf | perl -pe 's/\.fq\.trim\.bt2//g' | perl -pe 's/\.trim\.bt2//g' | perl -pe 's/\.m//g' | perl -pe 's/^chrom/chr/' >round2.names.vcf
```
----------------------------
### Variant Quality Score Recalibration (VQSR)
Make .tab file with replicates, best90_reps.tab:
```
HC1-2B-L5B	HC1-2B-L5
HC2-15A-L5	HC2-15-L5B
HC2-17A	HC2-17B
HC3-5-L5B HC3-5-L5
HC4-4A	HC4-4B
NF1-14A	NF1-14B
NF2-6A	NF2-6B
NF2-6B	NF2-6C
NF5-3B	NF5-3
SS2-12A	SS2-12B
SS2-4A	SS2-4B
SS3-5B	SS3-5C
SS3-5C	SS3-5D
SS4-10A	SS4-10B
```
extracting SNPs that are consistently genotyped in replicates 
and have not too many heterozygotes:
```sh
../z0on_2brad/replicatesMatch.pl vcf=round2b90.names.vcf replicates=b90_reps.tab > vqsr.vcf
```
Output:
```
8423 total SNPs
1242 pass hets and match filters
756 show non-reference alleles
615 have alterantive alleles in at least 2 replicate pair(s)
615 have matching heterozygotes in at least 0 replicate pair(s)
263 polymorphic
615 written
```
Determining transition-transversion ratio for true snps (will need it for tranche calibration)
```sh
vcftools --vcf vqsr.vcf --TsTv-summary
```
Output: Ts/Tv ratio: 1.449
Recalibrating genotype calls: VQSR
 * step one - creating recalibration models
```sh
#PBS -l nodes=1:ppn=8
#PBS -o vqsr.out
#PBS -e vqsr.err
cd /scratch/t.cri.ksilliman/CommonG/best90_map/
module load gatk/3.6

export GENOME_REF=/scratch/t.cri.ksilliman/Ostrea_lurida-Scaff-10k.fa

java -jar $GATK -T VariantRecalibrator -R $GENOME_REF -input round2b90.names.vcf \
-nt 8 -resource:repmatch,known=true,training=true,truth=true,prior=10 \
vqsr.vcf -an QD -an MQ -an FS -mode SNP --maxGaussians 4 --target_titv 1.449 \
-tranche 85.0 -tranche 90.0 -tranche 95.0 -tranche 99.0 -tranche 100 -recalFile\
round2b90.recal -tranchesFile recalibrate_SNP.tranche
```
non-parametric quantile-based recalibration a-la de novo pipeline
```
$ ../z0on_2brad/replicatesMatch.pl vcf=round2b90.names.vcf replicates=b90_reps.tab polyonly=1 >vqsr.poly.vcf

8423 total SNPs
1207 pass hets and match filters
740 show non-reference alleles
607 have alterantive alleles in at least 2 replicate pair(s)
607 have matching heterozygotes in at least 0 replicate pair(s)
266 polymorphic
266 written

$ ../z0on_2brad/recalibrateSNPs_gatk.pl vcf=round2b90.names.vcf true=vqsr.poly.vcf > gatk_after_vqsrNP.vcf
------------------------
65.19%	at qual <1 (64.19% gain)
70.53%	at qual <5 (65.53% gain)
76.84%	at qual <10 (66.84% gain)
80.13%	at qual <15 (65.13% gain)
82.07%	at qual <20 (62.07% gain)
85.58%	at qual <30 (55.58% gain)
------------------------
```
qual less than 10 had highest gain

Filtering SNPs with excess heterozygotes using dDocent lab's script, filter_hetexc_by_pop.pl
```sh
$  ../filter_hetexc_by_pop.pl -v gatk_after_vqsrNP.vcf -p three.pop -h 0.001 -c 0.1 -o hetfiltNP_3
Processing population: HC (35 inds)
Processing population: NF (34 inds)
Processing population: SS (37 inds)
Outputting results of HWE test for filtered loci to 'filtered.hetexc'
Kept 8199 of a possible 8423 loci (filtered 224 loci)

$ vcftools --vcf hetfiltNP_3.recode.vcf --max-missing 0.5 --recode --recode-INFO-all --out hetfiltNP_3_m50

After filtering, kept 106 out of 106 Individuals
After filtering, kept 3627 out of a possible 8199 Sites

#Using z0on script
$ ../z0on_2brad/hetfilter.pl vcf=gatk_after_vqsrNP.vcf > hetfiltNP_z.vcf

8423 total loci
4572 dropped because fraction of missing genotypes exceeded 0.5
152 dropped because fraction of heterozygotes exceeded 0.75
3699 written
```
Thinning dataset to 1 SNP. Did not do hetfilter
```sh
../z0on_2brad/thinner.pl vcf=hetfiltNP_3_m50.recode.vcf > thinNP_m50.vcf

3627 total loci
233 loci skipped because they were closer than 40
2044 loci selected
```
Thinning for tajima's D down the road
```sh
$ vcftools --vcf hetfiltNP_3_m50.recode.vcf --remove clones2remove.all --minQ 10 --max-missing 0.75  --maf 0.025 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out b90finalmaf025_Unthin

After filtering, kept 91 out of 106 Individuals
Outputting VCF file...
After filtering, kept 1141 out of a possible 3627 Sites
Run Time = 2.00 seconds
#No maf filtering
After filtering, kept 1600 out of a possible 3627 Sites
```
Filtering based on quality from non parametric recalibration. Need to fill in minQ with value from above.
```sh
$ vcftools --vcf thinNP_m50.vcf --minQ 10 --max-missing 0.75  --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out filt0NP

After filtering, kept 1008 out of a possible 2045 Sites
```
genotypic match between pairs of replicates 
```sh
$ ../z0on_2brad/repMatchStats.pl vcf=filt0NP.recode.vcf replicates=b90_reps.tab 
```
pair |	gtyped	|match	|[ 00	01	11 ]	HetMatch	HomoHetMismatch	HetNoCall	|HetsDiscoveryRate
---|---|---|-----|---
HC1-2B-L5B:HC1-2B-L5|	1193|	1155|(96.8%)	 [60%	15%	25% ]	169	33	|	0.91	
HC2-15A-L5:HC2-15-L5B|	1242|	867|(69.8%)	 [60%	15%	25% ]	129	320	8|	0.44	
HC2-17A:HC2-17B	|1205	|1151|(95.5%)	 [59%	14%	26% ]	164	43	1	|0.88	
HC3-5-L5B:HC3-5-L5|	1184|	1141|(96.4%)	 [59%	16%	25% ]	186	35	|	0.91	
HC4-4A:HC4-4B	|1212|	1161|(95.8%)	 [60%	14%	26% ]	163	39	2	|0.89	
NF1-14A:NF1-14B	|1211|	1178|(97.3%)	 [58%	15%	26% ]	181	28		|0.93	
NF2-6A:NF2-6B	|1190|	1141|(95.9%)	 [56%	18%	25% ]	209	36		|0.92	
NF2-6B:NF2-6C	|1192|	973|(81.6%)	 [59%	15%	26% ]	146	62	24	|0.77	
NF5-3B:NF5-3	|1211|	1160|(95.8%)	 [60%	14%	26% ]	162	32	2|	0.91	
SS2-12A:SS2-12B	|1234|	966|(78.3%)	 [62%	11%	26% ]	110	225	|	0.49	
SS2-4A:SS2-4B	|1215|	1167|(96.0%)	 [55%	21%	24% ]	244	47	2|	0.91	
SS3-5B:SS3-5C	|1196|	1122|(93.8%)	 [60%	13%	27% ]	149	44	5|	0.86	
SS3-5C:SS3-5D	|1194|	1133|(94.9%)	 [59%	14%	26% ]	162	31	3|	0.91	
SS4-10A:SS4-10B	|1192|	1079|(90.5%)	 [62%	12%	26% ]	131	46	5|	0.84	

All of these look pretty good except HC2-15-L5B:HC2-15-L5 and SS2-12A:SS2-12B. I'm going to leave these replicates in for now, and then plot with PCA to see if they look weird.

Edit: checked PCA and they still cluster correctly. I think it is due to these samples havingf a much greater difference in missing data than others. EWill remove the copes at the end for SS2-12B and HC2-15-L5B.

Make clones2remove.some file with replicates to remove
```HC1-2B-L5B
HC2-17A
HC3-5-L5B
HC4-4B
NF1-14A
NF2-6A
NF2-6C
NF5-3
SS2-4A
SS3-5C
SS3-5D
SS4-10B
```
Final file b90finalNP-m75.recode.vcf
```
$ vcftools --vcf thinNP_m50.vcf --remove clones2remove.all --minQ 10 --max-missing 0.75  --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out b90finalNP-m75

Excluding individuals in 'exclude' list
After filtering, kept 91 out of 106 Individuals
Outputting VCF file...
After filtering, kept 945 out of a possible 2044 Sites
```
After filtering, kept 91 out of 106 Individuals
Outputting VCF file...
After filtering, kept 711 out of a possible 2044 Sites
Run Time = 2.00 seconds

## Analysis
Running bayescan

WC FST
```
$ vcftools --vcf ../Inputs/b90finalNP-m75maf025.recode.vcf --weir-fst-pop ../SS.pop --weir-fst-pop ../NF.pop --weir-fst-pop ../HC.pop

Outputting Weir and Cockerham Fst estimates.
Weir and Cockerham mean Fst estimate: 0.056388
Weir and Cockerham weighted Fst estimate: 0.063903
```
MAking Structure
```sh
fixing chr names
 awk '{gsub(/^chr/,""); print}' chromb90maf025.recode.vcf > nochrb90maf025.vcf

$ ../z0on_2brad/thinner.pl vcf=hetfiltNP_3_m50.recode.vcf interval=5000 | perl -pe 's/scaffold/chr/g' >thinchromb90.vcf

3627 total loci
56 loci skipped because they were closer than 5000
1863 loci selected
$ vcftools --vcf thinchromb90.vcf --remove clones2remove.all --minQ 10 --max-missing 0.75  --min-alleles 2 --max-alleles 2 --maf 0.025 --recode --out chromb90maf025
After filtering, kept 91 out of 106 Individuals
Outputting VCF file...
After filtering, kept 662 out of a possible 1863 Sites
#Non-maf filtering
After filtering, kept 879 out of a possible 1863 Sites

# reformatting VCF into plink binary BED format
$ plink --vcf chromb90maf025.recode.vcf --make-bed --allow-extra-chr --out b90maf025.admix

for K in 1 2 3 4 5 6; do ../admixture/admixture --cv b90maf025.admix.bed $K | tee log$K.out; done

grep -h CV log*.out
CV error (K=1): 0.51591
CV error (K=2): 0.49729
CV error (K=3): 0.47835
CV error (K=4): 0.47789
CV error (K=5): 0.47185
CV error (K=6): 0.48451

#non-maf filtering
CV error (K=1): 0.41895
CV error (K=2): 0.40288
CV error (K=3): 0.39039
CV error (K=4): 0.38908
CV error (K=5): 0.38560
CV error (K=6): 0.39963

```

## Moving Dec2015 2bRAD over to tarbell  

Extracting the total number of reads and the filtered number of reads from cutadapt's output and putting it in 1 tab-separated file. 
-bash: BASH_FUNC_module(): line 0: syntax error near unexpected token `)'
-bash: BASH_FUNC_module(): line 0: `BASH_FUNC_module() () {  eval `/usr/bin/modulecmd bash $*`'
-bash: error importing function definition for `BASH_FUNC_module'
1612570 reads; of these:
  1612570 (100.00%) were unpaired; of these:
    1052638 (65.28%) aligned 0 times
    303072 (18.79%) aligned exactly 1 time
    256860 (15.93%) aligned >1 times
34.72% overall alignment rate
```python
import re
import locale
import pandas as pd
import sys

def extractReads(infile,outfile):    
    locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
    sample_dict = {}
    IN = open(infile,"r")
    #Extract number of reads from cut_adapt output file 
    name_re = re.compile('(\w+-\d+\D*)(_)')
    reads_re = re.compile('\d+,*\d*,*\d*')
    for line in IN:
        if "parameters" in line:
            mO = re.search(name_re,line)
            name = mO.group(1)
            sample_dict[name] = {'Total_Reads':0,'Filtered_Reads':0}
        if "Total reads processed" in line:
            t_reads = locale.atoi(reads_re.search(line).group())
            sample_dict[name]['Total_Reads'] = t_reads
        if "Reads written" in line:
            w_reads = locale.atoi(reads_re.search(line).group())
            sample_dict[name]['Filtered_Reads'] = w_reads
    IN.close()
    #Write to new file using pandas
    order = ['Total_Reads','Filtered_Reads']
    df = pd.DataFrame(sample_dict).T.reindex(columns=order)
    df.to_csv(outfile,sep='\t', na_rep="none")
#get infile and outfile name arguments from command line
def main(argv):
    infile = argv[1]
    outfile = argv[2]
    extractReads(infile, outfile)

if __name__ == "__main__":
    status = main(sys.argv)
    sys.exit(status)
```
Did for ca_2brad_201512.out and saved in ca_2b_201512.summary. Moved over samples with at least 500,000 filtered reads to tarbell. Note: The last cutadapt output for the last file, SS5-9, was cut off, so needed to use grep to get the number of filtered reads. Edited the ca_2b_201512.summary by hand for that sample.
Sample | Total_Reads | Filtered_Reads | Mapped
---|---|---|---|
HC1-1|1268772|1263491|y
HC1-10|1225170|1219587| y
HC1-11|2535757|2526024| y
HC1-12	|1717040|	1709933| y
HC1-13	|3828645|	3812983| y
HC1-14	|2291	|2283| n
HC1-15	|2737811|	2726748| y
HC1-17	|5293446|	5272078| y
HC1-18	|2680752|	2668705| y
HC1-19	|1984978|	1974182| y
HC1-2	|1010	|1006| n
HC1-3	|1511800|	1504925| y
HC1-4	|1546873|	1540114| y
HC1-9	|1732814|	1725537| y
HC2-1	|2645537|	2634778| y
HC2-10	|2693857|	2683936| y
HC2-11	|1657	|1649| n
HC2-12	|2885460|	2872374| y
HC2-13	|2059669|	2050077| y
HC2-14	|1942224|	1933549| y
HC2-15	|1640	|1632| n
HC2-16	|1996899|	1987311| y
HC2-17A	|1352044|	1345996| y
HC2-17B	|2799603|	2788009| y
HC2-17C	|1588	|1581| n
HC2-17D	|79	|79 | n
HC2-18	|1602743|	1596093| y
HC2-19	|4081416|	4064304| y
HC2-2	|3113264|	3101304| y
HC2-20	|1544820|	1537734| y
HC2-3	|1936742|	1929163| y
HC2-4	|2783094|	2773083| y
HC2-5	|59|	59| n
HC2-6	|2857|	2851| n
HC2-7	|5072537|5054136| y
HC2-8	|1035482|	1031698|y
HC2-9	|3078436|	3067001|y
HC3-1	|1218091|	1212847|y
HC3-10	|2683090|	2670752|y
HC3-11	|1908913|	1901118|y
HC3-12	|1683125|	1675065|y
HC3-13	|1623968|	1617312|y
HC3-14	|2077291|	2067995|y
HC3-15	|20815|	20711|n
HC3-16	|1633372|	1626123|y
HC3-17	|1196412|	1191443|y
HC3-18	|1676072|	1668818|y
HC3-2	|1568096|	1561163|y
HC3-3	|1667571|	1660960|y
HC3-4	|2494888|	2484001|y
HC3-5	|1215|	1208|n
HC3-6	|1541293|	1534810|y
HC3-7	|1562722|	1556367|y
HC3-8	|1693813|	1686649|y
HC3-9	|2389283|	2379381|y
HC4-10	|983269 |979069|y
HC4-11	|1169618|	1164870|y
HC4-12A	|857374| 853423|y
HC4-13A |1135990|	1131538|y
HC4-14	|1273198|	1267576|y
HC4-15	|915148	|911175|y
HC4-16	|1187286|	1182640|y
HC4-17	|1482	|1476|n
HC4-18	|1091726|	1087118|y
HC4-19	|1283	|1277|n
HC4-1A	|5315133|	5293397|y
HC4-2	|4586311|	4567570|y
HC4-20	|1087970|	1083213|y
HC4-3	|6872868|	6843355|y
HC4-4A	|2988443|	2977343|y
HC4-4B	|2337983|	2328610|y
HC4-5	|4304093|	4289068|y
HC4-6	|1752516|	1744669|y
HC4-7	|1168808|	1163684|y
HC4-8	|1343076|	1337453|y
HC4-9	|1290278|	1284837|y
HC5-1	|1305084|	1298808|y
HC5-10	|1467710|	1460071|y
HC5-11	|2273381|	2262220|y
HC5-12	|2427463|	2416220|y
HC5-13	|3298151|	3283166|y
HC5-14	|0	|0|y
HC5-15	|2130537|	2120787|y
HC5-16	|2246303|	2236585|y
HC5-2	|1632	|1627|n
HC5-3	|1200	|1196|n
HC5-4	|1030931|	1026097|y
HC5-5	|1331224|	1324981|y
HC5-6	|989152	|984663|y
HC5-7	|644228	|640837|y
HC5-8	|1350680|	1344491|y
HC5-9|	2524	|2517|n
NF1-1|	3027613	|3013832|y
NF1-10|	1812208	|1805680|y
NF1-11|	1963447	|1955687|y
NF1-12|	2182057	|2173776|y
NF1-13|	1992326	|1985121|y
NF1-14A|	2066735|	2059315|y
NF1-14B	|1809448	|1802475|y
NF1-15|	1903327	|1895668|y
NF1-16|	2091937	|2083350|y
NF1-17|	2254489	|2245290|y
NF1-18|	2464202|	2455372|y
NF1-19|	2099700|	2089552|y
NF1-2|	1369101	|1363753|y
NF1-20|	2539513	|2530153|y
NF1-21|	1920	|1910|n
NF1-3|	1596220	|1588315|y
NF1-7|	15849	|15790|n
NF1-8|	1838791	|1829916|y
NF1-9|	2193851	|2185759|y
NF2-1|	2031	|2024|n
NF2-10-SS|	2328|	2319|n
NF2-11-SS|	1608379|	1601165|n
NF2-12|	1886944	|1877302|y
NF2-13|	2132490|	2122566|y
NF2-15|	2431649	|2418588|y
NF2-16|	478723|	475988|n
NF2-17|	1549998|	1543662|y
NF2-18A|	411492|	409674|n
NF2-18B	|4375|	4360|n
NF2-19|	2098862|	2088477|y
NF2-2|	1377623|	1370978|y
NF2-20|	4523946|	4504573|y
NF2-3|	1425827|	1419477|y
NF2-4|	1453047|	1447073|y
NF2-5|	1769535	|1761232|y
NF2-6A|	1683099	|1674902|y
NF2-6B|	1176575	|1171530|y
NF2-6C|	954682	|950508|y
NF2-6D|	9831	|9781|n
NF2-9|	764124	|760032|y
NF3-1|	1787934	|1779156|y
NF3-10|	1994258	|1986821|y
NF3-11|	1917619	|1908608|y
NF3-12|	3704	|3689|n
NF3-13|	178939	|178199|n
NF3-14|	1439	|1431|n
NF3-15|	946707	|942561|y
NF3-16|	483527	|481631|n
NF3-17|	1264738	|1260094|y
NF3-18|	1082513	|1077934|y
NF3-19|	941167	|937583|y
NF3-2|	2505	|2494|n
NF3-20|	1121101	|1115618|y
NF3-3|	2842249	|2826903|y
NF3-4|	1634	|1629|n
NF3-5|	1453430	|1446748|y
NF3-6|	1533860	|1528338|y
NF3-7|	1574	|1565|n
NF3-8|	2230249	|2221723|y
NF3-9|	1990095	|1982165|y
NF4-1|	1236208	|1230986|y
NF4-10|	1010743	|1006044|y
NF4-11|	807395	|804316|y
NF4-12|	1176631	|1172212|y
NF4-13|	1060954	|1056214|y
NF4-14|	1461041	|1452129|y
NF4-15|	2029	|2018|n
NF4-16|	1331687	|1326299|y
NF4-17|	1130413	|1126313|y
NF4-18|	1874	|1863|n
NF4-19|	2533773	|2524293|y
NF4-2|	1147317	|1142958|y
NF4-3|	1150847	|1145819|y
NF4-4|	1151562	|1147156|y
NF4-5|	1118405	|1113683|y
NF4-6|	1070023	|1065253|y
NF4-7|	1100397	|1094684|y
NF4-8|	888449	|885167|y
NF4-9	618454	|615907|y
NF5-1|	3387	|3362|n
NF5-10|	1668229	|1661457|y
NF5-11|	1642580	|1635259|y
NF5-12|	2712807	|2700206|y
NF5-13|	2000767	|1992004|y
NF5-14|	2451966	|2441033|y
NF5-15|	2225049	|2214389|y
NF5-16|	1404	|1396|n
NF5-17|	889	|884|n
NF5-18|	1780626|	1772296|y
NF5-19|	2408893	|2398621|y
NF5-3|	2175715	|2166392|y
NF5-3B|	2738455	|2726728|y
NF5-4|	1663698	|1655075|y
NF5-5|	2025689	|2015343|y
NF5-6|	1971149	|1962788|y
NF5-7|	1431173	|1422791|y
NF5-8|	2212107	|2202355|y
NF5-9|	2201317	|2191586|y
SS1-1|	1611195	|1604300|y
SS1-10|	989	|989|n
SS1-11|	3656	|3635|n
SS1-12|	1523234|	1516558|y
SS1-13|	781287	|778245|y
SS1-14|	1362147	|1356326|y
SS1-15|	498974	|496864|n
SS1-16|	642006	|639322|y
SS1-17|	1240|	1236|n
SS1-18|	847078|	843506|y
SS1-19|	1549475|	1543188|y
SS1-2|	1627718	|1621168|y
SS1-20|	2712	|2702|n
SS1-3|	1249866|	1244449|y
SS1-4|	899499	|895883|y
SS1-5|	1628288	|1622061|y
SS1-6|	23	|23|n
SS1-7|	1178|	1173|n
SS1-8|	1835590|	1827356|y
SS1-9|	1459	|1452|n
SS2-1|	1272343	|1266731|y
SS2-10|	614	|608|n
SS2-11|	64	|62|n
SS2-12A|	2605101|	2593577|y
SS2-12B	|1481306	|1474501|y
SS2-13|	800867	|797506|y
SS2-13-NF|	1211157|	1206182|y
SS2-14|	19557	|19445|n
SS2-15|	2586316	|2576340|y
SS2-16|	2073918	|2067797|y
SS2-17|	2780261|	2768720|y
SS2-18|	2844	|2836|y
SS2-19|	2507913|	2496969|y
SS2-2|	685945|	682963|y
SS2-3|	884632|	881024|y
SS2-4A|	905029|	900941|y
SS2-4B|	3206755|	3192141|y
SS2-5|	3392	|3377|n
SS2-9|	501	|499|n
SS3-1|	2968054|	2957839|y
SS3-10|	2378558	|2368646|y
SS3-11|	2562879	|2553796|y
SS3-12|	1292889|	1287155|y
SS3-13|	1953036|	1945198|y
SS3-14|	2466006	|2455908|y
SS3-15|	2570524|	2559903|y
SS3-16|	1914066|	1905649|y
SS3-17|	1867838|	1859857|y
SS3-18|	1780|	1769|n
SS3-19|	2466897|	2457043|y
SS3-2|	2623623	|2613687|y
SS3-20|	766362	|763028|y
SS3-21|	956152	|952393|y
SS3-3|	1961959	|1953724|y
SS3-4|	2169	|2161|n
SS3-5A|	1561	|1553|n
SS3-5B|	3223513	|3209113|y
SS3-5C|	1264471	|1259090|y
SS3-5D|	1450307	|1443561|y
SS3-6|	2566944	|2556481|y
SS3-7|	1902180	|1894646|y
SS3-8|	2526682	|2516469|y
SS3-9|	2149285	|2140117|y
SS4-1|	1380260	|1374271|y
SS4-10A|	982128	|978009|y
SS4-10B	|559899	|557901|y
SS4-11|	1062737	|1058922|y
SS4-12|	1159439	|1154949|y
SS4-13|	1138140	|1133306|y
SS4-14|	1204901	|1199619|y
SS4-15|	1059439	|1055062|y
SS4-16|	1167428	|1162631|y
SS4-17|	1111396	|1106599|y
SS4-18|	1036292	|1032026|y
SS4-19|	937099	|933106|y
SS4-2|	550769	|548511|y
SS4-20|	628341	|625439|y
SS4-3|	1045200	|1041082|y
SS4-4|	1358670	|1353327|y
SS4-5|	1175	|1167|n
SS4-6|	1050914	|1046786|y
SS4-7|	1164666	|1160532|y
SS4-8|	1292049	|1286555|y
SS4-9|	753177	|750050|y
SS5-1|	1030929	|1025560|y
SS5-10|	970	|967|n
SS5-11|	1286812	|1280319|y
SS5-12|	1312319	|1305974|y
SS5-13A|	790865|	787174|y
SS5-13B	|4359006|	4338872|y
SS5-14	|1037668|	1032836|y
SS5-15|	2432902	|2421539|y
SS5-16|	502696	|500103|y
SS5-17|	1187948	|1182780|y
SS5-18|	1550948	|1543523|y
SS5-19|	1475170	|1468186|y
SS5-2|	358620	|356384|n
SS5-3|	4652	|4621|n
SS5-4|	1784333	|1775193|y
SS5-5|	1987909	|1979082|y
SS5-6|	1981802	|1973400|y
SS5-7|	2397777	|2386765|y
SS5-8|	108	|108|n
SS5-9|	0	|1255239|y

## Quality filter using fastx
```sh
#PBS -N qualDec15
#PBS -o qualDec15.out
#PBS -e qualDec15.err
module load fastx_toolkit/0.0.13 
cd /scratch/t.cri.ksilliman/CommonG/2bRAD_Dec2015/trim_Dec2015/

cat HC1-10.tr0 | fastq_quality_filter -Q33 -q 20 -p 90 >HC1-10.trim
#ect...

#Change name so samples have a fastq file extension
for file in *trim;
   do mv $file ${file/.trim/_trim.fastq};
   done
```
Use cutadapt to trim off any adapter sequences left on reads and filter reads shorter than 34 bp:
```sh
module load cutadapt/1.2.1
for file in *.trim;
    do cutadapt -a AGATCG -m 34 --too-short-output ${file/_trim.fastq/toosh.fastq} -o ${file/_trim.fastq/ca.fastq} $file;
done
```
## Map to genome with bowtie2
```sh
export REF=/scratch/t.cri.ksilliman/Ostrea_lurida-Scaff-10k.fa

for ffile in *.fastq; do 
    echo "#!/bin/bash"  >> ${ffile/ca.fastq/.sh}
    echo "#PBS -N ${ffile/ca.fastq/m}" >> ${ffile/ca.fastq/.sh}
    echo "#PBS -S /bin/bash" >> ${ffile/ca.fastq/.sh}
    echo "#PBS -l nodes=1:ppn=1"  >> ${ffile/ca.fastq/.sh}
    echo "#PBS -o ${ffile/ca.fastq/m.out}"  >>  ${ffile/ca.fastq/.sh}
    echo "#PBS -e ${ffile/ca.fastq/m.err}" >> ${ffile/ca.fastq/.sh}
    echo "#PBS -l walltime=120:00:00" >> ${ffile/ca.fastq/.sh}
    echo "#PBS -V" >> ${ffile/ca.fastq/.sh}
    echo "cd /scratch/t.cri.ksilliman/CommonG/2bRAD_Mar2017/ca_Mar2017/" >> ${ffile/ca.fastq/.sh}
    echo "module load bowtie2/2.1.0" >> ${ffile/ca.fastq/.sh}
    echo "bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x $REF -U $ffile -S ${ffile/ca.fastq/.bt2.sam}" >> ${ffile/ca.fastq/.sh}
done
#Submit to cluster
for file in *L5.sh; do qsub $file; done
mv *sm ../map_Dec2015
```
Convert from sam to bam
```sh
for ffile in *.bt2.sam; do 
    echo "#!/bin/bash"  >> ${ffile/bt2.sam/s2b.sh}
    echo "#PBS -N ${ffile/bt2.sam/s2b}" >> ${ffile/bt2.sam/s2b.sh}
    echo "#PBS -S /bin/bash" >> ${ffile/bt2.sam/s2b.sh}
    echo "#PBS -l mem=5gb,nodes=1:ppn=1"  >> ${ffile/bt2.sam/s2b.sh}
    echo "#PBS -o ${ffile/bt2.sam/s2b.out}"  >>  ${ffile/bt2.sam/s2b.sh}
    echo "#PBS -e ${ffile/bt2.sam/s2b.err}" >> ${ffile/bt2.sam/s2b.sh}
    echo "#PBS -l walltime=50:00:00" >> ${ffile/bt2.sam/s2b.sh}
    echo "cd /scratch/t.cri.ksilliman/CommonG/2bRAD_Mar2017/map_Mar2017/" >> ${ffile/bt2.sam/s2b.sh}
    echo "module load samtools/1.3.1" >> ${ffile/bt2.sam/s2b.sh} 
    echo "module load picard-tools/1.92" >> ${ffile/bt2.sam/s2b.sh} 
    echo "samtools view -bS $ffile > ${ffile/.bt2.sam/_unsorted.bam}" >> ${ffile/bt2.sam/s2b.sh}
    echo "samtools sort ${ffile/.bt2.sam/_unsorted.bam} -o ${ffile/.bt2.sam/_sorted.bam}" >> ${ffile/bt2.sam/s2b.sh}
    echo "java -Xmx5g -jar \$PICARD/AddOrReplaceReadGroups.jar INPUT=${ffile/.bt2.sam/_sorted.bam} OUTPUT=${ffile/bt2.sam/bam} RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${ffile/bt2.sam/m}" >> ${ffile/bt2.sam/s2b.sh} 
    echo "samtools index ${ffile/bt2.sam/bam}" >> ${ffile/bt2.sam/s2b.sh}
done

#Submit scripts
for file in *s2b.sh; do qsub $file; done
rm *sorted*
```
## GATK
Realigning around indels
```sh
#finding places to realign
ls *.bam > bams
cat bams | perl -pe 's/(\S+)\.bam/java -Xmx5g -jar \$GATK -T RealignerTargetCreator -R \$GENOME_REF -I $1\.bam -o $1\.intervals/' >intervals

#edit intervals so it can be submitted to cluster 
qsub intervals. sh
#realigning
cat bams | perl -pe 's/(\S+)\.bam/java -Xmx5g -jar \$GATK -T IndelRealigner -R \$GENOME_REF -targetIntervals $1\.intervals -I $1\.bam -o $1\.real.bam -LOD 0\.4/' >realign.sh
#edit to submit to cluster

```
1st round of GATK
```sh
```
March, remove less than 500k reads
```sh
vcftools --vcf round2.names.vcf --remove-indv HC3-18-L5 --remove-indv SS3-18-L5 --recode --recode-INFO-all --out Mar.round2.over5k
mv Mar.round2.over5k.recode.vcf Mar.round2.over5k.vcf
```

  - Import a HTML file and watch it magically convert to Markdown
  - Drag and drop images (requires your Dropbox account be linked)


You can also:
  - Import and save files from GitHub, Dropbox, Google Drive and One Drive
  - Drag and drop markdown and HTML files into Dillinger
  - Export documents as Markdown, HTML and PDF

Markdown is a lightweight markup language based on the formatting conventions that people naturally use in email.  As [John Gruber] writes on the [Markdown site][df1]

> The overriding design goal for Markdown's
> formatting syntax is to make it as readable
> as possible. The idea is that a
> Markdown-formatted document should be
> publishable as-is, as plain text, without
> looking like it's been marked up with tags
> or formatting instructions.

This text you see here is *actually* written in Markdown! To get a feel for Markdown's syntax, type some text into the left window and watch the results in the right.

### Tech

Dillinger uses a number of open source projects to work properly:

* [AngularJS] - HTML enhanced for web apps!
* [Ace Editor] - awesome web-based text editor
* [markdown-it] - Markdown parser done right. Fast and easy to extend.
* [Twitter Bootstrap] - great UI boilerplate for modern web apps
* [node.js] - evented I/O for the backend
* [Express] - fast node.js network app framework [@tjholowaychuk]
* [Gulp] - the streaming build system
* [Breakdance](http://breakdance.io) - HTML to Markdown converter
* [jQuery] - duh

And of course Dillinger itself is open source with a [public repository][dill]
 on GitHub.

### Installation

Dillinger requires [Node.js](https://nodejs.org/) v4+ to run.

Install the dependencies and devDependencies and start the server.

```sh
$ cd dillinger
$ npm install -d
$ node app
```

For production environments...

```sh
$ npm install --production
$ npm run predeploy
$ NODE_ENV=production node app
```

### Plugins

Dillinger is currently extended with the following plugins. Instructions on how to use them in your own application are linked below.

| Plugin | README |
| ------ | ------ |
| Dropbox | [plugins/dropbox/README.md] [PlDb] |
| Github | [plugins/github/README.md] [PlGh] |
| Google Drive | [plugins/googledrive/README.md] [PlGd] |
| OneDrive | [plugins/onedrive/README.md] [PlOd] |
| Medium | [plugins/medium/README.md] [PlMe] |
| Google Analytics | [plugins/googleanalytics/README.md] [PlGa] |


### Development

Want to contribute? Great!

Dillinger uses Gulp + Webpack for fast developing.
Make a change in your file and instantanously see your updates!

Open your favorite Terminal and run these commands.

First Tab:
```sh
$ node app
```

Second Tab:
```sh
$ gulp watch
```

(optional) Third:
```sh
$ karma test
```
#### Building for source
For production release:
```sh
$ gulp build --prod
```
Generating pre-built zip archives for distribution:
```sh
$ gulp build dist --prod
```
### Docker
Dillinger is very easy to install and deploy in a Docker container.

By default, the Docker will expose port 80, so change this within the Dockerfile if necessary. When ready, simply use the Dockerfile to build the image.

```sh
cd dillinger
docker build -t joemccann/dillinger:${package.json.version}
```
This will create the dillinger image and pull in the necessary dependencies. Be sure to swap out `${package.json.version}` with the actual version of Dillinger.

Once done, run the Docker image and map the port to whatever you wish on your host. In this example, we simply map port 8000 of the host to port 80 of the Docker (or whatever port was exposed in the Dockerfile):

```sh
docker run -d -p 8000:8080 --restart="always" <youruser>/dillinger:${package.json.version}
```

Verify the deployment by navigating to your server address in your preferred browser.

```sh
127.0.0.1:8000
```

#### Kubernetes + Google Cloud

See [KUBERNETES.md](https://github.com/joemccann/dillinger/blob/master/KUBERNETES.md)


### Todos

 - Write MOAR Tests
 - Add Night Mode

License
----

MIT


**Free Software, Hell Yeah!**

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)


   [dill]: <https://github.com/joemccann/dillinger>
   [git-repo-url]: <https://github.com/joemccann/dillinger.git>
   [john gruber]: <http://daringfireball.net>
   [df1]: <http://daringfireball.net/projects/markdown/>
   [markdown-it]: <https://github.com/markdown-it/markdown-it>
   [Ace Editor]: <http://ace.ajax.org>
   [node.js]: <http://nodejs.org>
   [Twitter Bootstrap]: <http://twitter.github.com/bootstrap/>
   [jQuery]: <http://jquery.com>
   [@tjholowaychuk]: <http://twitter.com/tjholowaychuk>
   [express]: <http://expressjs.com>
   [AngularJS]: <http://angularjs.org>
   [Gulp]: <http://gulpjs.com>

   [PlDb]: <https://github.com/joemccann/dillinger/tree/master/plugins/dropbox/README.md>
   [PlGh]: <https://github.com/joemccann/dillinger/tree/master/plugins/github/README.md>
   [PlGd]: <https://github.com/joemccann/dillinger/tree/master/plugins/googledrive/README.md>
   [PlOd]: <https://github.com/joemccann/dillinger/tree/master/plugins/onedrive/README.md>
   [PlMe]: <https://github.com/joemccann/dillinger/tree/master/plugins/medium/README.md>
   [PlGa]: <https://github.com/RahulHP/dillinger/blob/master/plugins/googleanalytics/README.md>
