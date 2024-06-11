# [FDRBench](https://doi.org/10.1101/2024.06.01.596967)
**FDRBench** is a tool for FDR control evaluation in proteomics. It provides two main functions: (1) build entrapment databases using randomly shuffled target sequences or using sequences from foreign species, and (2) estimate the FDP using the lower bound, combined, and paired methods.

## Installation

FDRBench is written using Java. Only Java is required to be installed to run FDRBench.

## Usage

```
$ java -jar fdrbench-0.0.1.jar
usage: Options
 -db <arg>           Protein database file
 -crux <arg>         The peptide list generated by Crux
 -target_pep <arg>   The target peptide list to consider
 -decoy              Add decoy or not
 -clip_n_m           When digesting a protein starting with amino acid M, two copies of the leading
                     peptides (with and without the N-terminal M) are considered or not. Default is
                     false.
 -fix_nc <arg>       Fix N/C terminal amino acid. n/N: only N terminal, c/C: only C terminal and
                     nc/NC/cn/CN: both (default)
 -o <arg>            Output file
 -I2L                Convert I to L
 -diann              For DIA-NN
 -uniprot            For Uniprot
 -prosit             Generate Prosit input file
 -charge <arg>       For prosit input: charge range, 2,3,4
 -nce <arg>          For prosit input: NCE
 -enzyme <arg>       Enzyme used for protein digestion. 0:Non enzyme, 1:Trypsin, 2:Trypsin (no P
                     rule) (default), 3:Arg-C, 4:Arg-C (no P rule), 5:Arg-N, 6:Glu-C, 7:Lys-C
 -miss_c <arg>       The max missed cleavages, default is 1
 -minLength <arg>    The minimum length of peptide to consider, default is 7
 -maxLength <arg>    The maximum length of peptide to consider, default is 35
 -export_db          Export protein database or not
 -seed <arg>         Random seed for generating decoy peptides
 -fix_seed           Use a fixed random seed for all decoy peptides generation
 -fold <arg>         The number of folds for generating entrapment proteins/peptides
 -r <arg>            For FDP calculation: #entrapment/#target
 -use_v1 <arg>       Use the first version of FDP calculation for 1-fold
 -pick <arg>         If a group has multiple proteins, how to peak one protein: first
                     (default),last,random
 -a                  Generate entrapment protein(s) for each target protein independently
 -check              Checking for duplicates and random shuffling
 -method <arg>       The method to generate a random peptide: 0:shuffle (default), 1:swap
 -ms <arg>           Multiple species entrapment: Fasta files of foreign species
 -ns                 no shared peptides between entrapment and target protein
 -swap               Reverse the order of generated random peptide sequences
 -i <arg>            PSM/peptide/precursor/protein file
 -score <arg>        The score name for ranking precursor/peptide/protein for FDP calculation
 -level <arg>        PSM, peptide, precursor or protein
 -pep <arg>          peptide/protein pair file
 -debug              Print detailed information for debugging
 -h                  Print this usage information
```

#### Build entrapment databases

##### Build entrapment databases using randomly shuffled target sequences - peptide level:

Generate a peptide level entrapment database using the target database `UP000005640_9606.fasta`:
```shell
java -jar fdrbench-0.0.1.jar -I2L -level peptide -db UP000005640_9606.fasta -o UP000005640_9606_entrapment_pep.txt -uniprot -diann -fix_nc c
```
Using the above command line, a peptide level entrapment database in FASTA format "UP000005640_9606_entrapment_pep.fasta" and a peptide tsv format file "UP000005640_9606_entrapment_pep.txt" will be generated:

The format of the Fasta file "UP000005640_9606_entrapment_pep.fasta" looks like below:
```
>sp|PSLDQLAAHPWMLGADGGVPESCDLR_target|PSLDQLAAHPWMLGADGGVPESCDLR_target
PSLDQLAAHPWMLGADGGVPESCDLR
>sp|PALLAVGGADSLLEDGHQPCSWDMPR_p_target|PALLAVGGADSLLEDGHQPCSWDMPR_p_target
PALLAVGGADSLLEDGHQPCSWDMPR
>sp|QLQGASWELQSLR_target|QLQGASWELQSLR_target
QLQGASWELQSLR
>sp|QQSWLSLQGLEAR_p_target|QQSWLSLQGLEAR_p_target
QQSWLSLQGLEAR
>sp|YPERDNR_target|YPERDNR_target
YPERDNR
```

The format of the tsv file "UP000005640_9606_entrapment_pep.txt" looks like below. This file is used in FDP calculation.

```
sequence                    decoy  proteins                                 peptide_type  peptide_pair_index
PSLDQLAAHPWMLGADGGVPESCDLR  No     sp|Q86V86|PIM3_HUMAN                     target        0
PALLAVGGADSLLEDGHQPCSWDMPR  No     sp|Q86V86_p_target|PIM3_HUMAN_p_target   p_target      0
QLQGASWELQSLR               No     sp|Q96N95|ZN396_HUMAN                    target        1
QQSWLSLQGLEAR               No     sp|Q96N95_p_target|ZN396_HUMAN_p_target  p_target      1
YPERDNR                     No     sp|Q5T200|ZC3HD_HUMAN                    target        2
YDRPENR                     No     sp|Q5T200_p_target|ZC3HD_HUMAN_p_target  p_target      2
SYKALADQMNLLLSK             No     sp|Q9UBJ2|ABCD2_HUMAN                    target        3
YSALSNMDLQKLALK             No     sp|Q9UBJ2_p_target|ABCD2_HUMAN_p_target  p_target      3
```

##### Build entrapment databases using randomly shuffled target sequences - protein level:

Generate a peptide level entrapment database using the target database `UP000005640_9606.fasta`:
```shell
java -jar fdrbench-0.0.1.jar -level protein -db UP000005640_9606.fasta -o UP000005640_9606_I2L_entrapment.fasta -I2L -diann -uniprot -fix_nc c -check
```
Using the above command line, a protein level entrapment database in FASTA format "UP000005640_9606_I2L_entrapment.fasta" will be generated:

The format of the Fasta file "UP000005640_9606_I2L_entrapment.fasta" looks like below:
```
>sp|A0A0G2JMI3|HV692_HUMAN Immunoglobulin heavy variable 1-69-2 OS=Homo sapiens OX=9606 GN=IGHV1-69-2 PE=3 SV=2
MDCTWRLLLLVAAATGTHAEVQLVQSGAEVKKPGATVKLSCKVSGYTFTDYYMHWVQQAPGKGLEWMGLVDPEDGETLYAEKFQGRVTLTADTSTDTAYMELSSLRSEDTAVYYCAT
>sp|A0A0G2JMI3_p_target|HV692_HUMAN_p_target Immunoglobulin heavy variable 1-69-2 OS=Homo sapiens OX=9606 GN=IGHV1-69-2_p_target PE=3 SV=2
CWTMDRTGAVEAVQVQAESAGLALTVLHLLKKTAGPVKCSLKYYMATPVDWTQQFSYVGHGKDAPGEEELDVMEGLGYLWTKQFGREAMLTLADYVLDSTTSTSTRDVTACESAYYT
```

#### FDP calculation

For protein level FDP calculation, an example input is shown below. Only the three columns **protein**, **q_value** and **score** are required.

```
Protein.Group  PG.Q.Value  q_value     protein  score
P62857         1.64447e-4  1.64447e-4  P62857   1
P29590         1.64447e-4  1.64447e-4  P29590   2
O95989         1.64447e-4  1.64447e-4  O95989   3
P61006         1.64447e-4  1.64447e-4  P61006   4
Q70J99         1.64447e-4  1.64447e-4  Q70J99   5
P36404         1.64447e-4  1.64447e-4  P36404   6
Q71SY5         1.64447e-4  1.64447e-4  Q71SY5   7
Q9BPW8         1.64447e-4  1.64447e-4  Q9BPW8   8
Q8WUP2         1.64447e-4  1.64447e-4  Q8WUP2   9
```

Below is an example command line to run protein FDP calculation. The database used to generate the identification result must be generated using FDRBench.

```shell
java -jar fdrbench-0.0.1.jar -i protein-fdp_protein_input.tsv -level protein -o protein-diann_fdp_protein.csv -score 'score:0' -fold 1 -pick first
```

The output file "protein-diann_fdp_protein.csv" contains estimated FDP values using different methods.

## How to cite:

Wen, Bo, et al. "[Assessment of false discovery rate control in tandem mass spectrometry analysis using entrapment.](https://doi.org/10.1101/2024.06.01.596967)" bioRxiv (2024): 2024-06.



