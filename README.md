# [FDRBench](https://doi.org/10.1101/2024.06.01.596967)
**FDRBench** is a tool for false discovery rate (FDR) control evaluation in proteomics. It provides two main functions: (1) build entrapment databases using randomly shuffled target sequences or using sequences from foreign species, and (2) estimate the false discovery proportion (FDP) using the lower bound, combined, and paired methods.

![Downloads](https://img.shields.io/github/downloads/Noble-Lab/FDRBench/total.svg) ![Release](https://img.shields.io/github/release/Noble-Lab/FDRBench.svg)![Downloads](https://img.shields.io/github/downloads/Noble-Lab/FDRBench/latest/total)

## Installation

FDRBench is written using Java and can be run on Windows, Mac OS and Linux. Only Java is required to be installed to run FDRBench. If java is not installed, please install Java by following the instruction at https://openjdk.org/install/ or https://www.oracle.com/java/technologies/downloads/. After java is installed, FDRBench can be downloaded at https://github.com/Noble-Lab/FDRBench/releases.

## Usage

```
$ java -jar fdrbench-0.0.2.jar
usage: Options
 -db <arg>                 Protein database file
 -crux <arg>               The peptide list generated by Crux
 -target_pep <arg>         The target peptide list to consider
 -decoy                    Add decoy or not
 -clip_n_m                 When digesting a protein starting with amino acid M, two copies of the
                           leading peptides (with and without the N-terminal M) are considered or
                           not. Default is false.
 -fix_nc <arg>             Fix N/C terminal amino acid. n/N: only N terminal, c/C: only C terminal
                           and nc/NC/cn/CN: both (default)
 -o <arg>                  Output file
 -I2L                      Convert I to L
 -diann                    For DIA-NN
 -uniprot                  For Uniprot
 -prosit                   Generate Prosit input file
 -charge <arg>             For prosit input: charge range, 2,3,4
 -nce <arg>                For prosit input: NCE
 -enzyme <arg>             Enzyme used for protein digestion. 0:Non enzyme, 1:Trypsin, 2:Trypsin (no
                           P rule) (default), 3:Arg-C, 4:Arg-C (no P rule), 5:Arg-N, 6:Glu-C,
                           7:Lys-C
 -miss_c <arg>             The max missed cleavages, default is 1
 -minLength <arg>          The minimum length of peptide to consider, default is 7
 -maxLength <arg>          The maximum length of peptide to consider, default is 35
 -export_db                Export protein database or not
 -seed <arg>               Random seed for generating decoy peptides
 -fix_seed                 Use a fixed random seed for all decoy peptides generation
 -fold <arg>               The number of folds for generating entrapment proteins/peptides
 -r <arg>                  For FDP calculation: #entrapment/#target
 -use_v1 <arg>             Use the first version of FDP calculation for 1-fold
 -pick <arg>               If a group has multiple proteins, how to peak one protein: first
                           (default),last,random
 -a                        Generate entrapment protein(s) for each target protein independently
 -check                    Checking for duplicates and random shuffling
 -method <arg>             The method to generate a random peptide: 0:shuffle (default), 1:swap
 -ms <arg>                 Multiple species entrapment: Fasta files of foreign species
 -ns                       no shared peptides between entrapment and target protein
 -swap                     Reverse the order of generated random peptide sequences
 -i <arg>                  PSM/peptide/precursor/protein file
 -score <arg>              The score name for ranking precursor/peptide/protein for FDP calculation.
                           The format could be "score", "score:0" or "score:1". The second part is 0
                           or 1, 0: lower is better, 1: higher is better
 -level <arg>              PSM, peptide, precursor or protein
 -pep <arg>                peptide/protein pair file
 -decoy_label <arg>        Label for decoy: rev_ in default
 -decoy_pos <arg>          Position of decoy label: 0 (start, in default); 1 (end)
 -entrapment_label <arg>   Label for entrapment: _p_target in default
 -entrapment_pos <arg>     Position of entrapment label: 0 (start); 1 (end, in default)
 -debug                    Print detailed information for debugging
 -h                        Print this usage information
```

#### Build entrapment databases

##### Build entrapment databases using randomly shuffled target sequences - peptide level:

Generate a **peptide level** entrapment database using the human target database `UP000005640_9606.fasta` (~20k human proteins). The database can be downloaded from [UniProt](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz).

```shell
java -jar fdrbench-0.0.1.jar -I2L -level peptide -db UP000005640_9606.fasta -o UP000005640_9606_entrapment_pep.txt -uniprot -diann -fix_nc c
```
Using the above command line, a peptide level entrapment database in FASTA format "UP000005640_9606_entrapment_pep.fasta" and a peptide tsv format file "UP000005640_9606_entrapment_pep.txt" will be generated:

The format of the Fasta file "UP000005640_9606_entrapment_pep.fasta" looks like below. The Fasta header of each peptide sequence in the file contains "_target" when the sequence is an original target or "_p_target" when the sequence is an entrapment peptide.

<code>
<pre>
>sp|PSLDQLAAHPWMLGADGGVPESCDLR<mark>_target</mark>|PSLDQLAAHPWMLGADGGVPESCDLR<mark>_target</mark>
PSLDQLAAHPWMLGADGGVPESCDLR
>sp|PALLAVGGADSLLEDGHQPCSWDMPR<mark>_p_target</mark>|PALLAVGGADSLLEDGHQPCSWDMPR<mark>_p_target</mark>
PALLAVGGADSLLEDGHQPCSWDMPR
>sp|QLQGASWELQSLR<mark>_target</mark>|QLQGASWELQSLR<mark>_target</mark>
QLQGASWELQSLR
>sp|QQSWLSLQGLEAR<mark>_p_target</mark>|QQSWLSLQGLEAR<mark>_p_target</mark>
QQSWLSLQGLEAR
>sp|YPERDNR<mark>_target</mark>|YPERDNR<mark>_target</mark>
YPERDNR
</code>
</pre>


The format of the tsv file "UP000005640_9606_entrapment_pep.txt" looks like below. This file is used in FDP calculation.

<code>
<pre>
sequence                    decoy  proteins                                 peptide_type  peptide_pair_index
PSLDQLAAHPWMLGADGGVPESCDLR  No     sp|Q86V86|PIM3_HUMAN                     target        0
PALLAVGGADSLLEDGHQPCSWDMPR  No     sp|Q86V86<mark>_p_target</mark>|PIM3_HUMAN<mark>_p_target</mark>   p_target      0
QLQGASWELQSLR               No     sp|Q96N95|ZN396_HUMAN                    target        1
QQSWLSLQGLEAR               No     sp|Q96N95<mark>_p_target</mark>|ZN396_HUMAN<mark>_p_target</mark>  p_target      1
YPERDNR                     No     sp|Q5T200|ZC3HD_HUMAN                    target        2
YDRPENR                     No     sp|Q5T200<mark>_p_target</mark>|ZC3HD_HUMAN<mark>_p_target</mark>  p_target      2
SYKALADQMNLLLSK             No     sp|Q9UBJ2|ABCD2_HUMAN                    target        3
YSALSNMDLQKLALK             No     sp|Q9UBJ2<mark>_p_target</mark>|ABCD2_HUMAN<mark>_p_target</mark>  p_target      3
</code>
</pre>

Below please find the description of each column in the output tsv file "UP000005640_9606_entrapment_pep.txt":

| Column name  | Description |
| ------------ | ----------- |
| sequence | peptide sequence (original target peptide or entrapment peptide) |
| decoy | decoy peptide (**Yes**) or not (**No**)  |
| proteins | protein ID (multiple IDs are separated by ";") |
| peptide_type | peptide type, original target (**target**) or entrapment (**p_target**) |
| peptide_pair_index | peptide pair index: an original target peptide and its paired entrapment have the same index |


The above example command line took about 30 seconds to run on a Mac MacBook computer.

##### Build entrapment databases using randomly shuffled target sequences - protein level:

Generate a **protein level** entrapment database using the human target database `UP000005640_9606.fasta` (~20k human proteins). The database can be downloaded from [UniProt](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz):
```shell
java -jar fdrbench-0.0.1.jar -level protein -db UP000005640_9606.fasta -o UP000005640_9606_I2L_entrapment.fasta -I2L -diann -uniprot -fix_nc c -check
```
Using the above command line, a protein level entrapment database in FASTA format "UP000005640_9606_I2L_entrapment.fasta" will be generated. For each target human protein, an entrapment protein is generated.

The format of the Fasta file "UP000005640_9606_I2L_entrapment.fasta" looks like below. The first protein as shown below is a target human protein while the second protein is its paired entrapment protein. Both command line parameters `-diann` and `-uniprot` are set in the command line, so an entrapment label "_p_target" is added in three positions of the Fasta protein header. This format is compatible with [DIA-NN](https://github.com/vdemichev/DiaNN) search.

<code>
<pre>
>sp|A0A0G2JMI3|HV692_HUMAN Immunoglobulin heavy variable 1-69-2 OS=Homo sapiens OX=9606 GN=IGHV1-69-2 PE=3 SV=2
MDCTWRLLLLVAAATGTHAEVQLVQSGAEVKKPGATVKLSCKVSGYTFTDYYMHWVQQAPGKGLEWMGLVDPEDGETLYAEKFQGRVTLTADTSTDTAYMELSSLRSEDTAVYYCAT
>sp|A0A0G2JMI3<mark>_p_target</mark>|HV692_HUMAN<mark>_p_target</mark> Immunoglobulin heavy variable 1-69-2 OS=Homo sapiens OX=9606 GN=IGHV1-69-2<mark>_p_target</mark> PE=3 SV=2
CWTMDRTGAVEAVQVQAESAGLALTVLHLLKKTAGPVKCSLKYYMATPVDWTQQFSYVGHGKDAPGEEELDVMEGLGYLWTKQFGREAMLTLADYVLDSTTSTSTRDVTACESAYYT
</code>
</pre>

The above example command line took about 10 seconds to run on a Mac MacBook computer.

#### FDP calculation

##### Peptide (or precursor) level FDR control evaluation

For **precursor level** FDP calculation, an example input is shown below:

```
run                                                peptide                             mod_peptide                                     charge  q_value  PEP  protein                                    score
20230406_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_I_1  AAAPAPEEEMDECEQALAAEPK              AAAPAPEEEMDEC(UniMod:4)EQALAAEPK                3       1e-9     0    AAAPAPEEEMDECEQALAAEPK_target              1
20230406_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_I_1  AAASGQPRPEMQCPAEHEEDMYR             AAASGQPRPEMQC(UniMod:4)PAEHEEDMYR               4       1e-9     0    AAASGQPRPEMQCPAEHEEDMYR_target             2
20230406_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_I_1  AAEVWMDEYKNFYYAAVPSAR               AAEVWMDEYKNFYYAAVPSAR                           3       1e-9     0    AAEVWMDEYKNFYYAAVPSAR_target               3
20230406_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_I_1  AAFTVSLDPGPLEQFPHSMEPQLR            AAFTVSLDPGPLEQFPHSMEPQLR                        3       1e-9     0    AAFTVSLDPGPLEQFPHSMEPQLR_target            4
20230406_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_I_1  AAHVLMPHESTVEHTHVDLNEMESPLATR       AAHVLMPHESTVEHTHVDLNEMESPLATR                   4       1e-9     0    AAHVLMPHESTVEHTHVDLNEMESPLATR_target       5
20230406_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_I_1  AALLSPGDPALWAGLMAACHADDKLALVNNTQPK  AALLSPGDPALWAGLMAAC(UniMod:4)HADDKLALVNNTQPK    4       1e-9     0    AALLSPGDPALWAGLMAACHADDKLALVNNTQPK_target  6
20230406_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_I_1  AALQQKENLPVSSDGNLPQQAASAPSR         AALQQKENLPVSSDGNLPQQAASAPSR                     3       1e-9     0    AALQQKENLPVSSDGNLPQQAASAPSR_target         7
20230406_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_I_1  ADMLCNSQSNDLLQHQGSNCGGTSNK          ADMLC(UniMod:4)NSQSNDLLQHQGSNC(UniMod:4)GGTSNK  3       1e-9     0    ADMLCNSQSNDLLQHQGSNCGGTSNK_target          8
20230406_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_I_1  AFLADPSAFVAAAPVAAATTAAPAAAAAPAK     AFLADPSAFVAAAPVAAATTAAPAAAAAPAK                 3       1e-9     0    AFLADPSAFVAAAPVAAATTAAPAAAAAPAK_target     9
```

The required columns are described below. Decoy hits shouldn't be included in the input file for FDP calcualtion.
| Column name  | Description |
| ------------ | ----------- |
| peptide | peptide sequence |
| mod_peptide | peptide modification (optional for peptide level FDR control evaluation) |
| charge | precursor charge (optional for peptide level FDR control evaluation) |
| q_value | FDR reported by the search engine used to generate the result |
| score | precursor or peptide score used to rank precursor or peptide for FDR calculation |


Below is an example command line to run precursor FDP calculation:

```shell
java -jar fdrbench-0.0.1.jar  -i peptide-fdp_precursor_input.tsv -fold 1 -pep UP000005640_9606_entrapment_pep.txt -level precursor -o peptide-diann_fdp_precursor.csv -score 'score:0' 

```

The output file "peptide-diann_fdp_precursor.csv" contains estimated FDP values using different methods. The file looks like below:

```
run                                                peptide                             mod_peptide                                     charge  q_value  PEP  protein                                    score  combined_fdp  n_t  n_p  paired_fdp  n_p_t_s  n_p_s_t  vt  lower_bound_fdp
20230406_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_I_1  AAAPAPEEEMDECEQALAAEPK              AAAPAPEEEMDEC(UniMod:4)EQALAAEPK                3       1.0E-9   0.0  AAAPAPEEEMDECEQALAAEPK_target              1      0.0           393  0    0.0         0        0        0   0.0
20230406_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_I_1  AAASGQPRPEMQCPAEHEEDMYR             AAASGQPRPEMQC(UniMod:4)PAEHEEDMYR               4       1.0E-9   0.0  AAASGQPRPEMQCPAEHEEDMYR_target             2      0.0           393  0    0.0         0        0        0   0.0
20230406_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_I_1  AAEVWMDEYKNFYYAAVPSAR               AAEVWMDEYKNFYYAAVPSAR                           3       1.0E-9   0.0  AAEVWMDEYKNFYYAAVPSAR_target               3      0.0           393  0    0.0         0        0        0   0.0
20230406_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_I_1  AAFTVSLDPGPLEQFPHSMEPQLR            AAFTVSLDPGPLEQFPHSMEPQLR                        3       1.0E-9   0.0  AAFTVSLDPGPLEQFPHSMEPQLR_target            4      0.0           393  0    0.0         0        0        0   0.0
20230406_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_I_1  AAHVLMPHESTVEHTHVDLNEMESPLATR       AAHVLMPHESTVEHTHVDLNEMESPLATR                   4       1.0E-9   0.0  AAHVLMPHESTVEHTHVDLNEMESPLATR_target       5      0.0           393  0    0.0         0        0        0   0.0
20230406_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_I_1  AALLSPGDPALWAGLMAACHADDKLALVNNTQPK  AALLSPGDPALWAGLMAAC(UniMod:4)HADDKLALVNNTQPK    4       1.0E-9   0.0  AALLSPGDPALWAGLMAACHADDKLALVNNTQPK_target  6      0.0           393  0    0.0         0        0        0   0.0
20230406_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_I_1  AALQQKENLPVSSDGNLPQQAASAPSR         AALQQKENLPVSSDGNLPQQAASAPSR                     3       1.0E-9   0.0  AALQQKENLPVSSDGNLPQQAASAPSR_target         7      0.0           393  0    0.0         0        0        0   0.0
20230406_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_I_1  ADMLCNSQSNDLLQHQGSNCGGTSNK          ADMLC(UniMod:4)NSQSNDLLQHQGSNC(UniMod:4)GGTSNK  3       1.0E-9   0.0  ADMLCNSQSNDLLQHQGSNCGGTSNK_target          8      0.0           393  0    0.0         0        0        0   0.0
20230406_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_I_1  AFLADPSAFVAAAPVAAATTAAPAAAAAPAK     AFLADPSAFVAAAPVAAATTAAPAAAAAPAK                 3       1.0E-9   0.0  AFLADPSAFVAAAPVAAATTAAPAAAAAPAK_target     9      0.0           393  0    0.0         0        0        0   0.0
```

Below please find the description of the main columns in the output tsv file "peptide-diann_fdp_precursor.csv":

| Column name  | Description |
| ------------ | ----------- |
| peptide | peptide sequence |
| mod_peptide | peptide modification |
| charge | precursor charge |
| q_value | FDR reported by the search engine used to generate the result |
| score | precursor score used to rank precursor for FDR calculation |
| combined_fdp | the FDP estimated using the combined method |
| paired_fdp | the FDP estimated using the paired method |
| lower_bound_fdp | the lower bound FDP |

The above example command line took about 30 seconds to run on a Mac MacBook computer.

##### Protein level FDR control evaluation

For **protein level** FDP calculation, an example input is shown below. Only the three columns **protein**, **q_value** and **score** are required. Decoy hits shouldn't be included in the input file for FDP calcualtion.

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

| Column name  | Description |
| ------------ | ----------- |
| protein | protein ID (multiple IDs are separated by ";"). Each row is a protein group |
| q_value | FDR reported by the search engine used to generate the result  |
| score | protein score used to rank protein for FDR calculation |



Below is an example command line to run protein FDP calculation. The database used to generate the identification result must be generated using FDRBench. The command line parameter **-score** is set to use the column "score" as a secondary ranking score in FDP calculation. The number **0** indicates lower score is better (more confident). The example file **protein-fdp_protein_input.tsv** is available in the [example](https://github.com/Noble-Lab/FDRBench/tree/master/example) folder.

```shell
java -jar fdrbench-0.0.1.jar -i example/protein-fdp_protein_input.tsv -level protein -o protein-diann_fdp_protein.csv -score 'score:0' -fold 1 -pick first
```

The output file "protein-diann_fdp_protein.csv" contains estimated FDP values using different methods. The file looks like below:

```
Protein.Group  PG.Q.Value  q_value     protein  score  combined_fdp  n_t   n_p  paired_fdp  n_p_t_s  n_p_s_t  vt  lower_bound_fdp
P62857         1.64447E-4  1.64447E-4  P62857   1      0.0           6079  0    0.0         0        0        0   0.0
P29590         1.64447E-4  1.64447E-4  P29590   2      0.0           6079  0    0.0         0        0        0   0.0
O95989         1.64447E-4  1.64447E-4  O95989   3      0.0           6079  0    0.0         0        0        0   0.0
P61006         1.64447E-4  1.64447E-4  P61006   4      0.0           6079  0    0.0         0        0        0   0.0
Q70J99         1.64447E-4  1.64447E-4  Q70J99   5      0.0           6079  0    0.0         0        0        0   0.0
```

Below please find the description of the main columns in the output tsv file "protein-diann_fdp_protein.csv":

| Column name  | Description |
| ------------ | ----------- |
| protein | protein ID (multiple IDs are separated by ";"). Each row is a protein group |
| q_value | FDR reported by the search engine used to generate the result  |
| score | protein score used to rank protein for FDR calculation |
| combined_fdp | the FDP estimated using the combined method |
| paired_fdp | the FDP estimated using the paired method |
| lower_bound_fdp | the lower bound FDP |

The above example command line took about 5 seconds to run on a Mac MacBook computer.



## How to cite:

Wen, Bo, et al. "[Assessment of false discovery rate control in tandem mass spectrometry analysis using entrapment.](https://doi.org/10.1101/2024.06.01.596967)" bioRxiv (2024): 2024-06.



