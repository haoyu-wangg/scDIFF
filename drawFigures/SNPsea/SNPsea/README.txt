# SNPsea Data Files

Kamil Slowikowski <slowikow@broadinstitute.org>

May 20, 2014

This archive contains data for:

    Slowikowski, et al. SNPsea: an algorithm to identify cell types, tissues,
    and pathways that affect risk loci. Bioinformatics (2014).
    doi:10.1093/bioinformatics/btu326

Please see the SNPsea documentation for more details:

    http://www.broadinstitute.org/mpg/snpsea

## Data

```bash
cd snpsea
curl -LOk http://files.figshare.com/1382662/SNPsea_data_20140212.zip
unzip SNPsea_data_20140212.zip
```

Download the compressed archive with data required to perform this analysis
(138M). The direct link to the zip shown above may be out of date and fail
to load. If so, please visit the link below instead:

<http://dx.doi.org/10.6084/m9.figshare.871430>

Contents of the compressed archive with data:

    Celiac_disease-Trynka2011-35_SNPs.gwas
    HDL_cholesterol-Teslovich2010-46_SNPs.gwas
    Multiple_sclerosis-IMSGC-51_SNPs.gwas
    Red_blood_cell_count-Harst2012-45_SNPs.gwas

    GeneAtlas2004.gct.gz  # Gene Atlas 2004 Affymetrix expression matrix
    ImmGen2012.gct.gz     # ImmGen 2012 Affymetrix expression matrix
    FANTOM2014.gct.gz     # FANTOM5 2014 CAGE matrix
    GO2013.gct.gz         # Gene Ontology 2013 binary annotation matrix

    NCBIgenes2013.bed.gz  # NCBI gene intervals
    Lango2010.txt.gz      # LD-pruned SNPs
    TGP2011.bed.gz        # 1000 Genomes Project SNP linkage intervals


### SNP sets

Phenotype             SNPs  Loci  Reference
---------             ----  ----  ---------
Celiac disease        35    34    Table 2 ([Trynka, *et al.* 2011][Trynka2011])
HDL cholesterol       46    46    Supp. Table 2 ([Teslovich, *et al.* 2010][Teslovich2010])
Multiple sclerosis    51    47    Supp. Table A ([IMSGC WTCCC 2011][IMSGC2011])
Red blood cell count  45    45    Table 1 ([Harst *et al.* 2012][Harst2012])


#### Celiac_disease-Trynka2011-35_SNPs.gwas

35 SNPs associated with Celiac disease taken from Table 2.
Positions are on hg19. All SNPs have $P \le 5e-8$.

> Trynka G, Hunt KA, Bockett NA, et al. Dense genotyping identifies and
> localizes multiple common and rare variant association signals in celiac
> disease. Nat Genet. 2011;43(12):1193-201. [Pubmed][Trynka2011]


#### HDL_cholesterol-Teslovich2010-46_SNPs.gwas

46 SNPs associated with HDL taken from Supplementary Table 2.
Positions are on hg19. All SNPs have $P \le 5e-8$.

> Teslovich TM, Musunuru K, Smith AV, et al. Biological, clinical and
> population relevance of 95 loci for blood lipids. Nature.
> 2010;466(7307):707-13. [Pubmed][Teslovich2010]


#### Multiple_sclerosis-IMSGC-51_SNPs.gwas

51 SNPs associated with Multiple Sclerosis taken from Supplementary Table A.
Positions are on hg19. All SNPs have $P \le 5e-8$.

> Sawcer S, Hellenthal G, Pirinen M, et al. Genetic risk and a primary role
> for cell-mediated immune mechanisms in multiple sclerosis. Nature.
> 2011;476(7359):214-9. [Pubmed][IMSGC2011]


#### Red_blood_cell_count-Harst2012-45_SNPs.gwas

45 SNPs associated with red blood cell count (RBC) taken from Table 1.
Positions are on hg19. All SNPs have $P \le 5e-8$.

> van der Harst P, Zhang W, Mateo leach I, et al. Seventy-five genetic loci
> influencing the human red blood cell. Nature. 2012;492(7429):369-75.
> [Pubmed][Harst2012]


### Gene matrices

Type    Genes  Conditions  Species                     Reference
----    -----  ----------  -------                     ---------
Affy    17581  79 tissues  homo sapiens                [GeneAtlas 2004][GeneAtlas]
Affy    15139  249 cells   mus musculus                [ImmGen 2012][ImmGen]
CAGE    18502  533 cells   homo sapiens                [FANTOM5 2014][FANTOM5]
Binary  19111  1751 terms  homo sapiens, mus musculus  [Gene Ontology] 2013, [Homologene]


#### GeneAtlas2004.gct.gz

Gene expression data for 79 human tissues from [GSE1133][GeneAtlas]. We
averaged the expression values for tissue replicates. For each gene, we
selected the single probe with the largest minimum value. Finally, we
converted the file to [GCT] format.

> Su AI et al. A gene atlas of the mouse and human protein-encoding
> transcriptomes. Proc Natl Acad Sci U S A, 2004 Apr 9;101(16):6062-7.
> [PubMed][GeneAtlas]


#### GO2013.gct.gz

A [GCT] formatted gene matrix with 1,751 annotation terms (1s and 0s
indicating presence or absence of the gene in a Gene Ontology term).

We downloaded the OBO file from [Gene Ontology] (data-version: 2013-06-29,
CVS revision: 9700).

For each gene, we climbed the hierarchy of ontology terms and applied
parental terms. If a gene is annotated with some term $T$, we also add
all of the terms that are parents of $T$. We copy terms between
homologous genes using [Homologene] data. If a mouse gene is annotated
with some term and the human homolog is not, then we copy the term to
the human gene. We discard all GO terms assigned to fewer than 100 or to
more than 1000 genes. This leaves us with a matrix of 19,111 genes and
1,751 terms.


#### ImmGen2012.gct.gz

Gene expression data for 249 blood cell types from [GSE15907][ImmGen]. We
averaged cell type replicates. For each gene, we selected the single probe
with the largest minimum.


#### FANTOM2014.gct.gz

CAGE data for 533 human cell types from [FANTOM5]. We averaged cell type
replicates. We discarded CAGE entries with 0 or multiple corresponding NCBI
Entrez IDs. Then, we summed the CAGE entries for each gene.


### LD-pruned SNPs and Genomic Intervals

#### Lango2010.txt.gz

A list of SNPs that span the whole genome, pruned by linkage disequilibrium
(LD). SNPsea samples null SNP sets matched on the number of genes in the
user's SNP set from this list. See this paper for more information:

> Lango allen H, Estrada K, Lettre G, et al. Hundreds of variants clustered
> in genomic loci and biological pathways affect human height. Nature.
> 2010;467(7317):832-8. [Pubmed][Lango2010]


#### NCBIgenes2013.bed.gz

All human start and stop positions taken from:

> <ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz>


#### TGP2011.bed.gz

Linkage intervals for a filtered set of SNPs from the [1000 Genomes
Project][tgp] Phase 1 (May 21, 2011). We downloaded a filtered (diallelic and
5 or more copies of the minor allele) set of markers from the [BEAGLE] website
and calculated pairwise LD (EUR) for all SNPs in a 1 Mb sliding window. The
linkage intervals were extended to the nearest [HapMap] recombination hotspot
with >3 cM/Mb recombination rate.

[Harst2012]: http://www.ncbi.nlm.nih.gov/pubmed/23222517
[IMSGC2011]: http://www.ncbi.nlm.nih.gov/pubmed/21833088
[Lango2010]: http://www.ncbi.nlm.nih.gov/pubmed/20881960
[Trynka2011]: http://www.ncbi.nlm.nih.gov/pubmed/22057235
[Teslovich2010]: http://www.ncbi.nlm.nih.gov/pubmed/20686565

[GeneAtlas]: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1133
[ImmGen]: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15907
[FANTOM5]: http://fantom.gsc.riken.jp/5/data/

[Gene Ontology]: http://www.geneontology.org
[tgp]: http://www.1000genomes.org/
[HapMap]: http://hapmap.ncbi.nlm.nih.gov/downloads/
[NCBI]: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/
[Homologene]: http://www.ncbi.nlm.nih.gov/homologene

[BEAGLE]: http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes.phase1_release_v3
[Bash]: http://www.gnu.org/software/bash/manual/bashref.html
[GCT]: http://www.broadinstitute.org/cancer/software/genepattern/gp_guides/file-formats/sections/gct


