options=(
    --snps              My.gwas
    --gene-matrix       GeneAtlas2004.gct.gz
    --gene-intervals    NCBIgenes2013.bed.gz
    --snp-intervals     TGP2011.bed.gz
    --null-snps         Lango2010.txt.gz
    --out               out
    --slop              10e3
    --threads           4
    --null-snpsets      0
    --min-observations  100
    --max-iterations    1e7
)
snpsea ${options[*]}