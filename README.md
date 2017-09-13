## Usage


### Running the script maunally
Paired-end
```
bash salmonella_SNP_PCR.sh \
    -q -m \
    -p SE_clades_primers_21-mer.fasta \
    -o output/ \
    sample1_L001_R1_001.fastq.gz \
    sample1_L001_R2_001.fastq.gz
```
Single-end
```
bash salmonella_SNP_PCR.sh \
    -q \
    -p SE_clades_primers_21-mer.fasta \
    -o output/ \
    sample1_L001_R1_001.fastq.gz .fastq.gz
```
Short version (16 SNPs instead of 60; must enable the short list in the code)
```
bash salmonella_SNP_PCR.sh \
    -q -m \
    -p primers_21-mers_short-version.txt \
    -o output/ \
    sample1_L001_R1_001.fastq.gz \
    sample1_L001_R2_001.fastq.gz
```
