# SV detection workflow

## Delly2 

Call SV in each sample

```sh
delly call -g $reference  -o $name.bcf $name.bam
```

Merging all results into a single file

```sh
delly merge -o all.sites.bcf $name1.bcf $name2.bcf $name3.bcf ...
```

Genotyping samples at each merged sv sites

```sh
delly call -g $reference -v all.sites.bcf -o $name1.geno.bcf $name1.bam
```

Merge genotyping results 

```sh
bcftools merge -m id -O b -o merge.final.bcf $name1.geno.bcf $name2.geno.bcf $name3.geno.bcf ...
```

Filter genotyping results using default QC field

```sh
for file in *geno.bcf.vcf;  do python3 fillter.py $file; done
```

filter.py:

```python
#!/usr/bin/python3
import sys
import re

infile = sys.argv[1]
name = re.search(r'(.*)\.vcf',infile)[1]
outfile = name + '.filter.vcf'
o = open(outfile, 'w')
with open(infile, 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('#'):
            # print(line)
            o.write(line + '\n')
            continue
        content = line.split()
        if content[6] != "PASS" or re.match(r'IMPRECISE', content[7]):
            continue
        # print(line)
        o.write(line + '\n')
```

## Lumpy/Smoove

Call SV

```sh
smoove call --outdir results-smoove/ --name $sample --fasta $reference -p 1 --genotype $name1.bam
```

Merging results into a single file

```sh
smoove merge --name merged -f $reference --outdir ./ results-smoove/*.genotyped.vcf.gz
```

Genotyping samples at each merged sv sites

```sh
smoove genotype -d -x -p 1 --name $name1-joint --outdir results-genotped/ --fasta $reference --vcf merged.sites.vcf.gz $name1.bam
```

## Manta

SV calling and genotyping of all samples

```sh
${MANTA_INSTALL_PATH}/bin/configManta.py --bam $name1.bam --bam $name2.bam --bam $name3.bam ... --referenceFasta $reference --runDir Manta_dir
```

## Merging SVs from different methods

```sh
svimmer --threads 10 all.vcf.list LG01 LG02 LG03 LG04 <......> # 染色体名称
all.vcf.list
```

## Re-genotyping across all samples

```sh
./graphtyper genotype_sv reference.fasta input.vcf.gz --sams=bam.list --region=scaffold:start-end
```