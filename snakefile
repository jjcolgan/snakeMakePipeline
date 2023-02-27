input = open("log", "r")
input = input.read()
SAMPLE = input.split("\n")
del SAMPLE[len(SAMPLE)-1]
rule all:
        input:
                geneFamiles = expand("/project/blekhman/jjcolgan/testFull/humannOutput/{sample}_genefamilies.tsv", sample = SAMPLE),
                pathAbudance = expand("/project/blekhman/jjcolgan/testFull/humannOutput/{sample}_pathabundance.tsv", sample = SAMPLE),
                pathCoverage = expand("/project/blekhman/jjcolgan/testFull/humannOutput/{sample}_pathcoverage.tsv", sample = SAMPLE)
rule trim:
        input:
                R1 ="/project/blekhman/sjarif/konzo/raw_data/Blekhman2_Project_002/{sample}_R1_001.fastq.gz",
                R2 ="/project/blekhman/sjarif/konzo/raw_data/Blekhman2_Project_002/{sample}_R2_001.fastq.gz"
        output:
                R1Trimmed = "/project/blekhman/jjcolgan/testFull/dehosted/{sample}_R1_001_trimmed.fastq.gz",
                R1Untrimmed = "/project/blekhman/jjcolgan/testFull/dehosted/{sample}_R1_001_untrimmed.fastq.gz",
                R2Trimmed = "/project/blekhman/jjcolgan/testFull/dehosted/{sample}_R2_001_trimmed.fastq.gz",
                R2Untrimmed = "/project/blekhman/jjcolgan/testFull/dehosted/{sample}_R2_001_untrimmed.fastq.gz"
        threads:12
        shell:
                """
                module load java
                java -jar /project/blekhman/jjcolgan/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads {threads} -phred33 \
                {input.R1} \
                {input.R2} \
                {output.R1Trimmed} \
                {output.R1Untrimmed} \
                {output.R2Trimmed} \
                {output.R2Untrimmed} \
                ILLUMINACLIP:/project/blekhman/jjcolgan/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:10 MINLEN:25
                """
rule bowtie_dehost:
        input:
                R1Trimmed = "/project/blekhman/jjcolgan/testFull/dehosted/{sample}_R1_001_trimmed.fastq.gz",
                R1Untrimmed = "/project/blekhman/jjcolgan/testFull/dehosted/{sample}_R1_001_untrimmed.fastq.gz",
                R2Trimmed = "/project/blekhman/jjcolgan/testFull/dehosted/{sample}_R2_001_trimmed.fastq.gz",
                R2Untrimmed = "/project/blekhman/jjcolgan/testFull/dehosted/{sample}_R2_001_untrimmed.fastq.gz"

        output:
                sam = "/project/blekhman/jjcolgan/testFull/dehosted/{sample}MappedAndUnmapped.sam"
        threads: 8
        shell:
                """
                bowtie2 -p {threads} -x /project/blekhman/jjcolgan/bowtieIndex/grch38_1kgmaj \
                -1 {input.R1Trimmed} \
                -2 {input.R2Trimmed} \
                -S {output.sam}
                rm {input.R1Untrimmed}
                rm {input.R2Untrimmed}
                """
rule convert_to_bam_dehost:
        input:
                sam = "/project/blekhman/jjcolgan/testFull/dehosted/{sample}MappedAndUnmapped.sam"
        output:
                bam = "/project/blekhman/jjcolgan/testFull/dehosted/{sample}MappedAndUnmapped.bam"
        shell:
                """
                samtools view -bS {input.sam} > {output.bam}
                rm {input.sam}
                """
rule filter_unmapped_dehost:
        input:
                both = "/project/blekhman/jjcolgan/testFull/dehosted/{sample}MappedAndUnmapped.bam"

        output:
                filtered = "/project/blekhman/jjcolgan/testFull/dehosted/{sample}Unmapped.bam"
        shell:
                """
                samtools view -b -f 12 -F 256 {input.both} -o {output.filtered}
                rm {input.both}
                """
rule sort_dehost:
        input:
                unsorted = "/project/blekhman/jjcolgan/testFull/dehosted/{sample}Unmapped.bam"
        output:
                sorted = "/project/blekhman/jjcolgan/testFull/dehosted/{sample}UnmappedSorted.bam"
        shell :
                """
                samtools sort -n {input.unsorted} -o {output.sorted}
                rm {input.unsorted}
                """
rule covert_dehost:
        input:
                sorted = "/project/blekhman/jjcolgan/testFull/dehosted/{sample}UnmappedSorted.bam"
        output:
                R1 = "/project/blekhman/jjcolgan/testFull/dehosted/{sample}_R1_001.fastq",
                R2 = "/project/blekhman/jjcolgan/testFull/dehosted/{sample}_R2_001.fastq"
        shell:
                """
                bedtools bamtofastq -i {input.sorted} -fq {output.R1} -fq2 {output.R2}
                rm {input.sorted}
                """
rule zip_dehost:
        input:
                R1 = "/project/blekhman/jjcolgan/testFull/dehosted/{sample}_R1_001.fastq",
                R2 = "/project/blekhman/jjcolgan/testFull/dehosted/{sample}_R2_001.fastq"
        output:
                R1 = "/project/blekhman/jjcolgan/testFull/dehosted/{sample}_R1_001.fastq.gz",
                R2 = "/project/blekhman/jjcolgan/testFull/dehosted/{sample}_R2_001.fastq.gz"
        shell:
                """
                gzip {input.R1}
                gzip {input.R2}
                """
rule bowtie_derrna:
        input:
                R1 = "/project/blekhman/jjcolgan/testFull/dehosted/{sample}_R1_001.fastq.gz",
                R2 = "/project/blekhman/jjcolgan/testFull/dehosted/{sample}_R2_001.fastq.gz"
        output:
                sam = "/project/blekhman/jjcolgan/testFull/derrna/{sample}.sam"
        threads: 8
        shell:
                "bowtie2 -p {threads} -x /project/blekhman/jjcolgan/rrnaDb/16sArcAndBacDb -1 {input.R1} -2 {input.R2} -S {output.sam}"
rule convert_to_bam_derrna:
        input:
                sam = "/project/blekhman/jjcolgan/testFull/derrna/{sample}.sam"
        output:
                bam = "/project/blekhman/jjcolgan/testFull/derrna/{sample}MappedAndUnmapped.bam"
        shell:
                """
                samtools view -bS {input.sam} > {output.bam}
                rm {input.sam}
                """
rule filter_unmapped_derrna:
        input:
                both = "/project/blekhman/jjcolgan/testFull/derrna/{sample}MappedAndUnmapped.bam"

        output:
                filtered = "/project/blekhman/jjcolgan/testFull/derrna/{sample}Filtered.bam"
        shell:
                """
                samtools view -b -f 12 -F 256 {input.both} -o {output.filtered}
                rm {input.both}
                """
rule sort_dehost_derrna:
        input:
                unsorted = "/project/blekhman/jjcolgan/testFull/derrna/{sample}Filtered.bam"
        output:
                sorted = "/project/blekhman/jjcolgan/testFull/derrna/{sample}UnmappedSorted.bam"
        shell :
                """
                samtools sort -n {input.unsorted} -o {output.sorted}
                rm {input.unsorted}
                """
rule covert_dehost_derrna:
        input:
                sorted = "/project/blekhman/jjcolgan/testFull/derrna/{sample}UnmappedSorted.bam"
        output:
                R1 = "/project/blekhman/jjcolgan/testFull/derrna/{sample}_R1_001.fastq",
                R2 = "/project/blekhman/jjcolgan/testFull/derrna/{sample}_R2_001.fastq"
        shell:
                """
                bedtools bamtofastq -i {input.sorted} -fq {output.R1} -fq2 {output.R2}
                rm {input.sorted}
                """
rule zip_dehost_derrna:
        input:
                R1 = "/project/blekhman/jjcolgan/testFull/derrna/{sample}_R1_001.fastq",
                R2 = "/project/blekhman/jjcolgan/testFull/derrna/{sample}_R2_001.fastq"
        output:
                R1 = "/project/blekhman/jjcolgan/testFull/derrna/{sample}_R1_001.fastq.gz",
                R2 = "/project/blekhman/jjcolgan/testFull/derrna/{sample}_R2_001.fastq.gz"
        shell:
                """
                gzip {input.R1}
                gzip {input.R2}
                """
rule metaphlan:
        input:
                R2= "/project/blekhman/jjcolgan/testFull/derrna/{sample}_R2_001.fastq.gz",
                R1 ="/project/blekhman/jjcolgan/testFull/derrna/{sample}_R1_001.fastq.gz"
        output:
                profile = "/project/blekhman/jjcolgan/testFull/profiles/{sample}Profile.txt",
                bowtie = "/project/blekhman/jjcolgan/testFull/profiles/{sample}.bowtie2.bz2}"
        shell:
                """
                metaphlan {input.R1},{input.R2} --bowtie2out {output.bowtie} \
                --index mpa_vJan21_CHOCOPhlAnSGB_202103 \
                --bowtie2db /project/blekhman/jjcolgan/humannDatabases \
                --nproc 8 \
                --input_type fastq -o {output.profile}
                """
rule mergeFastq:
        input:
                R1 ="/project/blekhman/jjcolgan/testFull/derrna/{sample}_R1_001.fastq.gz",
                R2 = "/project/blekhman/jjcolgan/testFull/derrna/{sample}_R2_001.fastq.gz"
        output:
                merged = "/project/blekhman/jjcolgan/testFull/merged/{sample}.fastq.gz"
        shell:
                """
                cat {input.R1} > {output.merged}
                cat {input.R2} >> {output.merged}
                """
rule humann3:
        input:
                merged = "/project/blekhman/jjcolgan/testFull/merged/{sample}.fastq.gz",
                profile ="/project/blekhman/jjcolgan/testFull/profiles/{sample}Profile.txt"

        output:
                geneFamiles = "/project/blekhman/jjcolgan/testFull/humannOutput/{sample}_genefamilies.tsv",
                pathAbudance = "/project/blekhman/jjcolgan/testFull/humannOutput/{sample}_pathabundance.tsv",
                pathCoverage = "/project/blekhman/jjcolgan/testFull/humannOutput/{sample}_pathcoverage.tsv"
        shell:
                """
                humann --nucleotide-database /project/blekhman/jjcolgan/humannDatabases/chocophlan --taxonomic-profile {input.profile} --protein-database \
                /project/blekhman/jjcolgan/humannDatabases/uniref50/uniref --threads 8 -i {input.merged} -o /project/blekhman/jjcolgan/testFull/humannOutput
                """
