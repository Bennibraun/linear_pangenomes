## Research Plan Outline

#### Part 1: Comparative analysis of pangenome methodologies

- Methodologies:
	- Our approach (linear pangenome consisting of SVs called from long read assemblies)
	- Conspecific linear reference genome
	- Heterospecific linear reference genome (selected from among closest available, ideally high quality and same family)
	- Minigraph-cactus graph pangenome built automatically from linear reference and long read assemblies
	- Polished graph pangenome based on high-quality assemblies (if available)
- Evaluation metrics:
	- Short read alignment rate and mapping quality
	- Allelic balance at heterozygous SNPs (quantifies reference bias, see Buyan et al. 2025)
	- Nucleotide diversity ($\pi$) (simple measure of polymorphism in a population)
	- Allele frequency spectra
	- $F_{st}$ scores (genetic differentiation) for short reads
	- $F_{ROH}$ scores (runs of homozygosity i.e. inbreeding coefficient) for short reads
	- Genome size (bp) indicating amount of novel sequence added
	- Number of SNPs called from short reads
		- Optionally, SNP-based PCA demonstrating effective grouping by geography, subspecies, etc.
	- Number of short reads aligned to non-reference sequence
- Timeline:
	- Paper submission in Summer 2026

#### Part 2: Variant effect prediction on genes in the linear pangenome
- A pipeline to annotate SNPs with VEP or SnpEff in genes
	- Comparative analysis similar to the above, demonstrating genotyping performance relative to linear reference and graph
	- Validation with RNA-seq when available
	- VEP may not work due to input requirements, potentially easier to make a custom database for SnpEff
- A pipeline for de novo re-annotation of genes overlapping SVs with a gene predictor e.g. AUGUSTUS when RNA-seq is available for validation
- Statistical outputs revealing common and rare variants affecting gene function in the population
	- To be used for hypothesis generation and downstream phenotype associations
	- e.g. $d_N/d_{S}$ ratio (balance between deleterious and beneficial mutations) for genes under selection
	- e.g. exploratory gene ontology / KEGG enrichments on SV genes
- Timeline:
	- Paper submission in 2027

#### General Obstacles
- Data quality
	- The project relies on very specific sets of sequencing runs, including multiple (10+) ONT or PacBio samples paired with a large enough set of Illumina samples (40+)
	- Most available wildlife sequencing is either Illumina or PacBio/ONT, but rarely both
	- Honey bee data is limited to ONT, with no matched Illumina and no high-quality pangenome (although one is in the works elsewhere)
	- Water strider data from Arild Husby includes plenty of Illumina but only one PacBio HiFi sample
	- Bivalve data from Dan MacQueen includes both WGS and RNA-seq, but again only on PacBio HiFi sample
	- Atlantic salmon from Dan MacQueen includes high-quality ONT assemblies and WGS (unknown how much); data likely hard to access due to commercial interests
		- Also mentioned a promising rainbow trout dataset with potentially fewer restrictions
- Managing Cluster Transition
	- May have some delays while transitioning to new Baffin cluster (setting up environment etc.)
	- Timelines are unclear, hoping that Fiji stays operational for long enough to complete part 1
- Handling repeats
	- May need to run RepeatModeler/Masker to combat the effects of species-specific repeat regions which could systematically affect SV calls and gene annotations
	- Can explore outputs to see if this ends up being a serious problem; Jeon et al. 2026 didn't address this issue


#### Experimental Outline

##### Comparing pangenome methodologies

- Datasets
	- Steps
		- Identify candidate sequencing datasets that include long and short read samples, ideally from the same population
		- Targeting 10-20 long-read samples and ~50 short-read samples
	- Output: at least 2 strong candidates in addition to Apis mellifera; data downloaded and ready for processing
- Data prep & quality control
	- Tools: FastQC, NanoPlot
	- Inputs: raw Illumina short reads, raw ONT/PacBio long reads
	- Steps
		- QC Illumina samples; trim adapters and low-quality bases if needed
		- QC long-read reads/assemblies, filter by length or quality thresholds if needed
		- Document coverage and N50 per sample
	- Outputs: QC reports for seq data; clean Illumina reads and filtered long reads; can decide whether to proceed with dataset
- Assembly and SV Calling
	- Tools: Flye or HiFiasm, QUAST, BUSCO, Minimap2, Sniffles2, CuteSV, SURVIVOR. CD-HIT
	- Inputs: Long reads for each sample, conspecific reference genome
	- Steps
		- Run assembler on each sample separately (haplotype-aware if possible)
		- Run QUAST and BUSCO to QC the assemblies using the closest family in the BUSCO db
		- Run script to align and call variants against the conspecific reference, then merge and deduplicate variants to make the linear pangenome
		- Visually inspect this output to ensure repeats are not too prevalent
	- Outputs: Per-sample assemblies (fasta); assembly QC metrics; per-sample SV VCFs; merged SV catalog; linear pangenome incorporating SV sequence
- Minigraph-cactus graph
	- Tools: Cactus
	- Inputs: Assemblies and conspecific reference
	- Steps
		- Run Cactus script through singularity with all input assemblies, specifying the reference as the initial genome
	- Outputs: Pangenome graph with the conspecific reference as a backbone, ready for alignment and genotyping
- Short read alignment
	- Tools: BWA, VG Giraffe
	- Inputs: Illumina reads and all references (linear & graph)
	- Steps
		- Align reads to linear refs with BWA
		- Align reads to graph with Giraffe
		- Compute alignment rate, mean mapping quality, and mean coverage depth
	- Outputs: BAM and GAM alignments; alignment stats table
- Variant calling
	- Tools: GATK, VG
	- Inputs: BAM and GAM alignments
	- Steps
		- Call SNVs per sample per reference with GATK HaplotypeCaller
		- Call SNVs for graph with VG pack and call or project to BAM for GATK
		- Compute total SNV count, transition/transversion ratio (QC for bias), per-sample heterozygosity
	- Outputs: VCF variant calls per sample per reference; SNV stats table
- Evaluation Metrics
	- Allelic balance at heterozygous SNVs (ref bias)
		- Tools: BCFTools, VCFTools, MIXALIME? (see Buyan et al. 2025)
		- Steps
			- Extract heterozygous SNV calls across all samples per ref
			- Compute allelic depth ratio at each site (ref allele depth / total depth) and summarize (ref-biased will have ratio > 0.5)
			- Statistical comparison across methodologies
		- Outputs: Allelic balance distributions; summary stats; comparison figure
	- $F_{ST}$ (genetic differentiation) and Allele Frequency Spectrum
		- Tools: Biopython, ANGSD
		- Steps
			- Compute $F_{ST}$ score for the population (pairwise and between subgroups) for each reference
			- Compute AFS for each reference
		- Outputs: Table of $F_{ST}$ scores and AFS; comparison figure
	- Nucleotide diversity and Tajima's *D*
		- Tools: VCFTools (may need additional tools)
		- Steps
			- Do a sliding-window scan of at least one chromosome (ideally the largest) for nucleotide diversity
			- Optionally, compute Tajima's *D* for nonrandom mutation patterns
		- Outputs: Table of nucleotide diversity (and Tajima's *D*) for each sample per reference
	- $F_{ROH}$ (inbreeding coefficient)
		- Tools: BCFTools, [script](https://github.com/jyj5558/theta/blob/master/bin/ROHparser.py) (see Jeon et al. 2026)
		- Steps
			- Compute ROH per sample per ref
			- Compute $F_{ROH}$ = sum(len(ROH)) / len(autosome)
		- Outputs: Table of $F_{ROH}$ per sample per ref



#### Datasets
- https://www.ncbi.nlm.nih.gov/bioproject/PRJNA382404
	- Rhesus macaque
	- 64 pacbio samples from a single population, along with 1000+ illumina
	- Only project on SRA with 10+ long-read and short read data for the same species (other than human)
- Dan MacQueen
	- Atlantic salmon pangenome - 11 populations, lots of nanopore assemblies but may be too diverse for this
	- Proprietary salmon sequencing - great dataset with matched short and long samples, but would be a process to get access
	- Rainbow trout - with a Chinese group, similar to above but light on details so far, supposedly no access issues
- Arild Husby
	- Water strider - lots of short and a single long, probably not a great case unless we could supplement the data; Arild was more interested in sharing this one
	- House sparrow - 20 long-read with some Hi-C, but seemingly no short; they want to use PGGB
- Other interesting species with sufficient data:
	- Erithacus rubecula (robin) (CUTE!)
	- Salmo salar (Atlantic salmon) - plenty of public data, could supplement the MacQueen dataset
	- Acropora tenuis (stony corral)
	- Bos grunniens (yak)
	- Oreochromis niloticus (nile tilapia)
