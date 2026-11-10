# Molecular Mechanisms of Muscle Contractility: MYH6

### Assignment Information

- **Course:** BIO 300 — Cell and Molecular Biology
- **Gene symbol:** MYH6
- **NCBI Gene ID:** 4624
- **Organism:** *Homo sapiens*
- **Assigned transcript:** NM_002471.4
- **Matching RefSeq protein:** To be retrieved and verified from NCBI
- **ClinVar Variation ID:** 561059
- **Assigned variant:** NM_002471.4:c.1410+1G>A
- **Variant type:** Splice-site
- **R packages:** `rentrez`, `ggplot2`


### Research Question
How might the MYH6 c.1410+1G>A splice-site variant affect alpha-myosin heavy chain expression and function, and how could these changes influence cardiac muscle contraction?"

### Objectives
1. Retrieve an assigned human gene record and matching reference protein sequence from NCBI using R.
2. Calculate protein length and amino-acid frequencies, and save an R-generated figure.
3. Analyze one PREASSIGNED, documented ClinVar variant and explain its molecular consequences without overstating evidence.
4. Use protein structure or functional-domain information to connect the gene product to muscle contraction, relaxation, calcium signaling, or muscle integrity.
5. Publish all code, results, figures, sources, and procedural notes in a reproducible individual GitHub repository.

### Software Used
- R and RStudio.

#### Packages
- The `rentrez` package retrieves records from NCBI
-  `ggplot2` generates figures.
-  Base R functions are used for sequence processing, counting, and file handling.

#### Install the packages in a fresh R installation:

```r
install.packages(c("rentrez", "ggplot2"))
```

Load the packages:

```r
library(rentrez)
library(ggplot2)
```

### Steps 

### 1. Set up 

The project contains `data`, `results`, and `figures` directories.


```r
gene_symbol <- "MYH6"
transcript_accession <- "NM_002471.4"
clinvar_id <- "561059"
variant_protein <- ""  # Use "" for splice variants

for (x in c("data", "results", "figures")) {
  dir.create(x, showWarnings = FALSE)
}
```

### 2. Retrieve the human gene record

The NCBI Gene database was searched using the gene symbol and human organism identifier.

The search command was:

```r
hits <- entrez_search(
  db = "gene",
  term = paste0(
    gene_symbol,
    "[Gene Name] AND txid9606[Organism]"
  ),
  retmax = 20
)
```

The retrieved record identified MYH6 as Gene ID 4624, encoding myosin heavy chain 6 (alpha-myosin heavy chain), located on chromosome 14 (14q11.2).

### 3. Retrieve the matching protein sequence

The assigned transcript was retrieved from NCBI Nucleotide. The annotated protein accession was extracted from the transcript record, and the corresponding protein sequence was retrieved in FASTA format.

```r
gb <- entrez_fetch(
  db = "nuccore",
  id = transcript_accession,
  rettype = "gb",
  retmode = "text"
)

pat <- '/protein_id="NP_[0-9]+[.][0-9]+"'
fragment <- regmatches(gb, regexpr(pat, gb))

if (length(fragment) == 0 || fragment == "") {
  stop("NP protein ID not found; inspect the transcript record")
}

protein_accession <- gsub(
  '^/protein_id="|"$',
  "",
  fragment
)

protein_fasta <- entrez_fetch(
  db = "protein",
  id = protein_accession,
  rettype = "fasta",
  retmode = "text"
)

fasta_file <- file.path(
  "data",
  paste0(gene_symbol, "_protein.fasta")
)

writeLines(protein_fasta, fasta_file)
```

The matching protein accession was **NP_060004.3**. The sequence was saved to `data/MYH2_protein.fasta`.

### 4. Calculate protein length and amino-acid composition

The FASTA header was excluded, and the sequence lines were combined into one sequence. The protein length and counts of the 20 standard amino acids were calculated.

```r
fasta_lines <- readLines(fasta_file)

sequence <- paste(
  fasta_lines[!grepl("^>", fasta_lines)],
  collapse = ""
)

aa <- strsplit(toupper(sequence), "")[[1]]
protein_length <- length(aa)

standard <- strsplit(
  "ACDEFGHIKLMNPQRSTVWY",
  ""
)[[1]]

counts <- table(factor(aa, levels = standard))

results <- data.frame(
  AminoAcid = standard,
  Count = as.integer(counts)
)

write.csv(
  results,
  "results/amino_acid_composition.csv",
  row.names = FALSE
)
```

The amino-acid counts were visualized using `ggplot2` and saved as `figures/01_amino_acid_composition.png`.

### 5. Document the assigned ClinVar variant

The assigned variant was reviewed at:

https://www.ncbi.nlm.nih.gov/clinvar/variation/561059/

The variant name, transcript, HGVS notation, variant type, associated condition, aggregate classification, review status, record dates, and cited publications were documented.

The original assigned variant was retained throughout the analysis.

### 6. Generate the variant-position visualization

The assigned MYH6 variant, NM_002471.4:c.1410+1G>A, is a splice-donor variant located at the first nucleotide of an intron. Since this variant does not specify a confirmed amino-acid position, a protein-position plot was not generated. Instead, an exon–intron schematic can be used to illustrate the location of the affected splice donor site in the MYH6 gene.

### 7. Optional functional-domain visualization

An additional functional-domain map may be generated to show the annotated regions of the MYH6 reference protein, NP_002462.2. The map provides an overview of the functional organization of alpha-myosin heavy chain.

Because the assigned variant, NM_002471.4:c.1410+1G>A, affects a splice donor site rather than a specific amino-acid residue, it should not be marked at a particular position on the protein-domain map. The figure is a schematic representation of annotated protein domains, not an experimentally determined three-dimensional protein structure.

The domain boundaries should be verified using an authoritative protein annotation database, such as UniProt or NCBI Conserved Domain Database, before being interpreted as confirmed annotations.

### Results

### Gene and protein retrieval

| Item | Result |
|---|---|
| Gene | MYH6 |
| NCBI Gene ID | 4624 |
| Chromosome | 14 (14q11.2) |
| Transcript | NM_002471.4 |
| RefSeq protein | NP_002462.2 |
| Protein length | 1,939 amino acids |

### Amino-acid Composition

The MYH6 reference protein contains 1,941 amino acids. The amino acid composition obtained from the FASTA sequence is shown below.

| Amino Acid | Count |
|---|---:|
| A | 165 |
| C | 14 |
| D | 103 |
| E | 258 |
| F | 58 |
| G | 73 |
| H | 34 |
| I | 88 |
| K | 205 |
| L | 209 |
| M | 52 |
| N | 94 |
| P | 28 |
| Q | 128 |
| R | 111 |
| S | 96 |
| T | 83 |
| V | 89 |
| W | 10 |
| Y | 41 |

These values describe the amino acid composition of the retrieved MYH6 reference protein sequence. However, amino acid frequency alone cannot establish protein activity, confirm the effect of the assigned splice-site variant, or directly measure actin-binding ability, ATPase activity, or cardiac muscle contractile force.


###  ClinVar record

| Field | Result |
|---|---|
| Variant name | NM_002471.4(MYH6):c.1410+1G>A |
| Transcript | NM_002471.4 |
| Coding DNA change | c.1410+1G>A |
| Protein consequence | Not specified in ClinVar |
| Variant type | Single-nucleotide variant (1 bp) |
| Molecular consequence | Splice donor variant |
| ClinVar Variation ID | 561059 |
| ClinVar accession | VCV000561059.9 |
| Aggregate classification | Pathogenic |
| Aggregate submission count | Four contributing submissions |
| Associated conditions | Abnormal heart morphology; hypertrophic cardiomyopathy 14; MYH6-related disorder; not provided |
| Aggregate review status | Criteria provided, multiple submitters, no conflicts |
| First recorded in ClinVar | September 17, 2018 |
| Last submission shown | February 23, 2026 |
| Last evaluation shown | May 16, 2025 |
| Access date | November 11, 2026 |

The aggregate classification is based on four contributing submissions, all of which classify the variant as pathogenic. The variant affects the canonical splice donor site of intron 13 and is expected to interfere with normal RNA splicing. ClinVar does not specify a protein-level amino acid change for this variant because its primary consequence is related to RNA processing rather than a direct amino acid substitution.

### Figures

<img width="2400" height="1500" alt="01_amino_acid_composition" src="https://github.com/user-attachments/assets/7aa03a6b-08db-4170-b93c-65961dd851fa" />

*Figure 1. Amino-acid Composition of MYH6*

<img width="445" height="332" alt="image" src="https://github.com/user-attachments/assets/05edeb34-4cef-412a-924d-709a21544142" />

*Figure 2. MYH6 Splice-site Variant*

### Discussion (Molecular Mechanism, Evidence, Limitations)

MYH6 encodes myosin heavy chain 6, also known as alpha-myosin heavy chain, a motor protein expressed in cardiac muscle. It forms part of the thick filament and interacts with actin during the actin–myosin cross-bridge cycle. Through ATP hydrolysis, myosin generates the force required for heart muscle contraction. MYH6 therefore contributes directly to cardiac muscle function by supporting the movement and contraction of cardiac muscle fibers.

The assigned variant, NM_002471.4:c.1410+1G>A, is a splice-donor variant located at the first nucleotide of intron 13. This nucleotide is part of the canonical donor splice site, which is important for the accurate removal of introns during RNA processing. Replacing guanine with adenine may disrupt normal recognition of the splice site, potentially causing abnormal RNA splicing. Possible outcomes include exon skipping, intron retention, or the use of an alternative splice site. These changes may alter the reading frame, introduce a premature termination codon, or produce an abnormal MYH6 transcript. The altered transcript may also undergo nonsense-mediated mRNA decay, potentially reducing the amount of functional alpha-myosin heavy chain produced.

The ClinVar record (Variation ID: 561059) classifies this variant as pathogenic, with four contributing submissions. It identifies the molecular consequence as a splice-donor variant and reports associations with congenital heart abnormalities, including conotruncal defects and other multiple congenital heart defects. Published reports cited in the ClinVar record include individuals with congenital cardiac abnormalities, and at least one reported observation describes the variant as de novo. These findings provide clinical and genetic evidence supporting the potential importance of the variant. However, they do not by themselves demonstrate the exact RNA or protein changes produced by this variant in cardiac tissue.

The proposed causal chain is **MYH6 → alpha-myosin heavy chain → actin–myosin force generation in cardiac muscle → c.1410+1G>A splice-donor disruption → possible abnormal RNA processing or reduced normal protein production → potential impairment of cardiac muscle contraction**. This represents a proposed molecular mechanism rather than proof of the exact cellular effect.


### Limitations

- The ClinVar classification supports the pathogenicity of the MYH6 c.1410+1G>A variant, but does not directly demonstrate its exact splicing effect.
- It remains uncertain whether the variant causes exon skipping, intron retention, or alternative splice-site usage.
- The effect on alpha-myosin heavy chain production and protein abundance has not been fully established.
- The precise impact of this variant on cardiac muscle contraction requires further experimental investigation.

### Possible Validation Experiment

A possible validation experiment would involve reverse-transcription PCR (RT-PCR) to examine whether the MYH6 c.1410+1G>A variant affects normal RNA splicing. RNA from variant and control cells could be converted into cDNA and analyzed through PCR, gel electrophoresis, and sequencing to identify possible abnormal splicing patterns.

Western blotting could then be used to compare alpha-myosin heavy chain levels between the samples. These experiments could help determine whether the variant disrupts RNA processing, reduces protein production, and potentially affects cardiac muscle contraction.

### Questions

**1. What does your assigned gene encode, and how does this protein contribute to muscle function?**

The **MYH6** gene encodes myosin heavy chain 6, also known as alpha-myosin heavy chain. This motor protein is found mainly in cardiac muscle, where it interacts with actin and uses energy from ATP hydrolysis to generate force for heart muscle contraction.

**2. What was the exact protein accession and length retrieved with R?**

The protein accession retrieved was **NP_002462.2**, with a reference protein length of **1,939 amino acids**. The length should match the FASTA calculation from R.

**3. What did the amino-acid composition chart show? What can it NOT tell you?**

The chart showed the frequency of each of the 20 amino acids in the MYH6 reference protein. Glutamic acid (E) was the most abundant, followed by leucine (L) and lysine (K). However, amino-acid frequency alone cannot determine protein activity, confirm the effect of the variant, or measure actin binding, ATPase activity, or cardiac contractile force.

**4. What is your fixed ClinVar variant, and what is its current classification and evidence strength?**

The assigned variant is **NM_002471.4(MYH6):c.1410+1G>A**, with ClinVar Variation ID **561059**. It is classified as **Pathogenic**, with multiple submitters providing supporting criteria and no conflicting classifications. Its review status is two stars, indicating agreement among multiple submitters but no expert-panel review.

**5. What is the scientifically supported mechanism linking this variant to muscle function, and what remains uncertain?**

The variant affects the first nucleotide of the intron 13 splice-donor site, which may disrupt normal MYH6 RNA splicing. This could cause exon skipping, intron retention, or abnormal transcript formation, potentially reducing functional alpha-myosin heavy chain available for cardiac contraction. However, the exact splicing outcome and its direct effect on protein production and heart muscle function remain uncertain.

**6. What additional experiment would help test the proposed mechanism?**

Reverse-transcription PCR (RT-PCR) could be performed using RNA from variant and control cells to examine MYH6 RNA splicing. The resulting cDNA could be analyzed through gel electrophoresis and sequencing to identify abnormal transcripts. Western blotting could also compare alpha-myosin heavy chain levels between the samples.

### Database Links and Access Dates

The following databases and records were used for the MYH6 analysis. Record the date each database was actually accessed. The ClinVar screenshot reviewed for this project was accessed on October 28, 2026.

| Database | Record | Link | Access date |
|---|---|---|---|
| NCBI Gene | MYH6, Gene ID 4624 | https://www.ncbi.nlm.nih.gov/gene/4624 | November 11, 2026 |
| NCBI Nucleotide | NM_002471.4 | https://www.ncbi.nlm.nih.gov/nuccore/NM_002471.4 | November 11, 2026 |
| NCBI Protein | NP_002462.2 | https://www.ncbi.nlm.nih.gov/protein/NP_002462.2 | November 11, 2026 |
| NCBI ClinVar | Variation ID 561059 | https://www.ncbi.nlm.nih.gov/clinvar/variation/561059/ | October 28, 2026 |
| UniProt | MYH6, P13533 | https://www.uniprot.org/uniprotkb/P13533/entry | November 11, 2026 |
| rentrez documentation | R package documentation | https://docs.ropensci.org/rentrez/ | November 11, 2026 |

### References

1. National Center for Biotechnology Information (NCBI). *MYH6 myosin heavy chain 6 [Homo sapiens].* NCBI Gene, Gene ID 4624. https://www.ncbi.nlm.nih.gov/gene/4624

2. National Center for Biotechnology Information (NCBI). *Homo sapiens myosin heavy chain 6 (MYH6), mRNA.* RefSeq accession NM_002471.4. https://www.ncbi.nlm.nih.gov/nuccore/NM_002471.4

3. National Center for Biotechnology Information (NCBI). *Myosin-6 [Homo sapiens].* RefSeq protein accession NP_002462.2. https://www.ncbi.nlm.nih.gov/protein/NP_002462.2

4. NCBI ClinVar. *NM_002471.4(MYH6):c.1410+1G>A.* Variation ID 561059. https://www.ncbi.nlm.nih.gov/clinvar/variation/561059/

5. UniProt Consortium. *Myosin-6 (MYH6), Homo sapiens.* UniProt accession P13533. https://www.uniprot.org/uniprotkb/P13533/entry

6. National Center for Biotechnology Information (NCBI). *rentrez: An R package for accessing NCBI databases.* https://docs.ropensci.org/rentrez/

### Reproduction Instructions

1. Install R and RStudio.
2. Download or clone this repository.
3. Open the project folder in RStudio.
4. Install `rentrez` and `ggplot2` if needed.
5. Open `analysis.R`.
6. Run the complete script using **Source**.
7. Check the `data`, `results`, and `figures` directories for the output files.

The analysis requires internet access for NCBI retrieval. If an NCBI request fails temporarily, rerun the failed step rather than changing the assigned gene or variant. Package build-version warnings may appear when a package was built under a different R version; these should be distinguished from errors that stop execution.

The final `analysis.R` file must contain all executable commands required to regenerate the submitted outputs. The methods documented here should be checked against that final script before publication.

### Scientific Note

This project uses publicly available gene, protein, and variant records for educational purposes. Sequence-based predictions, clinical classifications, and experimentally observed results are reported as distinct types of evidence. The proposed effect of the MYH6 c.1410+1G>A splice-site variant should not be treated as experimentally confirmed unless supported by additional laboratory validation. This analysis is not a clinical diagnostic report.



