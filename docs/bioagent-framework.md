# Viral Bioinformatics AI Agent Skill Library

This document contains a complete **skill library and orchestration
pipeline** for a bioinformatics AI agent capable of performing viral
genomics analysis.

The skills included:

-   viral_genome_fetch
-   viral_sequence_analysis_workflows
-   mutation_detection
-   protein_binding_analysis
-   protein_structure_fetch
-   drug_resistance_prediction
-   viral_analysis_pipeline (meta-skill orchestrator)

------------------------------------------------------------------------

# Skill: viral_genome_fetch

## Description

Retrieve viral genome sequences from NCBI using the E-utilities API.

### Inputs

-   organism (string)

### Outputs

-   genome_sequence (FASTA)
-   genome_metadata (JSON)

### Workflow

1.  Search for genome:

```{=html}
<!-- -->
```
    https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term={organism}+complete+genome

2.  Fetch genome record:

```{=html}
<!-- -->
```
    https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={genome_id}&rettype=fasta

------------------------------------------------------------------------

# Skill: viral_sequence_analysis_workflows

## Description

Perform core viral sequence analysis including multiple sequence
alignment and phylogenetic tree construction.

### Inputs

-   fasta_sequences
-   reference_sequence

### Outputs

-   alignment_file
-   phylogenetic_tree
-   conservation_scores

## Workflow 1: Multiple Sequence Alignment

Tools: - MAFFT - ClustalOmega - MUSCLE

Example command:

    mafft --auto --thread -1 input.fasta > output.aln

Key parameters:

-   `--auto` select best algorithm automatically
-   `--thread -1` use all CPU cores
-   `--maxiterate 1000` high accuracy mode

Post-processing:

-   convert to PHYLIP format
-   compute sequence identity matrix
-   visualize conservation scores

## Workflow 2: Variant Detection

Steps:

1.  Load alignment with Biopython
2.  Identify variable columns
3.  Classify mutation types
4.  Calculate allele frequencies
5.  Map to reference genome coordinates

## Workflow 3: Phylogenetic Tree Construction

Tools:

-   IQ-TREE
-   RAxML
-   FastTree

Pipeline:

1.  Clean alignment with trimal
2.  Model selection using IQ-TREE MFP
3.  Maximum likelihood tree inference
4.  Bootstrap analysis (1000 replicates)
5.  Tree rooting
6.  Metadata annotation

------------------------------------------------------------------------

# Skill: mutation_detection

## Description

Detect nucleotide mutations between reference and query genomes.

### Inputs

-   reference_genome
-   query_genome

### Outputs

-   mutation_list
-   mutation_summary

### Steps

1.  Align sequences
2.  Identify variant positions
3.  Classify mutations

Mutation types:

-   synonymous
-   missense
-   nonsense
-   frameshift

------------------------------------------------------------------------

# Skill: protein_binding_analysis

## Description

Determine whether mutations affect protein functional domains or binding
regions.

### Inputs

-   mutations
-   protein_annotations

### Outputs

-   binding_site_impacts

### Process

1.  Map mutations to proteins
2.  Check functional domains
3.  Identify potential binding site disruptions

------------------------------------------------------------------------

# Skill: protein_structure_fetch

## Description

Retrieve predicted protein structures from AlphaFold.

### Input

-   uniprot_id

### Output

-   pdb_structure

### API

Structure URL pattern:

    https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb

------------------------------------------------------------------------

# Skill: drug_resistance_prediction

## Description

Predict antiviral drug resistance based on observed mutations.

### Inputs

-   mutations
-   protein_structure

### Outputs

-   resistance_risk
-   affected_drugs

### Process

1.  Map mutations to drug binding sites
2.  Evaluate structural changes
3.  Predict drug efficacy changes

------------------------------------------------------------------------

# Meta Skill: viral_analysis_pipeline

## Description

End-to-end automated viral genomics analysis pipeline.

### Inputs

-   organism
-   gene
-   query_sequences

### Outputs

-   genome_data
-   alignment
-   mutation_report
-   phylogenetic_tree
-   binding_analysis
-   drug_resistance_assessment

### Execution Pipeline

    viral_genome_fetch
            ↓
    viral_sequence_analysis_workflows
            ↓
    mutation_detection
            ↓
    protein_binding_analysis
            ↓
    protein_structure_fetch
            ↓
    drug_resistance_prediction

------------------------------------------------------------------------

# Example Usage

User request:

    Analyze new SARS-CoV-2 variant sequences

Agent execution:

1.  Fetch reference genome
2.  Align sequences
3.  Detect mutations
4.  Build phylogenetic tree
5.  Analyze protein binding effects
6.  Predict drug resistance

Output example:

-   12 mutations detected
-   3 mutations in spike receptor binding domain
-   AlphaFold structure retrieved
-   potential reduced antibody binding

------------------------------------------------------------------------

# Notes

Best practices for viral phylogenetics:

-   Use time-aware evolutionary models for rapidly evolving viruses
-   Account for recombination in segmented viruses like influenza
-   Consider sampling bias when interpreting phylogenetic trees
