---
name: viral_bio_database_retrieval
description: Retrieve viral genomes, proteins, and predicted protein structures from NCBI, UniProt, and AlphaFold using their APIs.
version: 1.0
category: bioinformatics

inputs:
  organism:
    type: string
    description: Name or taxonomy ID of the virus
  gene:
    type: string
    description: Gene name (example spike or S)

outputs:
  genome_sequence:
    type: string
  protein_records:
    type: json
  protein_structure_file:
    type: file

tools:
- http_request
- xml_parser
- json_parser
- fasta_parser

rate_limits:
  ncbi:
    without_api_key: 3_per_second
    with_api_key: 10_per_second

apis:

  ncbi:
    base_url: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/

    endpoints:

      search_genome:
        path: esearch.fcgi
        method: GET
        params:
          db: nucleotide
          term: "{{organism}} complete genome"
          retmode: json

      fetch_genome:
        path: efetch.fcgi
        method: GET
        params:
          db: nucleotide
          id: "{{genome_id}}"
          rettype: gb
          retmode: text

  uniprot:

    base_url: https://rest.uniprot.org

    endpoints:

      search_protein:
        path: /uniprotkb/search
        method: GET
        params:
          query: "organism_id:{{organism_tax_id}} AND gene:{{gene}} AND reviewed:true"
          format: json

      get_entry:
        path: /uniprotkb/{{accession}}
        method: GET

  alphafold:

    base_url: https://alphafold.ebi.ac.uk

    endpoints:

      get_structure:
        path: /files/AF-{{uniprot_id}}-F1-model_v4.pdb
        method: GET

workflow:

  step1:
    description: Search NCBI for viral genome
    action: ncbi.search_genome
    save: genome_id

  step2:
    description: Retrieve genome record
    action: ncbi.fetch_genome
    parse: genbank
    extract:
- sequence
- annotations

  step3:
    description: Search UniProt for protein entry
    action: uniprot.search_protein
    extract:
- accession
- protein_name
- function

  step4:
    description: Retrieve predicted protein structure
    action: alphafold.get_structure
    save: protein_structure_file

error_handling:
- retry_on_http_error
- handle_empty_results
- validate_api_response

attribution:
  ncbi: "Data from NCBI Entrez"
  uniprot: "Protein data from UniProt"
  alphafold: "Structure prediction from AlphaFold DB"
---

# Skill: Viral Bio Database Retrieval

## Purpose

This skill enables an AI agent to retrieve viral biological data from major scientific databases using APIs.

The agent can:

- retrieve viral genomes
- identify proteins
- download predicted protein structures

## Supported Databases

### NCBI
Used for retrieving viral genome sequences.

### UniProt
Used for curated protein entries and annotations.

### AlphaFold
Used for predicted 3D protein structures.

## Example Task

Input:

organism: SARS-CoV-2  
gene: S

Expected Result:

- Wuhan reference genome
- spike protein entry
- AlphaFold predicted spike structure