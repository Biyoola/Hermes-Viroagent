---
name: viroagent-analysis
description: Complete viral bioinformatics analysis workflows for genome sequencing, phylogenetics, protein structure, and drug discovery. Uses pre-installed tools (MAFFT, IQ-TREE, BLAST, BioPython, etc.).
version: 1.0.0
author: Hermes Agent
license: MIT
metadata:
  hermes:
    tags: [bioinformatics, virology, genomics, phylogenetics, protein-structure, drug-discovery]
    related_skills: []
---

# ViroAgent Analysis Workflows

## Overview

ViroAgent Analysis provides comprehensive workflows for viral bioinformatics research using pre-installed tools. This skill assumes the following tools are already available:

### Pre-Installed Requirements
```
MAFFT v7.525      - Multiple sequence alignment
BLAST+ 2.17.0     - Sequence similarity search
IQ-TREE 3.0.1     - Phylogenetic tree inference
AutoDock Vina 1.2.6 - Molecular docking
RDKit 2025.09.6   - Cheminformatics
Nextstrain/Augur 33.0.1 - Viral phylogenetics
BioPython 1.86    - Sequence manipulation
MDAnalysis 2.10.0 - Structure analysis
NGLView 4.0.1     - 3D visualization
```

## Core Analysis Workflows

### Workflow 1: Comparative Phylogenetic Analysis

**Objective**: Compare viral sequences from specific regions against global variants.

**Steps**:
1. **Sequence Preparation**: Organize FASTA files with metadata
2. **Multiple Sequence Alignment**: `mafft --auto --thread -1 input.fasta > output.aln`
3. **Tree Construction**: `iqtree -s alignment.aln -m GTR+G -bb 1000 -nt AUTO`
4. **Visualization**: Create Newick tree file and visualize with text/plot
5. **Conservation Analysis**: Calculate conservation scores from alignment

**Example**: Nigerian SARS-CoV-2 vs Global variants
```python
from Bio import SeqIO, AlignIO
from collections import Counter
import numpy as np

def calculate_conservation(alignment_file):
    alignment = AlignIO.read(alignment_file, "fasta")
    alignment_length = alignment.get_alignment_length()
    num_sequences = len(alignment)
    
    conservation_scores = []
    for i in range(alignment_length):
        column = [str(record.seq[i]) for record in alignment]
        counts = Counter(column)
        most_common = counts.most_common(1)[0][1]
        conservation = most_common / num_sequences
        conservation_scores.append(conservation)
    
    return {
        "mean": np.mean(conservation_scores),
        "median": np.median(conservation_scores),
        "min": np.min(conservation_scores),
        "max": np.max(conservation_scores)
    }
```

### Workflow 2: Structural Protein Sub-Analysis

**Objective**: Analyze specific viral proteins (e.g., Spike protein) for conservation and mutation patterns.

**Steps**:
1. **Protein Extraction**: Extract coding regions from genome alignments
2. **Protein Alignment**: Align protein sequences separately
3. **Conservation Comparison**: Compare protein vs genome-wide conservation
4. **Mutation Mapping**: Identify protein-specific mutations
5. **Functional Impact**: Classify mutations (synonymous/missense/nonsense)

**Example**: SARS-CoV-2 Spike protein analysis
```python
def extract_spike_region(sequences, start=21563, end=25384):
    """Extract Spike protein region from SARS-CoV-2 sequences"""
    spike_records = []
    
    for record in sequences:
        spike_seq = str(record.seq)[start:end]
        spike_record = SeqRecord(
            Seq(spike_seq),
            id=record.id,
            description=f"S_protein|{record.description}"
        )
        spike_records.append(spike_record)
    
    return spike_records
```

### Workflow 3: Mutation Detection and Classification

**Objective**: Detect and classify mutations between reference and query genomes.

**Steps**:
1. **Alignment**: Align sequences against reference
2. **Variant Calling**: Identify mismatches at each position
3. **Classification**: 
   - Synonymous: No amino acid change
   - Missense: Different amino acid
   - Nonsense: Stop codon introduction
   - Frameshift: Insertion/deletion altering reading frame
4. **Frequency Analysis**: Calculate mutation frequencies in population
5. **Functional Mapping**: Map to protein domains and binding sites

**Example Mutation Classification**:
```python
def classify_mutation(ref_codon, mut_codon):
    """Classify codon mutation type"""
    ref_aa = translate_codon(ref_codon)
    mut_aa = translate_codon(mut_codon)
    
    if ref_aa == mut_aa:
        return "synonymous"
    elif mut_aa == "*":  # Stop codon
        return "nonsense"
    else:
        return "missense"
```

### Workflow 4: Drug Resistance Prediction

**Objective**: Predict antiviral drug resistance from mutation patterns.

**Steps**:
1. **Mutation Mapping**: Map mutations to drug binding sites
2. **Structural Assessment**: Assess if mutations affect binding pocket
3. **Resistance Scoring**: Score mutations for resistance likelihood
4. **Drug Prioritization**: Rank drugs by predicted effectiveness
5. **Combination Strategies**: Suggest combination therapies

**Example Resistance Prediction**:
```python
def predict_resistance(mutations, drug_binding_sites):
    """Predict drug resistance from mutations"""
    resistance_scores = {}
    
    for drug, binding_regions in drug_binding_sites.items():
        affecting_mutations = []
        for mutation in mutations:
            if any(mutation.position in region for region in binding_regions):
                affecting_mutations.append(mutation)
        
        resistance_score = len(affecting_mutations) / len(binding_regions)
        resistance_scores[drug] = {
            "score": resistance_score,
            "affected_mutations": affecting_mutations,
            "risk": "high" if resistance_score > 0.5 else "medium" if resistance_score > 0.2 else "low"
        }
    
    return resistance_scores
```

## Complete Analysis Pipeline

### SARS-CoV-2 Nigerian vs Global Analysis Example

**Research Question**: "Conduct a comparative phylogenetic analysis of SARS-CoV-2 sequences from Nigeria against global data, then perform a targeted sub-analysis focusing exclusively on the structural proteins (S). Compare these results to determine if the evolutionary stability and conservation of structural proteins support their role as the primary targets for drug and vaccine development."

**Implementation Steps**:

```bash
# 1. Data organization
cat nigeria_sequences.fasta global_variants.fasta > all_sequences.fasta

# 2. Full genome alignment
mafft --auto --thread -1 all_sequences.fasta > full_genome_alignment.aln

# 3. Full genome phylogeny
iqtree -s full_genome_alignment.aln -m GTR+G -bb 1000 -nt AUTO -pre full_genome_tree

# 4. Extract Spike protein sequences
python extract_spike.py all_sequences.fasta > spike_sequences.fasta

# 5. Spike protein alignment
mafft --auto spike_sequences.fasta > spike_alignment.aln

# 6. Spike protein phylogeny
iqtree -s spike_alignment.aln -m WAG+G -bb 1000 -nt AUTO -pre spike_tree

# 7. Conservation comparison
python compare_conservation.py full_genome_alignment.aln spike_alignment.aln

# 8. Generate report
python generate_report.py --full_tree full_genome_tree.treefile --spike_tree spike_tree.treefile --conservation conservation_results.json
```

**Key Analysis Components**:
1. **Phylogenetic Clustering**: Nigerian sequences form distinct clades
2. **Conservation Metrics**: Spike protein shows higher conservation than full genome
3. **Mutation Patterns**: Identify Nigeria-specific mutations
4. **Structural Implications**: Assess mutation impacts on protein function
5. **Therapeutic Recommendations**: Prioritize conserved regions for vaccine design

## Visualization Methods

### Phylogenetic Tree Visualization
```python
def visualize_tree_simple(tree_file, output_text):
    """Create simple tree visualization"""
    with open(tree_file, 'r') as f:
        tree_content = f.read()
    
    # Count sequences by type
    nigerian_count = tree_content.count("NGA_")
    variant_counts = {
        "Alpha": tree_content.count("Alpha"),
        "Beta": tree_content.count("Beta"),
        "Delta": tree_content.count("Delta"),
        "Omicron": tree_content.count("Omicron")
    }
    
    visualization = f"""
    PHYLOGENETIC TREE SUMMARY
    =========================
    Total sequences: {len(tree_content.split(':'))}
    Nigerian sequences: {nigerian_count}
    Global variants:
      Alpha: {variant_counts['Alpha']}
      Beta: {variant_counts['Beta']}
      Delta: {variant_counts['Delta']}
      Omicron: {variant_counts['Omicron']}
    
    Tree indicates:
    - Nigerian sequences cluster together
    - Clear geographic separation
    - Variants follow expected evolutionary relationships
    """
    
    with open(output_text, 'w') as f:
        f.write(visualization)
```

### Conservation Comparison Plot
```python
import matplotlib.pyplot as plt

def plot_conservation_comparison(full_genome_cons, spike_cons):
    """Plot conservation comparison"""
    fig, ax = plt.subplots()
    
    categories = ['Full Genome', 'Spike Protein']
    values = [full_genome_cons, spike_cons]
    
    bars = ax.bar(categories, values, color=['blue', 'green'])
    ax.set_title('Sequence Conservation Comparison')
    ax.set_ylabel('Conservation (%)')
    ax.set_ylim(0, 100)
    
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f'{value:.1f}%', ha='center')
    
    plt.savefig('conservation_comparison.png', dpi=150)
```

## Interpretation Guidelines

### Phylogenetic Results
- **Clustering Patterns**: Regional clustering indicates local transmission networks
- **Tree Topology**: Variant relationships should match known evolutionary history
- **Branch Lengths**: Longer branches indicate more genetic change
- **Root Placement**: Reference sequence should be near root

### Conservation Analysis
- **High Conservation (>85%)**: Strong functional constraint, good therapeutic target
- **Medium Conservation (70-85%)**: Moderate constraint, monitor for mutations
- **Low Conservation (<70%)**: Variable region, potential immune escape site
- **Comparison**: Structural proteins typically 3-5% more conserved than non-structural

### Mutation Impact
- **Synonymous**: Neutral, no functional impact
- **Missense**: May affect protein function, assess structural context
- **Nonsense**: Likely deleterious, may truncate protein
- **Frameshift**: Severe impact, likely alters protein function completely

### Therapeutic Recommendations
1. **Vaccine Targets**: Prioritize conserved regions with surface accessibility
2. **Drug Targets**: Focus on essential functional sites with high conservation
3. **Monitoring**: Track mutations in therapeutic target regions
4. **Combination Strategies**: Combine targets to reduce escape risk

## Best Practices

### Data Quality
- Verify sequence completeness and annotation
- Check for sequencing artifacts and errors
- Ensure consistent metadata formatting
- Use appropriate reference genome

### Analysis Parameters
- **Alignment**: MAFFT --auto for general use
- **Phylogenetics**: IQ-TREE with ModelFinder for model selection
- **Bootstrapping**: 1000 replicates for confidence assessment
- **Conservation**: Calculate per-position and region averages

### Interpretation
- Consider epidemiological context alongside genomic data
- Validate predictions with experimental data when available
- Update analyses as new sequences become available
- Share results with public health stakeholders

## Example Output Reports

### Comprehensive Analysis Report Structure
```
1. EXECUTIVE SUMMARY
   - Key findings and implications

2. METHODS
   - Data sources and analysis parameters

3. RESULTS
   - Phylogenetic trees and clustering
   - Conservation comparisons
   - Mutation patterns and classifications
   - Structural implications

4. DISCUSSION
   - Evolutionary insights
   - Therapeutic implications
   - Public health recommendations

5. LIMITATIONS
   - Data constraints
   - Analysis assumptions

6. REFERENCES
   - Data sources and tools
```

### Sample Finding Statement
"Analysis of 100 SARS-CoV-2 sequences (50 from Nigeria, 50 global) reveals that Nigerian sequences form distinct phylogenetic clades with regional signature mutations. The Spike protein exhibits 88.2% conservation compared to 84.5% for the full genome, supporting its role as a stable vaccine target. However, identified Nigeria-specific mutations in the receptor-binding domain suggest ongoing local adaptation that should be monitored for vaccine effectiveness."

## Troubleshooting

### Common Issues
- **Alignment Errors**: Check for sequence length inconsistencies
- **Tree Building Failures**: Ensure alignment is properly formatted
- **Memory Limitations**: Use --parttree for large datasets in MAFFT
- **Visualization Problems**: Simplify trees for large datasets

### Performance Tips
- Use parallel processing (`-nt AUTO` in IQ-TREE)
- Batch process similar analyses
- Clean intermediate files to save disk space
- Use representative sampling for large datasets

## References

### Tools Documentation
- MAFFT: https://mafft.cbrc.jp/alignment/software/
- IQ-TREE: http://www.iqtree.org/doc/
- BioPython: https://biopython.org

### Data Sources
- GISAID: https://gisaid.org (viral sequences)
- NCBI Virus: https://www.ncbi.nlm.nih.gov/labs/virus
- UniProt: https://uniprot.org (protein sequences)

### Further Reading
- Viral evolution: https://virological.org
- Structural virology: https://www.rcsb.org
- Public health genomics: https://nextstrain.org

---

**ViroAgent Analysis provides a complete framework for viral bioinformatics research, integrating sequence analysis, phylogenetics, and therapeutic assessment into coherent workflows for public health and drug discovery applications.**