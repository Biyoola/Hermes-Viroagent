# Skill: Sequence Analysis Pipeline

## Overview

This skill enables the AI agent to perform comprehensive sequence analysis workflows including multiple sequence alignment (MSA), variant calling, and phylogenetic tree construction. It provides the foundation for comparative genomics, mutation detection, and evolutionary analysis of viral genomes.

---

## Capabilities

- Retrieve and manipulate sequences from biological databases
- Perform multiple sequence alignment (MSA) using state-of-the-art algorithms
- Identify variants and mutations from alignments
- Construct and interpret phylogenetic trees
- Calculate sequence diversity and conservation metrics
- Prepare data for downstream structural and functional analysis

---

## Prerequisites

### Required Software

| Tool | Installation | Purpose |
|------|--------------|---------|
| **MAFFT** | `conda install -c bioconda mafft` | Multiple sequence alignment |
| **ClustalOmega** | `conda install -c bioconda clustalo` | Alternative MSA (large datasets) |
| **MUSCLE** | `conda install -c bioconda muscle` | Alternative MSA (high accuracy) |
| **BLAST+** | `conda install -c bioconda blast` | Sequence similarity search |
| **IQ-TREE** | `conda install -c bioconda iqtree` | Phylogenetic inference |
| **RAxML** | `conda install -c bioconda raxml` | Alternative phylogenetics |
| **FastTree** | `conda install -c bioconda fasttree` | Fast approximate phylogeny |
| **trimal** | `conda install -c bioconda trimal` | Alignment trimming |
| **Biopython** | `pip install biopython` | Sequence manipulation |
| **scikit-bio** | `pip install scikit-bio` | Diversity analysis |
| **DendroPy** | `pip install dendropy` | Phylogenetic tree handling |
| **Nextstrain/Augur** | `conda install -c bioconda augur` | Viral phylogenetics |

### Required Data Access

- **NCBI E-utilities**: For GenBank sequence retrieval
- **NCBI Datasets API**: Modern REST API for bulk downloads
- **GISAID**: For SARS-CoV-2 and influenza sequences (requires application)
- **UniProt**: For protein sequence retrieval

### Hardware Requirements

- **Minimum**: Standard CPU, 8GB RAM (small datasets <100 sequences)
- **Recommended**: Multi-core CPU, 16GB RAM (medium datasets 100-1000 sequences)
- **Optimal**: High-core count CPU, 32GB+ RAM, SSD storage (large datasets)

---

## Core Methods

### 1. Multiple Sequence Alignment (MSA)

#### MAFFT (Recommended for Most Cases)

**When to use**: General purpose alignment, medium to large datasets (up to thousands of sequences), balance of speed and accuracy

**Algorithm Selection**:
| Dataset Size | Recommended Algorithm | Command Flag |
|--------------|----------------------|--------------|
| <200 sequences, high accuracy | L-INS-i | `--localpair --maxiterate 1000` |
| <2,000 sequences | G-INS-i | `--globalpair --maxiterate 1000` |
| >2,000 sequences | FFT-NS-2 (default) | `--auto` |
| Ultra-large (>50,000) | PartTree | `--parttree` |

**Implementation Pattern**:
```python
import subprocess
from pathlib import Path
from typing import List, Optional
from Bio import AlignIO

def run_mafft_alignment(
    input_fasta: str,
    output_alignment: str,
    algorithm: str = "auto",
    threads: int = -1,
    max_iterate: Optional[int] = None,
    reorder: bool = True
) -> AlignmentResult:
    """
    Execute MAFFT multiple sequence alignment.

    Args:
        input_fasta: Path to input FASTA file
        output_alignment: Path for output alignment file
        algorithm: Alignment strategy ("auto", "localpair", "globalpair", "parttree")
        threads: Number of CPU threads (-1 = all available)
        max_iterate: Maximum refinement iterations (for iterative modes)
        reorder: Order sequences by similarity in output

    Returns:
        AlignmentResult with statistics and file paths
    """
    # Build command
    cmd = ["mafft"]

    if algorithm == "auto":
        cmd.append("--auto")
    elif algorithm == "localpair":
        cmd.extend(["--localpair", "--maxiterate", str(max_iterate or 1000)])
    elif algorithm == "globalpair":
        cmd.extend(["--globalpair", "--maxiterate", str(max_iterate or 1000)])
    elif algorithm == "parttree":
        cmd.append("--parttree")

    if threads > 0:
        cmd.extend(["--thread", str(threads)])
    elif threads == -1:
        cmd.extend(["--thread", "-1"])

    if reorder:
        cmd.append("--reorder")

    cmd.extend([input_fasta])

    # Execute
    print(f"Running MAFFT: {' '.join(cmd)}")
    with open(output_alignment, "w") as outfile:
        result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"MAFFT failed: {result.stderr}")

    # Parse results
    alignment = AlignIO.read(output_alignment, "fasta")

    return AlignmentResult(
        alignment_file=output_alignment,
        num_sequences=len(alignment),
        alignment_length=alignment.get_alignment_length(),
        average_identity=calculate_average_identity(alignment),
        execution_time=None  # Could add timing
    )

def calculate_average_identity(alignment) -> float:
    """Calculate average pairwise sequence identity."""
    from Bio.Align import PairwiseAligner

    aligner = PairwiseAligner()
    aligner.mode = "global"

    identities = []
    seqs = [str(record.seq) for record in alignment]

    # Sample pairwise comparisons if too many sequences
    import itertools
    if len(seqs) > 100:
        import random
        pairs = random.sample(list(itertools.combinations(range(len(seqs)), 2)), 100)
    else:
        pairs = itertools.combinations(range(len(seqs)), 2)

    for i, j in pairs:
        score = aligner.score(seqs[i], seqs[j])
        max_score = max(len(seqs[i]), len(seqs[j])) * aligner.match_score
        identity = score / max_score
        identities.append(identity)

    return sum(identities) / len(identities) if identities else 0.0
```

#### ClustalOmega (For Very Large Datasets)

**When to use**: >10,000 sequences, profile alignment, HMM-based alignment

```python
def run_clustalo_alignment(
    input_fasta: str,
    output_alignment: str,
    threads: int = 4,
    iterations: int = 3,
    full_iter: bool = True
) -> AlignmentResult:
    """
    Execute ClustalOmega alignment.
    Better for very large datasets than MAFFT.
    """
    cmd = [
        "clustalo",
        "-i", input_fasta,
        "-o", output_alignment,
        "--threads", str(threads),
        "-v",  # Verbose
        "--outfmt", "fasta"
    ]

    if iterations > 0:
        cmd.extend(["--iterations", str(iterations)])

    if full_iter:
        cmd.append("--full-iter")

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"ClustalOmega failed: {result.stderr}")

    alignment = AlignIO.read(output_alignment, "fasta")

    return AlignmentResult(
        alignment_file=output_alignment,
        num_sequences=len(alignment),
        alignment_length=alignment.get_alignment_length()
    )
```

#### MUSCLE (For High Accuracy on Small Datasets)

**When to use**: <500 sequences, when accuracy is critical and speed is secondary

```python
def run_muscle_alignment(
    input_fasta: str,
    output_alignment: str,
    maxiters: int = 16,
    diags: bool = True
) -> AlignmentResult:
    """
    Execute MUSCLE alignment.
    Higher accuracy than MAFFT for small datasets, slower.
    """
    cmd = [
        "muscle",
        "-in", input_fasta,
        "-out", output_alignment,
        "-maxiters", str(maxiters)
    ]

    if diags:
        cmd.append("-diags")  # Faster for diagonally dominant alignments

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"MUSCLE failed: {result.stderr}")

    alignment = AlignIO.read(output_alignment, "fasta")

    return AlignmentResult(
        alignment_file=output_alignment,
        num_sequences=len(alignment),
        alignment_length=alignment.get_alignment_length()
    )
```

### 2. Alignment Processing and Cleaning

#### Trim Poorly Aligned Regions

```python
def trim_alignment(
    input_alignment: str,
    output_alignment: str,
    method: str = "automated1",
    gap_threshold: float = 0.5,
    conservation_threshold: float = 0.0
) -> TrimReport:
    """
    Trim alignment using trimAl to remove unreliable regions.

    Methods:
    - automated1: Automatic selection based on similarity statistics
    - gappyout: Only remove gappy columns (good for phylogenetics)
    - strict: Strict trimming based on gap and similarity scores
    - nogaps: Remove all columns with any gaps

    Args:
        input_alignment: Path to input alignment
        output_alignment: Path for trimmed alignment
        method: Trimming strategy
        gap_threshold: Maximum allowed gap proportion (for manual methods)
        conservation_threshold: Minimum conservation score
    """
    cmd = [
        "trimal",
        "-in", input_alignment,
        "-out", output_alignment,
        "-" + method
    ]

    if method in ["gt", "gat"]:
        cmd.extend(["-gt", str(gap_threshold)])

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"trimAl failed: {result.stderr}")

    # Calculate statistics
    original = AlignIO.read(input_alignment, "fasta")
    trimmed = AlignIO.read(output_alignment, "fasta")

    return TrimReport(
        original_length=original.get_alignment_length(),
        trimmed_length=trimmed.get_alignment_length(),
        columns_removed=original.get_alignment_length() - trimmed.get_alignment_length(),
        percent_remaining=(trimmed.get_alignment_length() / original.get_alignment_length()) * 100,
        method_used=method
    )
```

#### Remove Duplicate Sequences

```python
def deduplicate_sequences(
    input_fasta: str,
    output_fasta: str,
    similarity_threshold: float = 1.0  # 1.0 = identical
) -> DedupReport:
    """
    Remove duplicate or near-duplicate sequences from dataset.
    Important for phylogenetics to avoid bias.
    """
    from Bio import SeqIO
    from collections import defaultdict

    sequences = list(SeqIO.parse(input_fasta, "fasta"))

    # Cluster by sequence identity
    clusters = []
    used = set()

    for i, seq1 in enumerate(sequences):
        if i in used:
            continue

        cluster = [seq1]
        used.add(i)

        for j, seq2 in enumerate(sequences[i+1:], start=i+1):
            if j in used:
                continue

            identity = calculate_sequence_identity(str(seq1.seq), str(seq2.seq))
            if identity >= similarity_threshold:
                cluster.append(seq2)
                used.add(j)

        clusters.append(cluster)

    # Keep one representative per cluster (longest or first)
    representatives = []
    for cluster in clusters:
        rep = max(cluster, key=lambda x: len(x.seq))
        representatives.append(rep)

    # Write output
    SeqIO.write(representatives, output_fasta, "fasta")

    return DedupReport(
        original_count=len(sequences),
        unique_count=len(representatives),
        clusters_formed=len(clusters),
        representatives_written=output_fasta
    )

def calculate_sequence_identity(seq1: str, seq2: str) -> float:
    """Calculate percent identity between two sequences."""
    from Bio.Align import PairwiseAligner

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.gap_score = 0

    score = aligner.score(seq1, seq2)
    max_len = max(len(seq1), len(seq2))

    return score / max_len
```

### 3. Variant Calling from Alignment

#### Identify Mutations

```python
from dataclasses import dataclass
from typing import List, Dict, Tuple
from collections import Counter
import numpy as np

@dataclass
class Variant:
    position: int  # 1-based position in alignment
    reference_base: str
    alternate_base: str
    reference_seq_id: str
    frequency: float
    affected_sequences: List[str]
    variant_type: str  # "SNP", "insertion", "deletion", "mixed"

def call_variants_from_alignment(
    alignment_file: str,
    reference_id: str,
    min_frequency: float = 0.01,
    ignore_gaps: bool = True
) -> List[Variant]:
    """
    Identify variants from multiple sequence alignment.

    Args:
        alignment_file: Path to aligned FASTA
        reference_id: ID of reference sequence for coordinate system
        min_frequency: Minimum frequency to report variant (0.01 = 1%)
        ignore_gaps: Whether to ignore gap characters in analysis

    Returns:
        List of Variant objects
    """
    alignment = AlignIO.read(alignment_file, "fasta")

    # Find reference sequence
    reference_idx = None
    for i, record in enumerate(alignment):
        if record.id == reference_id or reference_id in record.id:
            reference_idx = i
            break

    if reference_idx is None:
        raise ValueError(f"Reference {reference_id} not found in alignment")

    reference_seq = str(alignment[reference_idx].seq)
    num_sequences = len(alignment)

    variants = []

    # Analyze each column
    for col_idx in range(alignment.get_alignment_length()):
        column = [str(record.seq[col_idx]) for record in alignment]
        ref_base = reference_seq[col_idx]

        # Skip if reference is gap (unless tracking indels)
        if ref_base == "-":
            continue

        # Count frequencies
        counts = Counter(column)
        total = sum(counts.values())

        # Skip conserved positions
        if counts[ref_base] == total:
            continue

        # Identify alternate alleles
        for base, count in counts.items():
            if base == ref_base:
                continue
            if ignore_gaps and base == "-":
                continue

            frequency = count / total
            if frequency < min_frequency:
                continue

            # Determine variant type
            if ref_base == "-":
                vtype = "insertion"
            elif base == "-":
                vtype = "deletion"
            else:
                vtype = "SNP"

            # Get sequences with this variant
            affected = [
                alignment[i].id for i, b in enumerate(column) 
                if b == base
            ]

            variant = Variant(
                position=col_idx + 1,  # Convert to 1-based
                reference_base=ref_base,
                alternate_base=base,
                reference_seq_id=reference_id,
                frequency=frequency,
                affected_sequences=affected,
                variant_type=vtype
            )
            variants.append(variant)

    return sorted(variants, key=lambda x: x.position)
```

#### Amino Acid Translation and Impact

```python
def translate_variants_to_protein(
    nucleotide_variants: List[Variant],
    reference_dna: str,
    reading_frame: int = 0
) -> List[ProteinVariant]:
    """
    Translate DNA variants to protein impact.

    Returns synonymous, missense, and nonsense mutations.
    """
    from Bio.Seq import Seq

    protein_variants = []

    for variant in nucleotide_variants:
        if variant.variant_type != "SNP":
            continue  # Skip indels for simple translation

        # Calculate codon position
        codon_start = ((variant.position - 1) // 3) * 3
        codon_pos = (variant.position - 1) % 3

        # Get reference and mutant codons
        ref_codon = reference_dna[codon_start:codon_start+3]
        mut_codon = list(ref_codon)
        mut_codon[codon_pos] = variant.alternate_base
        mut_codon = "".join(mut_codon)

        # Translate
        ref_aa = str(Seq(ref_codon).translate())
        mut_aa = str(Seq(mut_codon).translate())

        # Determine impact
        if ref_aa == mut_aa:
            impact = "synonymous"
        elif mut_aa == "*":
            impact = "nonsense"
        elif ref_aa == "*":
            impact = "readthrough"
        else:
            impact = "missense"

        protein_variants.append(ProteinVariant(
            dna_position=variant.position,
            codon_number=codon_start // 3 + 1,
            reference_aa=ref_aa,
            alternate_aa=mut_aa,
            impact=impact,
            frequency=variant.frequency
        ))

    return protein_variants
```

### 4. Phylogenetic Tree Construction

#### IQ-TREE (Recommended)

```python
def build_phylogeny_iqtree(
    alignment_file: str,
    output_prefix: str,
    model: str = "MFP",  # ModelFinder Plus
    bootstrap_replicates: int = 1000,
    threads: int = 4,
    additional_options: List[str] = None
) -> PhylogenyResult:
    """
    Construct phylogenetic tree using IQ-TREE.

    Features:
    - Automatic model selection (MFP)
    - Ultrafast bootstrap (UFBoot)
    - SH-aLRT branch support
    - Tree topology tests

    Args:
        alignment_file: Path to trimmed alignment
        output_prefix: Prefix for output files
        model: Substitution model ("MFP" for auto-selection)
        bootstrap_replicates: Number of bootstrap replicates (1000 recommended)
        threads: CPU threads to use
        additional_options: Extra IQ-TREE flags
    """
    cmd = [
        "iqtree2",  # or "iqtree" depending on version
        "-s", alignment_file,
        "-pre", output_prefix,
        "-m", model,
        "-bb", str(bootstrap_replicates),  # Ultrafast bootstrap
        "-alrt", "1000",  # SH-aLRT test
        "-nt", str(threads),
        "-quiet"  # Reduce verbosity
    ]

    if additional_options:
        cmd.extend(additional_options)

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"IQ-TREE failed: {result.stderr}")

    # Parse results
    tree_file = f"{output_prefix}.treefile"
    log_file = f"{output_prefix}.log"

    # Extract best model from log
    best_model = parse_iqtree_model(log_file)

    # Load tree
    from dendropy import Tree
    tree = Tree.get(path=tree_file, schema="newick")

    return PhylogenyResult(
        tree_file=tree_file,
        best_model=best_model,
        log_likelihood=parse_iqtree_likelihood(log_file),
        num_taxa=len(tree.taxon_namespace),
        tree_object=tree
    )

def parse_iqtree_model(log_file: str) -> str:
    """Extract best-fit model from IQ-TREE log."""
    with open(log_file, "r") as f:
        for line in f:
            if "Best-fit model" in line:
                return line.split(":")[-1].strip()
    return "Unknown"
```

#### Time-Resolved Phylogeny (for Viruses)

```python
def build_timetree(
    alignment_file: str,
    metadata_file: str,  # TSV with strain, date columns
    output_prefix: str,
    clock_rate: float = None  # Optional clock rate prior
) -> TimetreeResult:
    """
    Build time-resolved phylogeny using TreeTime or LSD.
    Essential for viral phylogeography and dynamics.
    """
    # First build standard tree
    tree_result = build_phylogeny_iqtree(alignment_file, output_prefix)

    # Run TreeTime
    cmd = [
        "treetime",
        "--tree", tree_result.tree_file,
        "--aln", alignment_file,
        "--dates", metadata_file,
        "--outdir", f"{output_prefix}_treetime"
    ]

    if clock_rate:
        cmd.extend(["--clock-rate", str(clock_rate)])
    else:
        cmd.append("--clock-filter")  # Infer clock rate

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"TreeTime failed: {result.stderr}")

    return TimetreeResult(
        tree_file=f"{output_prefix}_treetime/timetree.nexus",
        divergence_times=parse_treetime_dates(f"{output_prefix}_treetime/dates.tsv"),
        clock_rate=parse_treetime_clock(f"{output_prefix}_treetime/clock_rate.txt")
    )
```

### 5. Conservation and Diversity Analysis

#### Calculate Conservation Scores

```python
def calculate_conservation(
    alignment_file: str,
    method: str = "shannon"
) -> List[ConservationScore]:
    """
    Calculate per-column conservation scores.

    Methods:
    - shannon: Shannon entropy (0 = conserved, high = variable)
    - simpson: Simpson diversity index
    - gini: Gini coefficient
    """
    from scipy.stats import entropy

    alignment = AlignIO.read(alignment_file, "fasta")
    scores = []

    for col_idx in range(alignment.get_alignment_length()):
        column = [str(record.seq[col_idx]) for record in alignment]

        # Remove gaps for conservation calculation
        residues = [c for c in column if c != "-"]
        if not residues:
            continue

        counts = Counter(residues)
        total = sum(counts.values())
        freqs = [c / total for c in counts.values()]

        if method == "shannon":
            score = entropy(freqs, base=2)
        elif method == "simpson":
            score = 1 - sum(f**2 for f in freqs)
        else:
            raise ValueError(f"Unknown method: {method}")

        # Normalize to 0-1 (1 = fully conserved)
        max_entropy = np.log2(len(set(residues))) if len(set(residues)) > 1 else 1
        normalized_score = 1 - (score / max_entropy) if max_entropy > 0 else 1

        scores.append(ConservationScore(
            position=col_idx + 1,
            score=normalized_score,
            dominant_residue=counts.most_common(1)[0][0],
            dominant_frequency=counts.most_common(1)[0][1] / total
        ))

    return scores
```

#### Detect Recombination

```python
def detect_recombination(
    alignment_file: str,
    method: str = "gubbins"  # or "3seq", "gard"
) -> RecombinationReport:
    """
    Detect recombination events in viral alignment.
    Critical for RNA viruses and segmented genomes.
    """
    if method == "gubbins":
        cmd = [
            "run_gubbins.py",
            alignment_file,
            "--threads", "4",
            "--verbose"
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)

        # Parse Gubbins output
        return parse_gubbins_output(alignment_file)

    elif method == "3seq":
        # Use HyPhy 3seq
        pass

    else:
        raise ValueError(f"Unknown method: {method}")
```

---

## Integration Workflows

### Complete Viral Variant Analysis Pipeline

```python
async def analyze_viral_variants(
    virus_name: str,
    gene: str,
    sample_ids: List[str],
    reference_id: str,
    output_dir: str
) -> VariantAnalysisReport:
    """
    End-to-end viral variant analysis workflow.

    Steps:
    1. Fetch sequences from database
    2. Align sequences
    3. Trim alignment
    4. Call variants
    5. Build phylogeny
    6. Annotate functional impact
    7. Generate report
    """
    from pathlib import Path
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Step 1: Data retrieval
    print("Fetching sequences...")
    fasta_file = f"{output_dir}/sequences.fasta"
    await fetch_sequences_from_database(virus_name, gene, sample_ids, fasta_file)

    # Step 2: Deduplication
    print("Removing duplicates...")
    dedup_file = f"{output_dir}/sequences_dedup.fasta"
    dedup_report = deduplicate_sequences(fasta_file, dedup_file)
    print(f"Removed {dedup_report.original_count - dedup_report.unique_count} duplicates")

    # Step 3: Alignment
    print("Aligning sequences...")
    alignment_file = f"{output_dir}/alignment.fasta"
    if dedup_report.unique_count > 1000:
        align_result = run_clustalo_alignment(dedup_file, alignment_file)
    else:
        align_result = run_mafft_alignment(dedup_file, alignment_file, algorithm="auto")

    # Step 4: Trimming
    print("Trimming alignment...")
    trimmed_alignment = f"{output_dir}/alignment_trimmed.fasta"
    trim_report = trim_alignment(alignment_file, trimmed_alignment, method="automated1")
    print(f"Trimmed {trim_report.columns_removed} columns ({100-trim_report.percent_remaining:.1f}% removed)")

    # Step 5: Variant calling
    print("Calling variants...")
    variants = call_variants_from_alignment(
        trimmed_alignment,
        reference_id,
        min_frequency=0.01
    )

    # Step 6: Phylogeny
    print("Building phylogeny...")
    tree_prefix = f"{output_dir}/phylogeny"
    tree_result = build_phylogeny_iqtree(
        trimmed_alignment,
        tree_prefix,
        bootstrap_replicates=1000
    )

    # Step 7: Conservation analysis
    print("Analyzing conservation...")
    conservation = calculate_conservation(trimmed_alignment)

    # Step 8: Generate report
    report = VariantAnalysisReport(
        input_sequences=dedup_report,
        alignment_stats=align_result,
        trimming_stats=trim_report,
        variants_found=len(variants),
        variant_list=variants,
        phylogeny=tree_result,
        conservation_scores=conservation,
        output_files={
            "alignment": trimmed_alignment,
            "tree": tree_result.tree_file,
            "variants": f"{output_dir}/variants.csv"
        }
    )

    # Save variants to CSV
    save_variants_to_csv(variants, f"{output_dir}/variants.csv")

    return report

def save_variants_to_csv(variants: List[Variant], output_file: str):
    """Export variants to CSV format."""
    import csv

    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Position", "Reference", "Alternate", "Type", 
            "Frequency", "Num_Sequences"
        ])
        for v in variants:
            writer.writerow([
                v.position, v.reference_base, v.alternate_base,
                v.variant_type, f"{v.frequency:.4f}", len(v.affected_sequences)
            ])
```

---

## Best Practices

### For Viral Genomics

1. **Reference Selection**: Use a well-characterized reference strain (e.g., Wuhan-Hu-1 for SARS-CoV-2, Makona for Ebola)

2. **Sampling Strategy**: 
   - Avoid over-sampling from single outbreaks (reduces diversity)
   - Include geographic and temporal diversity
   - Balance sample sizes across regions

3. **Alignment Parameters**:
   - Use `--auto` in MAFFT for most viral datasets
   - For highly divergent viruses (e.g., HIV), use `--localpair` with high iterations
   - For whole genomes, consider aligning genes separately if recombination suspected

4. **Variant Filtering**:
   - Remove variants in primer binding sites (ARTIFACTS)
   - Filter by quality score if available
   - Require minimum read depth for NGS data

5. **Phylogenetic Considerations**:
   - Always use trimmed alignments for tree building
   - Check for recombination before tree building
   - Use time-resolved trees for rapidly evolving viruses
   - Report bootstrap/SH-aLRT support values

### Quality Control Checklist

| Step | Check | Action if Failed |
|------|-------|------------------|
| Alignment | Average pairwise identity <50% | Split into subtypes/genes |
| Alignment | >50% gaps in reference | Re-align with different method |
| Variants | Stop codons in coding regions | Check reading frame |
| Variants | >5% Ns in variant positions | Filter low-quality calls |
| Phylogeny | Long branch attraction suspected | Use outgroup, check model |
| Phylogeny | Bootstrap <70% for key nodes | Add more data or sequences |

### Common Pitfalls

| Pitfall | Solution |
|---------|----------|
| Including too many identical sequences | Deduplicate at 99.9% identity |
| Aligning CDS with UTRs | Separate UTRs or use codon-aware aligner |
| Ignoring recombination | Run GARD or 3seq before tree building |
| Using DNA models for protein-coding | Use codon models or translate first |
| Bootstrapping with thousands of sequences | Use ultrafast bootstrap (UFBoot) |

---

## Troubleshooting

### Alignment Issues

**Problem**: MAFFT fails with "memory error" on large dataset
**Solutions**:
- Use `--parttree` algorithm for >10,000 sequences
- Split dataset into smaller batches
- Use ClustalOmega instead (more memory efficient)
- Increase system swap space

**Problem**: Poor alignment of highly variable regions
**Solutions**:
- Use `--localpair` for local alignment
- Trim problematic regions with trimAl
- Align conserved domains separately
- Use profile alignment (align to existing MSA)

### Variant Calling Issues

**Problem**: Too many variants called (false positives)
**Solutions**:
- Increase `min_frequency` threshold
- Filter by quality scores from sequencing
- Remove variants in homopolymer regions
- Require variants to be present in both directions (for NGS)

**Problem**: Missing known variants
**Solutions**:
- Check reference sequence matches expected
- Verify coordinates (0-based vs 1-based)
- Lower `min_frequency` threshold
- Check for strand bias

### Phylogeny Issues

**Problem**: IQ-TREE fails to converge
**Solutions**:
- Reduce bootstrap replicates (100 instead of 1000)
- Use simpler model (GTR instead of GTR+G+I)
- Check for identical sequences (remove duplicates)
- Use FastTree for quick approximate tree

**Problem**: Tree shows unexpected topology
**Solutions**:
- Check for recombination (use SplitsTree4)
- Verify outgroup is appropriate
- Check for long branch attraction
- Try different rooting methods

---

## Examples

### Example 1: SARS-CoV-2 Spike Protein Analysis

```python
# Analyze recent SARS-CoV-2 spike sequences
report = await analyze_viral_variants(
    virus_name="SARS-CoV-2",
    gene="S",
    sample_ids=["latest_1000_global"],  # Fetch from GISAID
    reference_id="Wuhan-Hu-1",
    output_dir="./sars2_spike_analysis"
)

# Expected outputs:
# - alignment_trimmed.fasta: Clean alignment
# - phylogeny.treefile: Time-resolved tree
# - variants.csv: All mutations with frequencies
# - Conservation scores highlighting key antigenic sites

# Key findings to look for:
# - Mutations in receptor binding domain (RBD: 319-541)
# - Changes in furin cleavage site
# - N-linked glycosylation site modifications
```

### Example 2: Influenza H3N2 Hemagglutinin

```python
# Seasonal influenza analysis
# Important: Account for egg-adapted changes if applicable

report = await analyze_viral_variants(
    virus_name="Influenza-A",
    gene="HA",
    sample_ids=["H3N2_2023_2024"],
    reference_id="A_Wisconsin_67_2005",
    output_dir="./flu_h3n2_analysis"
)

# Special considerations:
# - Check for egg-adapted mutations (T160K, Q226R)
# - Focus on antigenic sites (Sa, Sb, Ca1, Ca2, Cb)
# - Look for glycosylation site additions (immune evasion)
```

### Example 3: Ebola Virus GP Analysis

```python
# Ebola glycoprotein with recombination check
alignment = run_mafft_alignment("ebola_gp.fasta", "ebola_aligned.fasta")

# Check for recombination (important for filoviruses)
recomb = detect_recombination("ebola_aligned.fasta", method="gubbins")

if recomb.events_found > 0:
    print(f"Warning: {recomb.events_found} recombination events detected")
    # Use recombination-aware phylogeny methods

# Build tree
 tree = build_phylogeny_iqtree("ebola_aligned.fasta", "ebola_tree")
```

---

## References

1. **MAFFT**: Katoh & Standley, Nucleic Acids Res 2013 (doi:10.1093/nar/gkt010)
2. **ClustalOmega**: Sievers et al., Mol Syst Biol 2011 (doi:10.1038/msb.2011.75)
3. **IQ-TREE**: Minh et al., Mol Biol Evol 2020 (doi:10.1093/molbev/msaa015)
4. **trimAl**: Capella-Gutierrez et al., Bioinformatics 2009 (doi:10.1093/bioinformatics/btp348)
5. **Gubbins**: Croucher et al., Nucleic Acids Res 2015 (doi:10.1093/nar/gku1196)
6. **Nextstrain/Augur**: Hadfield et al., Bioinformatics 2018 (doi:10.1093/bioinformatics/bty407)

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2024-01-15 | Initial release with MAFFT, IQ-TREE support |
| 1.1 | 2024-04-10 | Added ClustalOmega for large datasets |
| 1.2 | 2024-07-20 | Integrated recombination detection |
| 1.3 | 2024-10-05 | Added time-resolved phylogeny methods |

---

## Related Skills

- **protein-structure-prediction.md**: Structure analysis of variants
- **molecular-docking.md**: Drug interaction analysis
- **phylogenetics.md**: Advanced evolutionary analysis
- **database-navigation.md**: Fetching sequences from NCBI/GISAID
- **mutation-impact-prediction.md**: Functional impact of variants
