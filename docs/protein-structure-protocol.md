# Skill: Protein Structure Prediction & Analysis

## Overview

This skill enables the AI agent to predict protein 3D structures from amino acid sequences and analyze the structural impact of mutations. It covers structure prediction methods, quality assessment, mutation impact analysis, binding site identification, and structure comparison techniques essential for viral protein research.

---

## Capabilities

- Predict protein structures using state-of-the-art methods (AlphaFold, ESMFold)
- Assess prediction quality and confidence metrics
- Analyze how mutations affect protein stability and function
- Identify druggable binding sites and functional domains
- Compare structures to detect conformational changes
- Integrate structural data with sequence-based analyses

---

## Prerequisites

### Required Software

| Tool | Installation | Purpose |
|------|--------------|---------|
| **ColabFold** | `pip install colabfold[alphafold]` | High-accuracy structure prediction |
| **ESMFold** | `pip install fair-esm` | Fast structure prediction |
| **FoldX** | Download from http://foldxsuite.crg.eu | Stability change prediction |
| **FPocket** | `conda install -c conda-forge fpocket` | Binding site detection |
| **Biopython** | `pip install biopython` | Structure manipulation |
| **MDAnalysis** | `pip install mdanalysis` | Structure analysis |
| **PyMOL/NGView** | `conda install -c conda-forge nglview` | Visualization |

### Required Data Access

- **AlphaFold Database**: https://alphafold.ebi.ac.uk/ (direct download, no API key)
- **PDB**: https://www.rcsb.org/ (REST API, no authentication required)
- **UniProt**: For protein sequence retrieval

### Hardware Requirements

- **Minimum**: CPU-only mode (slow, ~30 min per protein)
- **Recommended**: NVIDIA GPU with 8GB+ VRAM (CUDA-enabled)
- **Optimal**: Google Colab Pro or local RTX 3090/A100 for batch processing

---

## Core Methods

### 1. Structure Prediction

#### AlphaFold/ColabFold (Recommended)

**When to use**: High-accuracy prediction of single proteins or complexes; viral proteins up to ~2000 residues

**Input**: 
- Amino acid sequence (FASTA format)
- Optional: Multiple sequence alignment (MSA) for homology-rich proteins

**Output**:
- PDB format structure file
- JSON with pLDDT confidence scores per residue
- PAE (Predicted Aligned Error) matrix for domain packing confidence

**Confidence Interpretation**:
| pLDDT Range | Confidence Level | Interpretation |
|-------------|------------------|----------------|
| >90 | Very high | Structure suitable for drug design |
| 70-90 | High | Reliable backbone, side chains may vary |
| 50-70 | Low | Correct fold, uncertain details |
| <50 | Very low | Unstructured or disordered region |

**Implementation Pattern**:
```python
from colabfold.batch import get_queries, run

def predict_structure_alphafold(
    sequence: str, 
    job_name: str,
    use_templates: bool = True,
    num_models: int = 5,
    rank_by: str = "pLDDT"
) -> StructurePrediction:
    """
    Predict protein structure using ColabFold/AlphaFold.

    Args:
        sequence: Amino acid sequence (single letter code)
        job_name: Identifier for this prediction run
        use_templates: Whether to use PDB templates (recommended for viral proteins)
        num_models: Number of models to generate (1-5)
        rank_by: Metric to select best model ("pLDDT" or "pTMscore")

    Returns:
        StructurePrediction object with PDB path, scores, and metadata
    """
    # Write FASTA
    fasta_path = f"{job_name}.fasta"
    with open(fasta_path, "w") as f:
        f.write(f">{job_name}\n{sequence}\n")

    # Run prediction
    queries = get_queries(fasta_path)
    results = run(
        queries=queries,
        result_dir=f"./{job_name}_results",
        use_templates=use_templates,
        num_models=num_models,
        rank_mode=rank_by,
        model_type="auto"  # Will select AlphaFold2-multimer for complexes
    )

    return StructurePrediction(
        pdb_path=f"./{job_name}_results/{job_name}_unrelaxed_rank1.pdb",
        plddt_scores=results["plddt"],
        ptm_score=results["ptm"],
        job_name=job_name
    )
```

#### ESMFold (Fast Alternative)

**When to use**: Quick screening of many sequences, moderate accuracy acceptable, limited compute

**Input**: Amino acid sequence
**Output**: PDB format + per-residue confidence

**Implementation**:
```python
import torch
from esm import pretrained

def predict_structure_esmfold(sequence: str) -> StructurePrediction:
    """Fast structure prediction using ESMFold."""
    model = pretrained.esmfold_v1()
    model = model.eval().cuda()

    with torch.no_grad():
        output = model.infer_pdb(sequence)

    return StructurePrediction(
        pdb_content=output,
        method="ESMFold",
        confidence=model.predict_confidence(sequence)
    )
```

### 2. Quality Assessment

#### pLDDT Analysis

```python
def assess_structure_quality(pdb_path: str) -> QualityReport:
    """
    Analyze AlphaFold pLDDT scores and structural quality.

    Returns:
        QualityReport with overall score, problematic regions, and recommendations
    """
    from Bio.PDB import PDBParser
    import numpy as np

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)

    # Extract B-factors (AlphaFold stores pLDDT here)
    plddt_scores = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.has_id("CA"):
                    plddt_scores.append(residue["CA"].get_bfactor())

    plddt_array = np.array(plddt_scores)

    # Identify low-confidence regions
    low_confidence = np.where(plddt_array < 70)[0]
    very_low = np.where(plddt_array < 50)[0]

    return QualityReport(
        mean_plddt=float(np.mean(plddt_array)),
        min_plddt=float(np.min(plddt_array)),
        low_confidence_regions=low_confidence.tolist(),
        very_low_confidence_regions=very_low.tolist(),
        assessment="High quality" if np.mean(plddt_array) > 80 else "Use with caution"
    )
```

#### Structural Validation

```python
def validate_structure_geometry(pdb_path: str) -> ValidationReport:
    """
    Check for steric clashes, Ramachandran outliers, and bond geometry.
    Uses MolProbity criteria adapted for predicted structures.
    """
    # Implementation using Bio.PDB or external call to MolProbity
    pass
```

### 3. Mutation Impact Analysis

#### Stability Change Prediction (ΔΔG)

**FoldX Implementation**:
```python
import subprocess
import os

def calculate_stability_change(
    wild_type_pdb: str,
    mutation_position: int,
    wild_type_aa: str,
    mutant_aa: str,
    chain_id: str = "A"
) -> StabilityPrediction:
    """
    Predict stability change using FoldX.

    Args:
        wild_type_pdb: Path to wild-type structure
        mutation_position: Residue number (PDB numbering)
        wild_type_aa: Single letter wild-type amino acid
        mutant_aa: Single letter mutant amino acid
        chain_id: Chain identifier

    Returns:
        StabilityPrediction with ΔΔG and interpretation
    """
    # Prepare structure (repair PDB)
    repair_cmd = f"foldx --command=RepairPDB --pdb={wild_type_pdb}"
    subprocess.run(repair_cmd, shell=True, check=True)

    repaired_pdb = wild_type_pdb.replace(".pdb", "_Repair.pdb")

    # Create mutation file
    mutation_string = f"{wild_type_aa}{chain_id}{mutation_position}{mutant_aa};"
    with open("individual_list.txt", "w") as f:
        f.write(mutation_string)

    # Run FoldX BuildModel
    build_cmd = (
        f"foldx --command=BuildModel --pdb={repaired_pdb} "
        f"--mutant-file=individual_list.txt"
    )
    result = subprocess.run(build_cmd, shell=True, capture_output=True, text=True)

    # Parse ΔΔG from output
    ddg = parse_foldx_ddg("BuildModel_FoldXOutput.txt")

    interpretation = interpret_ddg(ddg)

    return StabilityPrediction(
        ddG=ddg,
        interpretation=interpretation,
        confidence="High" if abs(ddg) > 2 else "Low"
    )

def interpret_ddg(ddg: float) -> str:
    """Interpret FoldX ΔΔG values."""
    if ddg > 5:
        return "Strongly destabilizing - likely disrupts folding"
    elif ddg > 2:
        return "Destabilizing - may affect stability"
    elif ddg > -2:
        return "Neutral - minimal stability impact"
    elif ddg > -5:
        return "Stabilizing - may enhance stability"
    else:
        return "Strongly stabilizing - rare"
```

**Alternative: Rosetta ddG**:
```python
def calculate_rosetta_ddg(
    wild_type_pdb: str,
    mutant_sequence: str
) -> float:
    """
    Use Rosetta cartesian_ddg protocol for high-accuracy stability prediction.
    More accurate than FoldX but computationally expensive.
    """
    # Requires Rosetta installation
    pass
```

#### Structural Environment Analysis

```python
def analyze_mutation_environment(
    structure_path: str,
    position: int,
    radius: float = 5.0
) -> EnvironmentAnalysis:
    """
    Analyze local environment of mutation site.

    Returns:
        - Neighboring residues within radius
        - Secondary structure context
        - Solvent accessibility
        - Known functional annotations
    """
    from Bio.PDB import PDBParser, DSSP
    from Bio.PDB.SASA import ShrakeRupley

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", structure_path)
    model = structure[0]

    # Get target residue
    chain = model["A"]
    target_residue = chain[position]

    # Find neighbors
    neighbors = []
    target_coord = target_residue["CA"].get_coord()

    for residue in chain:
        if residue != target_residue and residue.has_id("CA"):
            dist = np.linalg.norm(residue["CA"].get_coord() - target_coord)
            if dist <= radius:
                neighbors.append({
                    "residue": residue.get_resname(),
                    "position": residue.get_id()[1],
                    "distance": float(dist)
                })

    # Calculate SASA
    sr = ShrakeRupley()
    sr.compute(structure, level="R")
    sasa = target_residue.sasa

    # Secondary structure (requires DSSP)
    # dssp = DSSP(model, structure_path)
    # ss = dssp[(chain.id, target_residue.id)]

    return EnvironmentAnalysis(
        neighbors=neighbors,
        solvent_accessibility=sasa,
        num_neighbors=len(neighbors),
        is_buried=sasa < 20.0,
        is_surface=sasa > 100.0
    )
```

### 4. Binding Site Identification

#### Geometric Detection (FPocket)

```python
def identify_binding_sites(pdb_path: str) -> List[BindingSite]:
    """
    Detect potential binding pockets using FPocket.

    Returns:
        List of binding sites with scores and druggability estimates
    """
    import subprocess
    import json

    # Run FPocket
    cmd = f"fpocket -f {pdb_path}"
    subprocess.run(cmd, shell=True, check=True)

    # Parse output
    output_dir = pdb_path.replace(".pdb", "_out")
    pockets = parse_fpocket_output(output_dir)

    return sorted(pockets, key=lambda x: x.score, reverse=True)

def parse_fpocket_output(output_dir: str) -> List[BindingSite]:
    """Parse FPocket output files."""
    pockets = []
    info_file = os.path.join(output_dir, "info.txt")

    with open(info_file, "r") as f:
        for line in f:
            if line.startswith("Pocket"):
                # Parse pocket details
                pass

    return pockets
```

#### Energy-Based Detection

```python
def predict_binding_sites_energetic(
    structure_path: str,
    probe_radius: float = 1.4
) -> List[BindingSite]:
    """
    Use SiteMap or similar energy-based methods.
    Identifies sites based on probe binding energy distributions.
    """
    # Implementation using OpenEye SiteMap or similar
    pass
```

### 5. Structure Comparison

#### TM-align for Global Comparison

```python
def compare_structures_tmalign(
    structure1_path: str,
    structure2_path: str
) -> StructuralAlignment:
    """
    Compare two structures using TM-align.

    Returns:
        TM-score (normalized), RMSD, and alignment mapping
    """
    import subprocess
    import re

    cmd = f"TMalign {structure1_path} {structure2_path} -o TM.sup"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    output = result.stdout

    # Parse TM-score
    tm_match = re.search(r"TM-score= ([0-9.]+)", output)
    tm_score = float(tm_match.group(1)) if tm_match else 0.0

    # Parse RMSD
    rmsd_match = re.search(r"RMSD=\s+([0-9.]+)", output)
    rmsd = float(rmsd_match.group(1)) if rmsd_match else 999.0

    interpretation = "Same fold" if tm_score > 0.5 else "Different fold"
    if tm_score > 0.8:
        interpretation = "Highly similar structures"

    return StructuralAlignment(
        tm_score=tm_score,
        rmsd=rmsd,
        alignment_file="TM.sup",
        interpretation=interpretation
    )
```

#### Local RMSD for Mutation Impact

```python
def calculate_local_rmsd(
    wild_type_pdb: str,
    mutant_pdb: str,
    radius: float = 10.0,
    center_position: int = None
) -> float:
    """
    Calculate RMSD in local region around mutation site.
    More sensitive than global RMSD for detecting local conformational changes.
    """
    from Bio.PDB import PDBParser, Superimposer

    parser = PDBParser(QUIET=True)
    wt_structure = parser.get_structure("wt", wild_type_pdb)
    mut_structure = parser.get_structure("mut", mutant_pdb)

    # Extract atoms within radius of mutation site
    wt_atoms = []
    mut_atoms = []

    if center_position:
        wt_center = wt_structure[0]["A"][center_position]["CA"].get_coord()
        mut_center = mut_structure[0]["A"][center_position]["CA"].get_coord()

        for residue in wt_structure[0]["A"]:
            if residue.has_id("CA"):
                dist = np.linalg.norm(residue["CA"].get_coord() - wt_center)
                if dist <= radius:
                    wt_atoms.append(residue["CA"])

        for residue in mut_structure[0]["A"]:
            if residue.has_id("CA"):
                dist = np.linalg.norm(residue["CA"].get_coord() - mut_center)
                if dist <= radius:
                    mut_atoms.append(residue["CA"])

    # Superimpose and calculate RMSD
    sup = Superimposer()
    sup.set_atoms(wt_atoms, mut_atoms)
    sup.apply(mut_structure[0].get_atoms())

    return sup.rms
```

---

## Integration Workflows

### Complete Mutation Analysis Pipeline

```python
async def analyze_mutation_structural_impact(
    wild_type_sequence: str,
    mutation_position: int,
    wild_type_aa: str,
    mutant_aa: str,
    protein_name: str = "protein"
) -> StructuralMutationReport:
    """
    Comprehensive structural analysis of a point mutation.

    Workflow:
    1. Predict wild-type structure
    2. Predict mutant structure  
    3. Calculate stability change (FoldX)
    4. Analyze local environment
    5. Compare structures (RMSD)
    6. Assess functional impact

    Returns:
        Complete structural analysis report
    """
    # Step 1: Predict wild-type structure
    print("Predicting wild-type structure...")
    wt_pred = predict_structure_alphafold(
        wild_type_sequence, 
        f"{protein_name}_WT"
    )

    # Step 2: Create mutant sequence and predict
    mutant_sequence = (
        wild_type_sequence[:mutation_position-1] + 
        mutant_aa + 
        wild_type_sequence[mutation_position:]
    )

    print("Predicting mutant structure...")
    mut_pred = predict_structure_alphafold(
        mutant_sequence,
        f"{protein_name}_{wild_type_aa}{mutation_position}{mutant_aa}"
    )

    # Step 3: Stability analysis
    print("Calculating stability change...")
    stability = calculate_stability_change(
        wt_pred.pdb_path,
        mutation_position,
        wild_type_aa,
        mutant_aa
    )

    # Step 4: Environment analysis
    print("Analyzing local environment...")
    environment = analyze_mutation_environment(
        wt_pred.pdb_path,
        mutation_position
    )

    # Step 5: Structural comparison
    print("Comparing structures...")
    local_rmsd = calculate_local_rmsd(
        wt_pred.pdb_path,
        mut_pred.pdb_path,
        center_position=mutation_position
    )

    global_alignment = compare_structures_tmalign(
        wt_pred.pdb_path,
        mut_pred.pdb_path
    )

    # Step 6: Functional assessment
    functional_impact = assess_functional_impact(
        environment,
        stability,
        local_rmsd
    )

    return StructuralMutationReport(
        wild_type=wt_pred,
        mutant=mut_pred,
        stability=stability,
        environment=environment,
        structural_changes={
            "local_rmsd": local_rmsd,
            "global_tm_score": global_alignment.tm_score,
            "global_rmsd": global_alignment.rmsd
        },
        functional_impact=functional_impact,
        recommendations=generate_recommendations(stability, environment)
    )

def assess_functional_impact(
    environment: EnvironmentAnalysis,
    stability: StabilityPrediction,
    local_rmsd: float
) -> str:
    """Assess likely functional impact based on structural analysis."""
    impacts = []

    if stability.ddG > 2:
        impacts.append("Destabilizing - may reduce protein levels")
    if stability.ddG < -2:
        impacts.append("Stabilizing - may increase protein half-life")

    if environment.is_buried and abs(stability.ddG) > 1:
        impacts.append("Core mutation - high impact on folding")
    if environment.is_surface:
        impacts.append("Surface mutation - may affect interactions")

    if local_rmsd > 2.0:
        impacts.append("Significant conformational change detected")
    elif local_rmsd > 1.0:
        impacts.append("Moderate local structural rearrangement")
    else:
        impacts.append("Minimal structural perturbation")

    if len(environment.neighbors) < 3:
        impacts.append("Flexible loop region - tolerant to change")
    elif len(environment.neighbors) > 8:
        impacts.append("Packed environment - sensitive to change")

    return "; ".join(impacts)
```

---

## Best Practices

### For Viral Proteins

1. **Use Templates When Available**: For SARS-CoV-2 spike, use experimental PDB structures as templates to improve accuracy in variable regions

2. **Account for Glycosylation**: Viral surface proteins are heavily glycosylated. Consider glycan shielding when analyzing antibody epitopes

3. **Domain Boundaries**: Viral proteins often have disordered linkers between domains. Analyze domains separately if pLDDT drops in linker regions

4. **Multimeric State**: Many viral proteins function as trimers (spike) or dimers. Use AlphaFold-Multimer for oligomeric predictions

### Quality Control

1. **Always Check pLDDT**: Never use residues with pLDDT < 50 for detailed analysis
2. **Validate with Experimental Data**: When available, compare predicted structures to cryo-EM/X-ray structures
3. **Cross-Validate Predictions**: Use both FoldX and Rosetta for critical stability predictions
4. **Consider Ensemble Effects**: Single structure may not capture dynamics; use MD simulations for high-stakes predictions

### Common Pitfalls

| Pitfall | Solution |
|---------|----------|
| Ignoring protonation states | Use H++ or PDB2PQR to assign correct protonation at pH 7.4 |
| Using apo structure for docking | Include bound ligands or ions if known to affect conformation |
| Neglecting crystal contacts | Filter out artifacts from crystal packing in experimental structures |
| Over-interpreting low-confidence regions | Restrict drug design to pLDDT > 90 regions |
| Ignoring entropy effects | FoldX/Rosetta estimate enthalpy; dynamics needed for entropy |

---

## Troubleshooting

### Low pLDDT Scores

**Problem**: Predicted structure has large regions with pLDDT < 50
**Solutions**:
- Check if region is intrinsically disordered (use IUPred or similar)
- Try including homologous sequences in MSA (increase MSA depth)
- For viral proteins, check if region is hypervariable (may genuinely lack defined structure)
- Split protein into domains and predict separately

### FoldX Failures

**Problem**: FoldX returns unrealistic ΔΔG values (>10 kcal/mol)
**Solutions**:
- Ensure PDB file has been repaired (`RepairPDB` command)
- Check that mutation position exists in structure
- Verify chain ID matches
- For large changes, use Rosetta as cross-check

### Memory Issues

**Problem**: Structure prediction runs out of GPU memory
**Solutions**:
- Reduce sequence length (predict domains separately)
- Use CPU-only mode (slower but unlimited memory)
- Use smaller batch sizes for MSA retrieval
- Consider ESMFold for initial screening (lower memory)

---

## Examples

### Example 1: SARS-CoV-2 Spike D614G Analysis

```python
# Wild-type sequence (Spike S1 domain excerpt)
wt_seq = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHV..."

# Analyze D614G (position 614, Asp -> Gly)
report = await analyze_mutation_structural_impact(
    wild_type_sequence=wt_seq,
    mutation_position=614,
    wild_type_aa="D",
    mutant_aa="G",
    protein_name="SARS-CoV-2_Spike"
)

# Expected findings:
# - Stability: Slightly stabilizing (ddG ~ -0.5 kcal/mol)
# - Environment: Surface exposed, near protomer interface
# - Impact: Increases trimer stability, enhances infectivity
# - Structural change: Minimal local RMSD, improved packing
```

### Example 2: Identifying Druggable Sites in Viral Polymerase

```python
# Predict structure of viral RNA-dependent RNA polymerase
rdRp_seq = "..."  # Full sequence
rdRp_structure = predict_structure_alphafold(rdRp_seq, "RdRp")

# Identify binding pockets
pockets = identify_binding_sites(rdRp_structure.pdb_path)

# Filter for drug-like pockets
druggable = [p for p in pockets if p.druggability_score > 0.7 and p.volume > 200]

# Check overlap with active site
active_site_residues = [456, 457, 458, 759, 760, 761]  # Example
for pocket in druggable:
    overlap = set(pocket.residues) & set(active_site_residues)
    if overlap:
        print(f"Pocket {pocket.id} overlaps with active site: {overlap}")
```

---

## References

1. **AlphaFold**: Jumper et al., Nature 2021 (doi:10.1038/s41586-021-03819-2)
2. **ColabFold**: Mirdita et al., Nature Methods 2022 (doi:10.1038/s41592-022-01488-1)
3. **ESMFold**: Lin et al., Science 2023 (doi:10.1126/science.ade2574)
4. **FoldX**: Schymkowitz et al., Nucleic Acids Res 2005 (doi:10.1093/nar/gki379)
5. **FPocket**: Le Guilloux et al., BMC Bioinformatics 2009 (doi:10.1186/1471-2105-10-168)
6. **TM-align**: Zhang & Skolnick, Nucleic Acids Res 2005 (doi:10.1093/nar/gki524)

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2024-01-15 | Initial release with AlphaFold2/ColabFold support |
| 1.1 | 2024-03-20 | Added AlphaFold-Multimer for complexes |
| 1.2 | 2024-06-10 | Added ESMFold fast prediction option |
| 1.3 | 2024-09-05 | Integrated FoldX stability prediction |

---

## Related Skills

- **sequence-analysis.md**: Multiple sequence alignment, variant calling
- **molecular-docking.md**: Protein-ligand interaction prediction
- **phylogenetics.md**: Evolutionary analysis of viral variants
- **drug-discovery.md**: Virtual screening and lead optimization
