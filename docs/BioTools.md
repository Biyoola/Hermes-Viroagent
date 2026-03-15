
---
name: setup_bioinformatics_applications_no_gpu
description: Install bioinformatics applications and tools that do not require GPU support for sequence analysis, phylogenetics, molecular docking, and visualization.
version: 1.0
category: environment_setup

inputs:
  install_method:
    type: string
    description: Package manager to use for installation (`conda` or `pip`)
    default: conda
  target_environment:
    type: string
    description: Name of the conda environment to install tools into
    default: bioinfo_env

outputs:
  installed_tools:
    type: list
    description: List of tools successfully installed
  errors:
    type: list
    description: Tools that failed to install

tools:
  - conda
  - pip
  - git
  - logging

execution_strategy: sequential
error_handling:
  - skip_failed_install
  - retry_failed_install
  - log_all_results

tags:
  - bioinformatics
  - environment
  - setup
---

# Skill: Setup Bioinformatics Applications (No GPU)

## Overview
Automates installation of bioinformatics applications that do not require GPU support.

## Applications to Install

| Application             | Type       | Purpose                      | Installation                           |
| ----------------------- | ---------- | ---------------------------- | -------------------------------------- |
| **MAFFT**               | CLI        | Multiple sequence alignment  | `conda install -c bioconda mafft`      |
| **BLAST+**              | CLI        | Sequence similarity search   | `conda install -c bioconda blast`      |
| **IQ-TREE**             | CLI        | Phylogenetic tree inference  | `conda install -c bioconda iqtree`     |
| **AutoDock Vina**       | CLI        | Molecular docking            | `conda install -c conda-forge vina`    |
| **RDKit**               | Python     | Cheminformatics              | `conda install -c conda-forge rdkit`   |
| **Nextstrain/Augur**    | Python/CLI | Viral phylogenetics          | `conda install -c bioconda augur`      |
| **BioPython**           | Python     | Sequence manipulation        | `pip install biopython`                |
| **MDAnalysis**          | Python     | Structure analysis           | `pip install mdanalysis`               |
| **PyMOL/NGLView**       | GUI/Python | 3D visualization             | `conda install -c conda-forge nglview` |

## Steps

1. Create or activate target environment

```bash
conda create -n {{target_environment}} python=3.11 -y
conda activate {{target_environment}}
```

2. Install CLI tools via conda

3. Install Python packages via pip if required

4. Verify installations

5. Log results: `installed_tools` and `errors`

## Example Input

```json
{
  "install_method": "conda",
  "target_environment": "bioinfo_env"
}
```

## Expected Output

```json
{
  "installed_tools": ["MAFFT","BLAST+","IQ-TREE","AutoDock Vina","RDKit","Nextstrain/Augur","BioPython","MDAnalysis","PyMOL/NGLView"],
  "errors": []
}
```

## Notes

- Excludes GPU-dependent tools such as ColabFold/AlphaFold and OpenMM
- Prepares environment for bioinformatics analyses without GPU requirements
