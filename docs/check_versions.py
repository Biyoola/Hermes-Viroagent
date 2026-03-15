#!/usr/bin/env python3
try:
    import Bio
    print(f"BioPython version: {Bio.__version__}")
except ImportError:
    print("BioPython not found")
except AttributeError:
    print("BioPython imported but version attribute missing")

try:
    import MDAnalysis
    print(f"MDAnalysis version: {MDAnalysis.__version__}")
except ImportError:
    print("MDAnalysis not found")
except AttributeError:
    print("MDAnalysis imported but version attribute missing")

print("All required packages checked.")