# autoDIAS 2.0

This repository now contains:

- `autoDIAS_1.0/`: the original script-based implementation preserved as a reference
- `src/autodias/`: the new modular implementation under development

The current `2.0` MVP focuses on replicating the core workflow with a cleaner architecture:

- TOML configuration
- RDKit-centered structure and fragment handling
- XYZ and Gaussian structure ingestion
- geometry extraction for distances, angles, and dihedrals
- optional frame reordering and reduction
- modular engine adapters
- Gaussian input preparation, execution, and energy analysis

ORCA-specific job generation and parsing will be added on top of the adapter layer next.

## Install

```bash
pip install -e .
```

## Quick Start

```bash
autodias prepare examples/autodias.toml
autodias run examples/autodias.toml
autodias analyze examples/autodias.toml
```

## Configuration

See [examples/autodias.toml](/c:/Users/Dennis/repos/autoDIAS/examples/autodias.toml) for a working template.
