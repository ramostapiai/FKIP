# Functional Knowledge Integration Pipeline (FKIP) v8.2

Deterministic and fully traceable harmonization of functional metagenomic results from COG, KEGG, and FOAM.

## Description
FKIP harmonizes and integrates functional metagenomic statistical results at the terminal-feature level. It decouples molecular functions from biological processes and applies specificity-first criteria to prioritize functional processes in microbiome studies.

## Installation

```bash
git clone https://github.com/ramostapiai/FKIP.git
cd FKIP
conda env create -f environment.yml
conda activate fkip_env
```

Alternatively, you can install the dependencies using `pip`:

```bash
pip install -r requirements.txt
```

## Quick Start / Usage

To run the pipeline with your input data and reference resources:

```bash
python scripts/fkip_pipeline.py --input-dir inputs --resource-dir resources --output-dir outputs
```

To run a test with the provided mock dataset:

```bash
python scripts/fkip_pipeline.py --input-dir examples --resource-dir resources --output-dir outputs
```

## License
MIT License. See [LICENSE](LICENSE) for details.