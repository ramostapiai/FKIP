# Functional Knowledge Integration Pipeline (FKIP) v8.2

Deterministic and fully traceable harmonization of functional metagenomic results from COG, KEGG, and FOAM.

## Descripción
FKIP integra resultados estadísticos de metagenómica a nivel de características terminales, desacoplando funciones moleculares de procesos biológicos y aplicando criterios de especificidad para priorización de procesos en estudios de microbioma.

## Instalación

```bash
git clone https://github.com/ramostapiai/FKIP.git
cd FKIP
conda env create -f environment.yml
conda activate fkip_env
```

## Uso

```bash
python scripts/fkip_pipeline.py --input-dir inputs --resource-dir resources --output-dir outputs
```

## Licencia
MIT License. Ver [LICENSE](LICENSE) para más detalles.
# FKIP
