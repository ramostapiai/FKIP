#!/usr/bin/env python3
"""
Functional Knowledge Integration Pipeline (FKIP) v8.2
Harmonization of functional metagenomic results from COG, KEGG and FOAM.
"""

import os
import re
import io
import csv
import json
import math
import time
import zipfile
import hashlib
import warnings
import unicodedata
import argparse
from pathlib import Path
from datetime import datetime, timezone
from typing import Dict, List, Optional, Tuple, Iterable
from collections import defaultdict, Counter

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

warnings.filterwarnings("ignore")
pd.set_option("display.max_columns", 200)
pd.set_option("display.max_colwidth", 180)

def detect_encoding(path: Path) -> str:
    sample = Path(path).read_bytes()[:200_000]
    for enc in ("utf-8-sig", "utf-8", "latin-1", "cp1252"):
        try:
            sample.decode(enc)
            return enc
        except UnicodeDecodeError:
            pass
    return "utf-8"

def read_text(path: Optional[Path]) -> str:
    if path is None or not Path(path).exists():
        return ""
    return Path(path).read_text(encoding=detect_encoding(path), errors="replace")

def read_table(path: Path, sep=",", required_columns=None, **kwargs) -> pd.DataFrame:
    params = dict(sep=sep, encoding=detect_encoding(path), low_memory=False)
    params.update(kwargs)
    try:
        df = pd.read_csv(path, engine="c", **params)
    except Exception:
        params.pop("low_memory", None)
        df = pd.read_csv(path, engine="python", on_bad_lines="warn", **params)
    df.columns = [str(c).strip() for c in df.columns]
    if required_columns:
        missing = set(required_columns).difference(df.columns)
        if missing:
            raise ValueError(f"{Path(path).name}: missing columns {sorted(missing)}")
    return df

def clean_text(value) -> str:
    return "" if pd.isna(value) else str(value).strip()

def canonical_text(value) -> str:
    text = unicodedata.normalize("NFKC", clean_text(value)).casefold()
    text = re.sub(r"\[(?:PATH|BR|MODULE|KO|RN|mD)\s*:[^\]]+\]", "", text, flags=re.I)
    text = re.sub(r"\([^)]*\b(?:PATH|BR|MODULE|KO|RN)\s*:[^)]*\)", "", text, flags=re.I)
    text = text.replace("_", " ")
    text = re.sub(r"[^a-z0-9+./'-]+", " ", text)
    return re.sub(r"\s+", " ", text).strip(" -:;,.")

def ordered_unique(values: Iterable) -> List[str]:
    seen, out = set(), []
    for value in values:
        value = clean_text(value)
        if value and value not in seen:
            seen.add(value)
            out.append(value)
    return out

def split_multi(value: str) -> List[str]:
    return ordered_unique(re.split(r"\s*\|\|\s*|\s*;\s*", clean_text(value)))

def normalize_database(value) -> str:
    text = clean_text(value).upper()
    if "COG" in text:
        return "COG"
    if "KEGG" in text or text == "KO":
        return "KEGG"
    if "FOAM" in text:
        return "FOAM"
    return text

def normalize_level(value) -> str:
    text = clean_text(value)
    m = re.search(r"level\s*[-_ ]?\s*(id|[1-4])", text, flags=re.I)
    return f"Level-{m.group(1).upper()}" if m else text

def extract_explicit_id(name: str, database: str) -> str:
    text = clean_text(name)
    patterns = {
        "COG": [r"\bCOG\d{4,5}\b"],
        "KEGG": [r"\bK\d{5}\b", r"\bR\d{5}\b", r"\bM\d{5}\b"],
        "FOAM": [r"\[\s*mD\s*:\s*([A-Za-z0-9_.-]+)\s*\]", r"\bK\d{5}\b", r"\bm\d{5}\b"],
    }
    for pattern in patterns.get(database, []):
        match = re.search(pattern, text, flags=re.I)
        if match:
            return (match.group(1) if match.lastindex else match.group(0)).upper()
    return ""

def extract_ecs(value: str) -> List[str]:
    return sorted(set(re.findall(r"\b\d+\.(?:\d+|-|n)\.(?:\d+|-|n)\.(?:\d+|-|n)\b", clean_text(value))))

UNINFORMATIVE = {
    "", "function unknown", "general function prediction only", "poorly characterized",
    "not included in pathway or brite", "not included in regular maps", "unclassified",
    "unknown", "other", "hypothetical protein", "uncharacterized protein",
}

HOST_CONTEXT_PATTERNS = [
    r"human disease", r"cancer", r"neurodegenerative disease", r"cardiovascular disease",
    r"endocrine and metabolic disease", r"infectious disease", r"immune system",
    r"nervous system", r"circulatory system", r"endocrine system", r"digestive system",
    r"meiosis", r"apoptosis", r"ribosome biogenesis in eukaryotes",
    r"protein processing in endoplasmic reticulum", r"autophagy - animal", r"autophagy - yeast",
    r"alzheimer", r"amyotrophic lateral sclerosis", r"huntington", r"parkinson",
]

THEME_RULES = [
    (r"translation|ribosom|trna|aminoacyl", "Translation and ribosome biogenesis"),
    (r"replication|recombination|dna repair|mismatch repair|homologous recombination", "DNA replication and repair"),
    (r"transcription|rna polymerase|rna processing|rna degradation", "Transcription and RNA regulation"),
    (r"methan|oxidative phosphorylation|hydrogenase|energy production|respiration", "Energy metabolism"),
    (r"amino acid|lysine|cysteine|methionine|arginine|proline|glutamate|aspartate|peptide", "Amino acid metabolism"),
    (r"carbohydrate|glycolysis|starch|sugar|glycan|pentose|butanoate|propanoate|lactose", "Carbohydrate metabolism"),
    (r"lipid|fatty acid|ether lipid|isoprenoid", "Lipid metabolism"),
    (r"purine|pyrimidine|nucleotide", "Nucleotide metabolism"),
    (r"cofactor|vitamin|folate|thiamine|riboflavin|porphyrin", "Cofactor and vitamin metabolism"),
    (r"polyketide|secondary metabolite|nonribosomal peptide|surfactin", "Secondary metabolism"),
    (r"transport|transporter|permease|channel|efflux|phosphotransferase system", "Transport and nutrient acquisition"),
    (r"two-component|signal transduction|quorum sensing|chemotaxis|regulator", "Environmental sensing and signaling"),
    (r"flagell|motility|biofilm|cellular community", "Cell motility and community behavior"),
    (r"defense|antimicrobial|drug resistance|crispr|stress response|detoxification|mercur", "Stress response and antimicrobial defense"),
    (r"cell cycle|cell division|chromosome partition|z-ring|ftsz", "Cell division and cellular organization"),
    (r"cell wall|membrane/envelope biogenesis|peptidoglycan", "Cell envelope biogenesis and maintenance"),
    (r"protein turnover|chaperone|posttranslational|proteolysis|peptidase", "Protein turnover and quality control"),
    (r"secretion|extracellular structure|intracellular trafficking", "Secretion and extracellular interaction"),
]

CORE_HOUSEKEEPING_PATTERNS = [
    r"^translation", r"ribosome", r"aminoacyl-trna", r"dna replication",
    r"rna polymerase", r"central carbon", r"glycolysis", r"oxidative phosphorylation",
    r"protein turnover", r"chaperone", r"basic transcription",
]

BROAD_PROCESS_PATTERNS = [
    r"^energy metabolism$", r"^amino acid metabolism$", r"^carbohydrate metabolism$",
    r"^lipid metabolism$", r"^nucleotide metabolism$", r"^transport and nutrient acquisition$",
    r"^environmental sensing and signaling$", r"^secondary metabolism$",
]

def theme_from_text(text: str) -> str:
    for pattern, theme in THEME_RULES:
        if re.search(pattern, clean_text(text), flags=re.I):
            return theme
    return ""

def is_host_text(text: str) -> Tuple[bool, str]:
    for pattern in HOST_CONTEXT_PATTERNS:
        if re.search(pattern, clean_text(text), flags=re.I):
            return True, pattern
    return False, ""

def first_informative(*values) -> str:
    for value in values:
        value = clean_text(value)
        if value and canonical_text(value) not in UNINFORMATIVE:
            return value
    return ""

def parse_function_fields(function: str, database: str, ec_value: str = "", gene_value: str = "") -> dict:
    native = clean_text(function)
    canonical = native
    gene_aliases = split_multi(gene_value)
    enzyme_aliases = []

    if database == "FOAM" and ";" in native:
        prefix, tail = native.split(";", 1)
        canonical = clean_text(tail)
        for token in [x.strip() for x in prefix.split(",")]:
            if re.fullmatch(r"E?\d+\.\d+\.\d+\.\d+", token, flags=re.I):
                enzyme_aliases.append(token.lstrip("Ee"))
            elif token and not token.lower().startswith("ec"):
                gene_aliases.append(token)

    canonical = re.sub(r"\[(?:mD|BR|PATH|MODULE|KO|RN)\s*:[^\]]+\]", "", canonical, flags=re.I)
    canonical = re.sub(r"^E\d+\.\d+\.\d+\.\d+\s*,?\s*", "", canonical, flags=re.I)
    canonical = re.sub(r"\s+", " ", canonical).strip(" ;,.-")
    normalized = canonical_text(canonical)

    low = normalized
    if not low or low in UNINFORMATIVE or re.search(r"uncharacterized|hypothetical|unknown function", low):
        semantic_class = "uncharacterized"
        informativeness = "uncharacterized"
    elif re.search(r"transport|transporter|permease|channel|efflux|importer|exporter", low):
        semantic_class = "specific_transporter"
        informativeness = "specific"
    elif re.search(r"regulator|transcription factor|response regulator|sensor kinase|quorum", low):
        semantic_class = "regulatory_protein"
        informativeness = "specific"
    elif re.search(r"cell division|z-ring|chromosome partition|ribosomal protein|chaperone", low):
        semantic_class = "cellular_system_protein"
        informativeness = "specific"
    elif extract_ecs(ec_value + " " + native) or re.search(r"ase\b|synthase|synthetase|dehydrogenase|kinase|hydrolase|reductase", low):
        semantic_class = "specific_enzyme"
        informativeness = "specific" if len(low.split()) >= 2 else "generic"
    elif low in {"atpase", "protein", "membrane protein", "enzyme"} or len(low.split()) <= 1:
        semantic_class = "generic_protein"
        informativeness = "generic"
    else:
        semantic_class = "described_protein"
        informativeness = "moderate"

    return {
        "native_function": native,
        "normalized_function": normalized,
        "canonical_function": canonical,
        "function_semantic_class": semantic_class,
        "function_informativeness": informativeness,
        "function_gene_aliases": " || ".join(ordered_unique(gene_aliases)),
        "function_enzyme_aliases": " || ".join(ordered_unique(enzyme_aliases + extract_ecs(ec_value))),
    }

def function_process_candidate(canonical_function: str, semantic_class: str) -> Tuple[str, str]:
    text = clean_text(canonical_function)
    low = canonical_text(text)
    if not low or semantic_class in {"uncharacterized", "generic_protein"}:
        return "", "FUNCTION_NOT_INFORMATIVE"

    if re.search(r"cell division protein zapb|interacts with ftsz|z-ring", low):
        return "Bacterial cell division and Z-ring organization", "FUNCTION_CELL_DIVISION_ZAPB"
    if re.search(r"mercuric reductase|mercury reductase|\bmera\b", low):
        return "Mercury detoxification", "FUNCTION_MERCURY_DETOXIFICATION"
    if "formaldehyde" in low:
        return "Formaldehyde detoxification and aldehyde oxidation", "FUNCTION_FORMALDEHYDE"
    if "quorum" in low or re.search(r"luxr|sdia", low):
        return "Quorum-sensing regulation", "FUNCTION_QUORUM_REGULATION"
    if re.search(r"crispr|cas protein", low):
        return "CRISPR-Cas adaptive defense", "FUNCTION_CRISPR"
    if re.search(r"trna.*modif|trna.*methyl|trm\w*", low):
        return "tRNA modification", "FUNCTION_TRNA_MODIFICATION"
    if re.search(r"dna repair|dna ligase|ribonuclease hi|dna polymerase", low):
        return "DNA replication and repair", "FUNCTION_DNA_REPAIR"
    if re.search(r"ribosomal protein|translation factor|initiation factor|elongation factor", low):
        return "Translation and ribosome biogenesis", "FUNCTION_TRANSLATION"
    if re.search(r"\bpts system\b|phosphotransferase system", low):
        substrate = re.sub(r"\b(?:pts|phosphotransferase) system.*$", "", low).strip(" ,;/-")
        if substrate:
            return substrate[:1].upper() + substrate[1:] + " phosphotransferase-system transport", "FUNCTION_PTS_TRANSPORT"
        return "Phosphotransferase-system transport", "FUNCTION_PTS_TRANSPORT"

    if semantic_class == "specific_transporter":
        cleaned = re.sub(r"\b(?:abc-type|abc type|pts|phosphotransferase system)\b", "", low)
        cleaned = re.sub(r"\b(?:transport system|transporter|transport|permease|channel|efflux pump|importer|exporter)\b.*$", "", cleaned)
        cleaned = re.sub(r"\b(?:atp-binding|atpase|substrate-binding|permease) component\b", "", cleaned)
        cleaned = cleaned.replace("-specific", "")
        substrate = re.sub(r"\s+", " ", cleaned).strip(" ,;/-")
        if substrate and 1 <= len(substrate.split()) <= 8:
            parts = [part.strip() for part in substrate.split("/") if part.strip()]
            if len(parts) > 1:
                substrate = ", ".join(parts[:-1]) + " and " + parts[-1]
            elif parts:
                substrate = parts[0]
            return substrate[:1].upper() + substrate[1:] + " transport", "FUNCTION_SPECIFIC_TRANSPORT"
        return "Specific membrane transport", "FUNCTION_GENERIC_TRANSPORT"

    if re.search(r"alcohol.*dehydrogenase|aldehyde.*dehydrogenase", low):
        return "Alcohol and aldehyde oxidation", "FUNCTION_ALCOHOL_ALDEHYDE"
    if re.search(r"surfactin synthetase|nonribosomal peptide", low):
        return "Nonribosomal peptide biosynthesis", "FUNCTION_NRP_BIOSYNTHESIS"
    if re.search(r"cell wall|peptidoglycan", low):
        return "Cell-envelope biogenesis and maintenance", "FUNCTION_CELL_ENVELOPE"

    return "", "FUNCTION_NO_CONTROLLED_PROCESS"

def process_specificity(process: str) -> Tuple[str, int, bool, bool]:
    text = canonical_text(process)
    if not text:
        return "unresolved", 0, False, False
    core = any(re.search(p, text, flags=re.I) for p in CORE_HOUSEKEEPING_PATTERNS)
    broad = any(re.search(p, text, flags=re.I) for p in BROAD_PROCESS_PATTERNS)
    if core:
        return "core_housekeeping", 1, True, broad
    if broad:
        return "broad_context", 1, False, True
    if re.search(r"specific membrane transport|abc-type transport|signal transduction|biosynthesis|degradation", text):
        return "functional_system", 3, False, False
    if re.search(r"transport|detoxification|quorum|cell division|z-ring|crispr|formaldehyde|mercury|phosphonate|sulfate|lactose", text):
        return "specific_process", 4, False, False
    return "intermediate_process", 3, False, False

CONFIDENCE_SCORE = {"High": 3, "Moderate": 2, "Low": 1, "Not assigned": 0, "Excluded": 0}

class UnionFind:
    def __init__(self, items):
        self.parent = {item: item for item in items}
    def find(self, item):
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item
    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra

def main():
    parser = argparse.ArgumentParser(description="FKIP v8.2 CLI Pipeline")
    parser.add_argument("--input-dir", type=Path, default=Path("inputs"), help="Directory with input matrices and master stats")
    parser.add_argument("--resource-dir", type=Path, default=Path("resources"), help="Directory with external resources")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"), help="Directory for output tables and summaries")
    parser.add_argument("--p-val-threshold", type=float, default=0.05, help="P-value threshold")
    parser.add_argument("--positive-group", type=str, default="Keratoconus", help="Name of positive group")
    parser.add_argument("--negative-group", type=str, default="Healthy siblings", help="Name of negative group")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    print(f"[{datetime.now().isoformat()}] Starting FKIP v8.2 pipeline...")
    print(f"Inputs: {args.input_dir.resolve()}")
    print(f"Resources: {args.resource_dir.resolve()}")
    print(f"Outputs: {args.output_dir.resolve()}")

    master_stats_path = args.input_dir / "LM_Functional_Master_Stats.csv"
    if not master_stats_path.exists():
        master_stats_path = args.input_dir / "inputs" / "LM_Functional_Master_Stats.csv"

    cog_rollup_path = args.input_dir / "native_rollups" / "combined_wide_COG_rollup.tsv"
    kegg_rollup_path = args.input_dir / "native_rollups" / "combined_wide_KEGG_rollup.tsv"
    foam_rollup_path = args.input_dir / "native_rollups" / "combined_wide_FOAM_rollup.tsv"

    print("Checking input tables...")
    for p, name in [(master_stats_path, "Master stats"), (cog_rollup_path, "COG rollup"), 
                    (kegg_rollup_path, "KEGG rollup"), (foam_rollup_path, "FOAM rollup")]:
        if p.exists():
            print(f"  [OK] {name}: {p}")
        else:
            print(f"  [WARN] Missing {name} at {p}")

    print("Pipeline ready for execution.")

if __name__ == "__main__":
    main()
