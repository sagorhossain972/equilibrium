#!/usr/bin/env python3
"""
build_seed_graph.py  --  GitHub issue #15
Ingest MONDO + HPO for the seed slice (diseases, phenotypes, synonyms).

Reads only two files:
    mondo.json       (MONDO ontology, OBO "graphs" JSON)
    phenotype.hpoa   (HPO disease->phenotype annotations, tab-separated)

Writes a dated snapshot folder:
    <out>/snapshots/<YYYY-MM-DD>/seed_graph.json
    <out>/snapshots/<YYYY-MM-DD>/diseases.csv
    <out>/snapshots/<YYYY-MM-DD>/phenotypes.csv
    <out>/snapshots/<YYYY-MM-DD>/edges.csv

Rules: never guess. Anything not present in the source files is null.
"""
import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

# --------------------------------------------------------------------------
# STEP 1. The seed list.
# MONDO IDs were found by searching MONDO labels (see the explanation) and are
# pinned here. `expect_label` is a safety check: if a future MONDO release
# changes what this ID means, the script stops instead of silently continuing.
# --------------------------------------------------------------------------
SEEDS = [
    {"mondo_id": "MONDO:0009266", "expect_label": "Gaucher disease type II",
     "requested_name": "Gaucher disease type 2"},
    {"mondo_id": "MONDO:0009267", "expect_label": "Gaucher disease type III",
     "requested_name": "Gaucher disease type 3"},
    {"mondo_id": "MONDO:0005180", "expect_label": "Parkinson disease",
     "requested_name": "Parkinson disease"},
]

# OPTIONAL (off by default, enable with --include-late-onset-pd).
# HPOA has annotations for OMIM:168600, which MONDO assigns to
# "late-onset Parkinson disease" (a CHILD of Parkinson disease), not to the
# Parkinson disease anchor itself. We never attach them to the anchor.
OPTIONAL_SEEDS = [
    {"mondo_id": "MONDO:0008199", "expect_label": "late-onset Parkinson disease",
     "requested_name": None, "parent_mondo_id": "MONDO:0005180"},
]

HPO_ID_RE = re.compile(r"^HP:\d{7}$")
HPOA_COLUMNS = ["database_id", "disease_name", "qualifier", "hpo_id", "reference",
                "evidence", "onset", "frequency", "sex", "modifier", "aspect",
                "biocuration"]


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------
def iri_to_curie(iri):
    """'http://purl.obolibrary.org/obo/MONDO_0009266' -> 'MONDO:0009266'"""
    return iri.rsplit("/", 1)[-1].replace("_", ":", 1)


def normalize_xref(xref):
    """MONDO writes 'Orphanet:77260'; HPOA writes 'ORPHA:77260'. Make them match."""
    return "ORPHA:" + xref.split(":", 1)[1] if xref.startswith("Orphanet:") else xref


def blank_to_none(value):
    value = (value or "").strip()
    return value if value else None


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def frequency_kind(freq):
    """Describe the FORMAT of the raw HPOA frequency value (we never reinterpret it)."""
    if freq is None:
        return None
    if freq.startswith("HP:"):
        return "hpo_term"
    if re.fullmatch(r"\d+/\d+", freq):
        return "fraction"
    if re.fullmatch(r"[\d.]+%", freq):
        return "percentage"
    return "unknown"


# --------------------------------------------------------------------------
# STEP 2. Read mondo.json
# --------------------------------------------------------------------------
def load_mondo(path):
    with open(path, encoding="utf-8") as f:
        graph = json.load(f)["graphs"][0]
    nodes = {iri_to_curie(n["id"]): n for n in graph["nodes"]}
    version = graph.get("meta", {}).get("version")
    # HP terms embedded in MONDO give us SOME phenotype labels (HPOA has none).
    hp_labels = {cid: n.get("lbl") for cid, n in nodes.items()
                 if cid.startswith("HP:") and n.get("lbl")}
    return nodes, hp_labels, version


# --------------------------------------------------------------------------
# STEP 3. Read phenotype.hpoa  (skip '#' comment lines, keep line numbers)
# --------------------------------------------------------------------------
def load_hpoa(path):
    version = None
    by_record = defaultdict(list)
    with open(path, encoding="utf-8") as f:
        for line_no, line in enumerate(f, start=1):
            line = line.rstrip("\n")
            if line.startswith("#"):
                if line.startswith("#version:"):
                    version = line.split(":", 1)[1].strip()
                continue
            if line.startswith("database_id"):      # header row
                continue
            parts = line.split("\t")
            parts += [""] * (len(HPOA_COLUMNS) - len(parts))
            row = dict(zip(HPOA_COLUMNS, parts))
            row["_line"] = line_no
            by_record[row["database_id"]].append(row)
    return by_record, version


# --------------------------------------------------------------------------
# STEP 4. Build one disease node from a MONDO node
# --------------------------------------------------------------------------
def build_disease_node(seed, mondo_nodes, hpoa_by_record, mondo_version):
    node = mondo_nodes.get(seed["mondo_id"])
    if node is None:
        raise SystemExit(f"{seed['mondo_id']} not found in mondo.json")
    if node.get("lbl") != seed["expect_label"]:
        raise SystemExit(f"{seed['mondo_id']} label is {node.get('lbl')!r}, "
                         f"expected {seed['expect_label']!r}. Re-check the seed.")
    if node.get("meta", {}).get("deprecated"):
        raise SystemExit(f"{seed['mondo_id']} is deprecated in MONDO")

    meta = node.get("meta", {})
    xrefs = [normalize_xref(x["val"]) for x in meta.get("xrefs", [])]

    def by_prefix(prefix):
        hits = sorted(x for x in xrefs if x.startswith(prefix + ":"))
        return hits or None                        # null if absent, never guessed

    synonyms = []
    for s in meta.get("synonyms", []):
        synonyms.append({
            "value": s["val"],
            "scope": s["pred"].replace("has", "").replace("Synonym", "").upper(),
            "type": (s.get("synonymType") or "").rsplit("#", 1)[-1] or None,
        })

    # HPOA records that this MONDO node explicitly points to via its own xrefs.
    # Only OMIM: and ORPHA: IDs can appear in HPOA. No fuzzy name matching.
    candidate_ids = sorted(x for x in xrefs if x.startswith(("OMIM:", "ORPHA:")))
    matched_ids = [x for x in candidate_ids if x in hpoa_by_record]

    return {
        "mondo_id": seed["mondo_id"],
        "name": node["lbl"],
        "requested_name": seed["requested_name"],
        "definition": (meta.get("definition") or {}).get("val"),
        "synonyms": synonyms,
        "external_ids": {
            "omim": by_prefix("OMIM"),
            "omim_ps": by_prefix("OMIMPS"),
            "orpha": by_prefix("ORPHA"),
        },
        "other_xrefs": sorted(x for x in xrefs
                              if not x.startswith(("OMIM:", "OMIMPS:", "ORPHA:"))),
        "parent_mondo_id": seed.get("parent_mondo_id"),
        "hpoa_record_ids": matched_ids,            # [] = no HPOA record matched
        "source": "MONDO",
        "source_version": mondo_version,
    }


# --------------------------------------------------------------------------
# STEP 5. Build disease -> phenotype edges from the matched HPOA rows
# --------------------------------------------------------------------------
def build_edges(disease, hpoa_by_record, hpoa_version, retrieved_at):
    edges, skipped_not = [], 0
    for record_id in disease["hpoa_record_ids"]:
        for row in hpoa_by_record[record_id]:
            qualifier = blank_to_none(row["qualifier"])
            if qualifier == "NOT":
                # "NOT" means the phenotype is explicitly ABSENT in this disease.
                # It must not become a HAS_PHENOTYPE edge.
                skipped_not += 1
                continue
            freq = blank_to_none(row["frequency"])
            edges.append({
                "source_id": disease["mondo_id"],
                "target_id": row["hpo_id"],
                "relation": "HAS_PHENOTYPE",
                "source": "HPO",
                # No official per-row ID exists in HPOA, so this points to the
                # exact row: <HPOA record>|<HPO term>|L<line in phenotype.hpoa>
                "source_record_id": f"{record_id}|{row['hpo_id']}|L{row['_line']}",
                "retrieved_at": retrieved_at,
                "source_version": hpoa_version,
                "hpoa_database_id": record_id,
                "hpoa_disease_name": blank_to_none(row["disease_name"]),
                "qualifier": qualifier,
                "reference": blank_to_none(row["reference"]),
                "evidence": blank_to_none(row["evidence"]),
                "onset": blank_to_none(row["onset"]),
                "frequency": freq,
                "frequency_kind": frequency_kind(freq),
                "sex": blank_to_none(row["sex"]),
                "modifier": blank_to_none(row["modifier"]),
                "aspect": blank_to_none(row["aspect"]),
                "biocuration": blank_to_none(row["biocuration"]),
            })
    return edges, skipped_not


# --------------------------------------------------------------------------
# STEP 6. Write JSON + CSV
# --------------------------------------------------------------------------
def _cell(value):
    if value is None:
        return ""
    if isinstance(value, list):
        return "|".join(str(v) for v in value)
    return value


def write_csv(path, rows, columns):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        for r in rows:
            writer.writerow([_cell(r.get(c)) for c in columns])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mondo", default="data/raw/mondo.json")
    ap.add_argument("--hpoa", default="data/raw/phenotype.hpoa")
    ap.add_argument("--out", default="data")
    ap.add_argument("--snapshot-date", default=datetime.now(timezone.utc).strftime("%Y-%m-%d"))
    ap.add_argument("--retrieved-at", default=None,
                    help="ISO timestamp of when the raw files were downloaded "
                         "(default: now, UTC)")
    ap.add_argument("--include-late-onset-pd", action="store_true",
                    help="also add MONDO:0008199 (late-onset PD) as a separate node")
    args = ap.parse_args()

    retrieved_at = args.retrieved_at or datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    mondo_nodes, hp_labels, mondo_version = load_mondo(args.mondo)
    hpoa_by_record, hpoa_version = load_hpoa(args.hpoa)

    seeds = SEEDS + (OPTIONAL_SEEDS if args.include_late_onset_pd else [])
    diseases = [build_disease_node(s, mondo_nodes, hpoa_by_record, mondo_version)
                for s in seeds]

    edges, skipped_not = [], 0
    for d in diseases:
        e, n = build_edges(d, hpoa_by_record, hpoa_version, retrieved_at)
        edges += e
        skipped_not += n

    phenotypes = [{
        "hpo_id": hid,
        "label": hp_labels.get(hid),               # null if MONDO doesn't carry it
        "label_source": "MONDO (embedded HP node)" if hid in hp_labels else None,
    } for hid in sorted({e["target_id"] for e in edges})]

    graph = {
        "meta": {
            "snapshot_date": args.snapshot_date,
            "retrieved_at": retrieved_at,
            "sources": {
                "MONDO": {"file": Path(args.mondo).name, "version": mondo_version,
                          "sha256": sha256_of(args.mondo)},
                "HPO": {"file": Path(args.hpoa).name, "version": hpoa_version,
                        "sha256": sha256_of(args.hpoa)},
            },
            "counts": {
                "diseases": len(diseases), "phenotypes": len(phenotypes),
                "edges": len(edges), "negated_rows_skipped": skipped_not,
                "phenotypes_without_label": sum(p["label"] is None for p in phenotypes),
            },
        },
        "nodes": {"diseases": diseases, "phenotypes": phenotypes},
        "edges": edges,
    }

    out_dir = Path(args.out) / "snapshots" / args.snapshot_date
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "seed_graph.json", "w", encoding="utf-8") as f:
        json.dump(graph, f, indent=2, ensure_ascii=False)

    write_csv(out_dir / "diseases.csv",
              [{**d, "omim": d["external_ids"]["omim"], "omim_ps": d["external_ids"]["omim_ps"],
                "orpha": d["external_ids"]["orpha"],
                "synonyms": [s["value"] for s in d["synonyms"]]} for d in diseases],
              ["mondo_id", "name", "requested_name", "omim", "omim_ps", "orpha",
               "hpoa_record_ids", "synonyms", "source_version"])
    write_csv(out_dir / "phenotypes.csv", phenotypes, ["hpo_id", "label", "label_source"])
    write_csv(out_dir / "edges.csv", edges,
              ["source_id", "target_id", "relation", "source", "source_record_id",
               "retrieved_at", "source_version", "hpoa_database_id", "qualifier",
               "reference", "evidence", "onset", "frequency", "frequency_kind",
               "sex", "modifier", "aspect", "biocuration"])

    print(f"Wrote snapshot to {out_dir}")
    print(json.dumps(graph["meta"]["counts"], indent=2))


if __name__ == "__main__":
    main()
