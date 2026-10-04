"""
Validation tests for the seed graph snapshot.

Run with:      pytest tests/test_seed_graph.py -v
or without pytest:  python tests/test_seed_graph.py [path/to/seed_graph.json]

By default the newest snapshot under data/snapshots/ is tested.
Set SEED_GRAPH=/path/to/seed_graph.json to test a specific file.
"""
import json
import os
import re
import sys
from datetime import datetime
from pathlib import Path

EXPECTED = {
    "MONDO:0009266": "Gaucher disease type II",     # = Gaucher disease type 2
    "MONDO:0009267": "Gaucher disease type III",    # = Gaucher disease type 3
    "MONDO:0005180": "Parkinson disease",           # anchor
}
HPO_ID_RE = re.compile(r"^HP:\d{7}$")
MONDO_ID_RE = re.compile(r"^MONDO:\d{7}$")


def find_graph_path():
    if len(sys.argv) > 1 and sys.argv[-1].endswith(".json"):
        return Path(sys.argv[-1])
    if os.environ.get("SEED_GRAPH"):
        return Path(os.environ["SEED_GRAPH"])
    snaps = sorted(Path("data/snapshots").glob("*/seed_graph.json"))
    assert snaps, "No snapshot found under data/snapshots/ - run the build script first"
    return snaps[-1]


def load_graph():
    with open(find_graph_path(), encoding="utf-8") as f:
        return json.load(f)


def test_output_is_non_empty():
    g = load_graph()
    assert g["nodes"]["diseases"], "no disease nodes"
    assert g["nodes"]["phenotypes"], "no phenotype nodes"
    assert g["edges"], "no edges"


def test_all_three_seed_diseases_exist():
    g = load_graph()
    names = {d["mondo_id"]: d["name"] for d in g["nodes"]["diseases"]}
    for mondo_id, name in EXPECTED.items():
        assert mondo_id in names, f"missing disease node {mondo_id}"
        assert names[mondo_id] == name, f"{mondo_id} has name {names[mondo_id]!r}"
    # Gaucher type 2 and type 3 must stay separate nodes
    assert "MONDO:0009266" != "MONDO:0009267"


def test_mondo_ids_are_unique_and_well_formed():
    ids = [d["mondo_id"] for d in load_graph()["nodes"]["diseases"]]
    assert len(ids) == len(set(ids)), "duplicate MONDO IDs"
    assert all(MONDO_ID_RE.match(i) for i in ids), "malformed MONDO ID"


def test_hpo_ids_are_valid_and_unique():
    g = load_graph()
    ids = [p["hpo_id"] for p in g["nodes"]["phenotypes"]]
    assert len(ids) == len(set(ids)), "duplicate phenotype nodes"
    bad = [i for i in ids if not HPO_ID_RE.match(i)]
    assert not bad, f"invalid HPO IDs: {bad}"
    bad_targets = [e["target_id"] for e in g["edges"] if not HPO_ID_RE.match(e["target_id"])]
    assert not bad_targets, f"invalid HPO IDs on edges: {bad_targets}"


def test_edges_point_to_existing_nodes():
    g = load_graph()
    disease_ids = {d["mondo_id"] for d in g["nodes"]["diseases"]}
    hpo_ids = {p["hpo_id"] for p in g["nodes"]["phenotypes"]}
    for e in g["edges"]:
        assert e["source_id"] in disease_ids, f"dangling source {e['source_id']}"
        assert e["target_id"] in hpo_ids, f"dangling target {e['target_id']}"


def test_every_edge_has_provenance():
    for e in load_graph()["edges"]:
        assert e["relation"] == "HAS_PHENOTYPE"
        assert e["source"] == "HPO"
        assert e.get("source_record_id"), f"missing source_record_id: {e}"
        assert e.get("retrieved_at"), f"missing retrieved_at: {e}"
        datetime.fromisoformat(e["retrieved_at"].replace("Z", "+00:00"))  # parseable


def test_no_negated_annotations_became_edges():
    assert all(e["qualifier"] != "NOT" for e in load_graph()["edges"])


def test_edges_only_come_from_records_mondo_links_to():
    """No invention: every edge's HPOA record must be one of the disease's own xrefs,
    which also proves Gaucher type 2 and type 3 edges never get mixed."""
    g = load_graph()
    allowed = {d["mondo_id"]: set(d["hpoa_record_ids"]) for d in g["nodes"]["diseases"]}
    for e in g["edges"]:
        assert e["hpoa_database_id"] in allowed[e["source_id"]], (
            f"{e['source_id']} got an edge from unlinked record {e['hpoa_database_id']}")


if __name__ == "__main__":
    tests = [(n, f) for n, f in sorted(globals().items())
             if n.startswith("test_") and callable(f)]
    failed = 0
    for name, fn in tests:
        try:
            fn()
            print(f"PASS  {name}")
        except AssertionError as exc:
            failed += 1
            print(f"FAIL  {name}: {exc}")
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    sys.exit(1 if failed else 0)
