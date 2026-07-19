#!/usr/bin/env python3
"""Build the isolated database snapshot used by the real-world SCI E2E matrix."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import yaml

from allenricher.database.animaltfdb_fetcher import AnimalTFDBFetcher
from allenricher.database.builder import DatabaseBuilder
from allenricher.database.chea3_fetcher import ChEA3Fetcher
from allenricher.database.custom_builder import CustomDatabaseBuilder
from allenricher.database.htftarget_fetcher import HTFtargetFetcher
from allenricher.database.kegg_fetcher import KEGGFetcher
from allenricher.database.manager import DatabaseManager
from allenricher.database.parsers.animaltfdb import AnimalTFDBParser
from allenricher.database.parsers.do import DOParser
from allenricher.database.parsers.reactome import ReactomeParser
from allenricher.database.trrust_fetcher import TRRUSTFetcher


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MATRIX = Path(__file__).with_name("case_matrix.yaml")
DEFAULT_DATABASE_ROOT = (
    PROJECT_ROOT / "test_e2e_2026" / "00_input_data" / "real_world_sci" / "database_snapshot"
)
DEFAULT_SOURCE_BASIC = PROJECT_ROOT / "database" / "basic"
V1_DISGENET_ROOT = PROJECT_ROOT.parent / "AllEnricher-v1" / "database" / "organism" / "v20190612" / "hsa"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def latest_version(root: Path, prefix: str) -> str:
    versions = sorted(path.name for path in root.iterdir() if path.is_dir() and path.name.startswith(prefix))
    if not versions:
        raise FileNotFoundError(f"no {prefix} source version under {root}")
    return versions[-1]


def _read_disc(path: Path) -> dict[str, str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return {
            parts[0]: parts[1]
            for parts in (line.rstrip("\r\n").split("\t") for line in handle)
            if len(parts) >= 2
        }


def _rewrite_gmt_descriptions(path: Path, names: dict[str, str]) -> None:
    rows = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) >= 3 and parts[0] in names:
                parts[1] = names[parts[0]]
                rows.append(parts)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for parts in rows:
            handle.write("\t".join(parts) + "\n")


def refresh_existing_hierarchy_metadata(
    database_root: Path,
    source_basic: Path,
    versions: dict,
    species: str,
) -> None:
    species_dirs = sorted((database_root / "organism").glob(f"v*/{species}"))
    if not species_dirs:
        return

    kegg_dir = next(
        (
            path for path in reversed(species_dirs)
            if (path / f"{species}.kegg2disc.gz").is_file()
        ),
        species_dirs[-1],
    )
    kegg_disc = kegg_dir / f"{species}.kegg2disc.gz"
    kegg_gmt = kegg_dir / f"{species}.KEGG.gmt.gz"
    if kegg_disc.is_file() and kegg_gmt.is_file():
        existing_terms = _read_disc(kegg_disc)
        fetcher = KEGGFetcher(str(source_basic / "kegg"))
        pathways = fetcher._list_pathways(species)
        categories = fetcher._get_brite_categories(species, pathways)
        clean_names = {
            f"{species}{pathway_id}": fetcher._clean_pathway_name(name)
            for pathway_id, name in pathways
        }
        with gzip.open(kegg_disc, "wt", encoding="utf-8") as handle:
            for term_id in sorted(existing_terms):
                pathway_id = term_id.removeprefix(species)
                name = clean_names.get(term_id, term_id)
                category, subcategory = categories.get(
                    pathway_id, ("Uncategorized", "Uncategorized")
                )
                hierarchy = ""
                if category != "Uncategorized":
                    hierarchy = "|".join(
                        part.replace("_", " ") for part in (category, subcategory, name)
                    )
                handle.write(f"{term_id}\t{name}\t{hierarchy}\n")
        _rewrite_gmt_descriptions(kegg_gmt, clean_names)

    reactome_dir = next(
        (
            path for path in reversed(species_dirs)
            if (path / f"{species}.Reactome2disc.gz").is_file()
        ),
        species_dirs[-1],
    )
    reactome_disc = reactome_dir / f"{species}.Reactome2disc.gz"
    reactome_gmt = reactome_dir / f"{species}.Reactome.gmt.gz"
    if reactome_disc.is_file() and reactome_gmt.is_file():
        source_dir = source_basic / "reactome" / versions["reactome"]
        fallback = _read_disc(reactome_disc)
        names, hierarchies = ReactomeParser._load_hierarchies(
            str(source_dir / "ReactomePathways.txt"),
            str(source_dir / "ReactomePathwaysRelation.txt"),
            set(fallback),
            fallback,
        )
        with gzip.open(reactome_disc, "wt", encoding="utf-8") as handle:
            for term_id in sorted(fallback):
                handle.write(
                    f"{term_id}\t{names.get(term_id, fallback[term_id])}\t"
                    f"{hierarchies.get(term_id, '')}\n"
                )
        _rewrite_gmt_descriptions(reactome_gmt, names)

    if species == "hsa":
        do_dir = next(
            (path for path in reversed(species_dirs) if (path / "hsa.DO2disc.gz").is_file()),
            species_dirs[-1],
        )
        do_disc = do_dir / "hsa.DO2disc.gz"
        do_gmt = do_dir / "hsa.DO.gmt.gz"
        if do_disc.is_file() and do_gmt.is_file():
            fallback = _read_disc(do_disc)
            names, hierarchies, obsolete = DOParser._load_ontology(
                str(source_basic / "do" / "doid.obo"), set(fallback)
            )
            valid_names = {
                term_id: names.get(term_id, fallback[term_id]).replace("_", " ")
                for term_id in fallback
                if term_id not in obsolete
                and names.get(term_id, fallback[term_id]).casefold() != term_id.casefold()
            }
            with gzip.open(do_disc, "wt", encoding="utf-8") as handle:
                for term_id in sorted(valid_names):
                    handle.write(
                        f"{term_id}\t{valid_names[term_id]}\t{hierarchies.get(term_id, '')}\n"
                    )
            _rewrite_gmt_descriptions(do_gmt, valid_names)


def fetch_tf_sources(
    basic_dir: Path,
    offline: bool,
    animal_species: set[str] | None = None,
) -> None:
    animal_species = animal_species or set()
    required = [
        basic_dir / "trrust" / "TRRUSTv2" / "trrust_rawdata.hsa.tsv",
        basic_dir / "htftarget" / "tf-target-information.txt",
        *(
            basic_dir / "animaltfdb" / "AnimalTFDBv4.0" / filename
            for latin_name in sorted(animal_species)
            for filename in (f"{latin_name}_TF", f"{latin_name}_ortholog_to_human")
        ),
    ]
    chea3_required = [
        basic_dir / "chea3" / "ChEA3v2024" / f"{name}_tf.gmt"
        for name in ChEA3Fetcher.CHEA3_GMT_LIBS
    ]
    missing = [path for path in [*required, *chea3_required] if not path.is_file() or not path.stat().st_size]
    if offline:
        if missing:
            raise FileNotFoundError("offline TF cache missing: " + ", ".join(map(str, missing)))
        return

    TRRUSTFetcher(str(basic_dir)).download_species("Homo sapiens")
    ChEA3Fetcher(str(basic_dir)).download_all_gmt_libraries()
    HTFtargetFetcher(str(basic_dir)).download()
    animal = AnimalTFDBFetcher(str(basic_dir))
    for latin_name in sorted(animal_species):
        outputs = animal.download_species_data(latin_name)
        if set(outputs) != {"tf_list", "ortholog"}:
            raise RuntimeError(f"AnimalTFDB download incomplete for {latin_name}: {outputs}")


def copy_v1_disgenet(outdir: Path) -> None:
    if not V1_DISGENET_ROOT.is_dir():
        raise FileNotFoundError(f"AllEnricher-v1 DisGeNET snapshot missing: {V1_DISGENET_ROOT}")
    for filename in ("hsa.CUI2gene.tab.gz", "hsa.CUI2disc.gz"):
        source = V1_DISGENET_ROOT / filename
        if not source.is_file():
            raise FileNotFoundError(source)
        shutil.copy2(source, outdir / filename)
    (outdir / "DISGENET_SOURCE_NOTICE.txt").write_text(
        "DisGeNET source: AllEnricher-v1 free snapshot v20190612.\n"
        "DisGeNET continues to update, but later full releases are not freely downloadable.\n",
        encoding="utf-8",
    )


def prepare_animaltfdb_gene_info(source: Path, target: Path, taxids: set[int]) -> Path:
    """Freeze only the taxa needed for AnimalTFDB ID conversion in one stream pass."""
    if target.is_file() and target.stat().st_size:
        found = set()
        with gzip.open(target, "rt", encoding="utf-8") as handle:
            for line in handle:
                if not line.startswith("#"):
                    found.add(line.partition("\t")[0])
        if {str(taxid) for taxid in taxids} <= found:
            return target
    target.parent.mkdir(parents=True, exist_ok=True)
    temp = target.with_suffix(target.suffix + ".part")
    wanted = {str(taxid) for taxid in taxids}
    with gzip.open(source, "rt", encoding="utf-8") as reader, gzip.open(
        temp, "wt", encoding="utf-8", newline=""
    ) as writer:
        for line in reader:
            if line.startswith("#") or line.partition("\t")[0] in wanted:
                writer.write(line)
    temp.replace(target)
    return target


def build_public_go_custom(database_root: Path, species: str, taxid: int, source_dir: Path) -> None:
    manager = DatabaseManager(str(database_root), species)
    manager.load_database("GO")
    terms = manager.get_all_term_data()["GO"]
    input_dir = database_root / "custom_inputs" / species
    input_dir.mkdir(parents=True, exist_ok=True)
    input_gmt = input_dir / f"{species}.PUBLIC_GO_CUSTOM.source.gmt.gz"
    with gzip.open(input_gmt, "wt", encoding="utf-8", newline="") as handle:
        for term_id, info in sorted(terms.items()):
            genes = sorted(set(map(str, info["genes"])))
            handle.write("\t".join([str(term_id), str(info.get("name") or term_id), *genes]) + "\n")
    CustomDatabaseBuilder(str(database_root)).build_from_gmt(
        gmt_file=str(input_gmt),
        species=species,
        taxid=taxid,
        db_name="PUBLIC_GO_CUSTOM",
    )


def validate_snapshot(database_root: Path, matrix: dict) -> dict:
    audit = {}
    for case_id, spec in matrix["datasets"].items():
        manager = DatabaseManager(str(database_root), spec["species"])
        manager.load_databases(spec["databases"])
        case_audit = {}
        for database, terms in manager.get_all_term_data().items():
            sizes = [len(set(info["genes"])) for info in terms.values()]
            passing = sum(10 <= size <= 500 for size in sizes)
            hierarchy_count = sum(bool(info.get("hierarchy")) for info in terms.values())
            case_audit[database] = {
                "gene_sets": len(terms),
                "passing_gsea_size_filter": passing,
                "genes": len({gene for info in terms.values() for gene in info["genes"]}),
                "hierarchy_terms": hierarchy_count,
            }
            if passing < 10:
                raise ValueError(f"{case_id}/{database} has only {passing} gene sets passing size filter")
            if database in {"GO", "KEGG", "Reactome", "DO"} and hierarchy_count == 0:
                raise ValueError(f"{case_id}/{database} is missing hierarchy metadata")
        audit[case_id] = case_audit
    return audit


def species_snapshot_ready(database_root: Path, spec: dict) -> bool:
    """Resume only when every configured database loads and passes the matrix gate."""
    try:
        validate_snapshot(database_root, {"datasets": {"resume_check": spec}})
    except (EOFError, FileNotFoundError, KeyError, OSError, ValueError):
        return False
    return True


def source_manifest(database_root: Path, source_basic: Path, versions: dict) -> dict:
    files = []
    source_paths = [
        source_basic / "go" / versions["go"] / "gene2go.gz",
        source_basic / "go" / versions["go"] / "gene_info.gz",
        source_basic / "go" / versions["go"] / "go-basic.obo",
        source_basic / "reactome" / versions["reactome"] / "NCBI2Reactome_All_Levels.txt.gz",
        source_basic / "reactome" / versions["reactome"] / "ReactomePathways.txt",
        source_basic / "reactome" / versions["reactome"] / "ReactomePathwaysRelation.txt",
        source_basic / "do" / "doid.obo",
        source_basic / "wikipathways" / f"WP{versions['wikipathways']}",
        database_root / "basic",
        database_root / "organism",
        database_root / "custom_inputs",
        V1_DISGENET_ROOT / "hsa.CUI2gene.tab.gz",
        V1_DISGENET_ROOT / "hsa.CUI2disc.gz",
    ]
    for source in source_paths:
        candidates = source.rglob("*") if source.is_dir() else [source]
        for path in candidates:
            if path.is_file():
                files.append({"path": str(path.resolve()), "size": path.stat().st_size, "sha256": sha256(path)})
    return {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "versions": versions,
        "files": files,
    }


def prepare(matrix_path: Path, database_root: Path, source_basic: Path, offline: bool) -> dict:
    matrix = yaml.safe_load(matrix_path.read_text(encoding="utf-8"))
    database_root.mkdir(parents=True, exist_ok=True)
    if offline:
        required = [
            database_root / "SOURCE_MANIFEST.json",
            database_root / "DATABASE_AUDIT.json",
            database_root / "organism",
        ]
        missing = [str(path) for path in required if not path.exists()]
        if missing:
            raise FileNotFoundError("offline database snapshot incomplete: " + ", ".join(missing))
        audit = validate_snapshot(database_root, matrix)
        (database_root / "DATABASE_AUDIT.json").write_text(
            json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )
        stored_manifest = json.loads(
            (database_root / "SOURCE_MANIFEST.json").read_text(encoding="utf-8")
        )
        refreshed_manifest = source_manifest(
            database_root, source_basic, stored_manifest.get("versions", {})
        )
        (database_root / "SOURCE_MANIFEST.json").write_text(
            json.dumps(refreshed_manifest, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        return audit

    tf_basic = database_root / "basic"
    animal_species = {
        spec["latin_name"].replace(" ", "_")
        for spec in matrix["datasets"].values()
        if "AnimalTFDB" in spec["databases"] and spec["species"] != "hsa"
    }
    fetch_tf_sources(tf_basic, offline, animal_species)
    versions = {
        "go": latest_version(source_basic / "go", "GO"),
        "reactome": latest_version(source_basic / "reactome", "Reactome"),
        "wikipathways": latest_version(source_basic / "wikipathways", "WP").removeprefix("WP"),
        "disgenet": "AllEnricher-v1-v20190612-free-snapshot",
        "trrust": "v2",
        "chea3": "v2024-frozen-local",
        "animaltfdb": "v4.0",
        "htftarget": "2020",
    }
    for spec in matrix["datasets"].values():
        refresh_existing_hierarchy_metadata(
            database_root,
            source_basic,
            versions,
            spec["species"],
        )
    animal_taxids = {
        9606,
        *(
            int(spec.get("annotation_taxid", spec["taxid"]))
            for spec in matrix["datasets"].values()
            if "AnimalTFDB" in spec["databases"] and spec["species"] != "hsa"
        ),
    }
    animal_gene_info = prepare_animaltfdb_gene_info(
        source_basic / "go" / versions["go"] / "gene_info.gz",
        tf_basic / "gene_info.animaltfdb_selected_taxa.gz",
        animal_taxids,
    )
    AnimalTFDBParser.load_external_id_symbol_maps(str(animal_gene_info), animal_taxids)

    for spec in matrix["datasets"].values():
        species = spec["species"]
        taxid = int(spec["taxid"])
        annotation_taxid = int(spec.get("annotation_taxid", taxid))
        if species_snapshot_ready(database_root, spec):
            print(f"[RESUME] {species}: All configuration databases can be loaded and jumped through reconstruction via size filter")
            continue
        builder = DatabaseBuilder(str(database_root))
        builder.basic_dir = source_basic
        generate_gmt_files = builder.generate_gmt_files
        builder.generate_gmt_files = lambda *_args, **_kwargs: {}
        builder.build_go(species, annotation_taxid, versions["go"])
        builder.build_kegg(species, annotation_taxid, versions["go"])
        builder.build_reactome(species, annotation_taxid, versions["reactome"])
        builder.build_wikipathways(species, annotation_taxid, versions["wikipathways"])
        outdir = Path(builder.organism_dir) / f"v{datetime.now().strftime('%Y%m%d')}" / species

        if species == "hsa":
            builder.build_do(taxid, versions["go"])
            copy_v1_disgenet(outdir)

        builder.basic_dir = tf_basic
        if species == "hsa":
            builder.build_trrust(species, taxid)
            builder.build_chea3(species, taxid)
            builder.build_htftarget(species, taxid)
            for stale_name in (
                "hsa.AnimalTFDB_2gene.tab.gz",
                "hsa.AnimalTFDB_mapped_2disc.gz",
                "ANIMALTFDB_HUMAN_SOURCE_NOTICE.txt",
            ):
                (outdir / stale_name).unlink(missing_ok=True)
        elif "AnimalTFDB" in spec["databases"]:
            builder.build_animaltfdb(
                species,
                taxid,
                species_latin=spec["latin_name"].replace(" ", "_"),
                gene_info_path=str(animal_gene_info),
            )

        generate_gmt_files(species, str(outdir))
        build_public_go_custom(database_root, species, annotation_taxid, outdir)

    audit = validate_snapshot(database_root, matrix)
    (database_root / "DATABASE_AUDIT.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    manifest = source_manifest(database_root, source_basic, versions)
    (database_root / "SOURCE_MANIFEST.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return audit


def repair_animaltfdb(matrix_path: Path, database_root: Path, source_basic: Path) -> dict:
    matrix = yaml.safe_load(matrix_path.read_text(encoding="utf-8"))
    tf_basic = database_root / "basic"
    animal_species = {
        spec["latin_name"].replace(" ", "_")
        for spec in matrix["datasets"].values()
        if "AnimalTFDB" in spec["databases"] and spec["species"] != "hsa"
    }
    fetch_tf_sources(tf_basic, offline=True, animal_species=animal_species)
    go_version = latest_version(source_basic / "go", "GO")
    specs = [
        spec for spec in matrix["datasets"].values()
        if "AnimalTFDB" in spec["databases"] and spec["species"] != "hsa"
    ]
    gene_info = prepare_animaltfdb_gene_info(
        source_basic / "go" / go_version / "gene_info.gz",
        tf_basic / "gene_info.animaltfdb_selected_taxa.gz",
        {9606, *(int(spec.get("annotation_taxid", spec["taxid"])) for spec in specs)},
    )
    AnimalTFDBParser.load_external_id_symbol_maps(
        str(gene_info),
        {9606, *(int(spec.get("annotation_taxid", spec["taxid"])) for spec in specs)},
    )
    for spec in specs:
        builder = DatabaseBuilder(str(database_root))
        builder.basic_dir = tf_basic
        outdir = builder.build_animaltfdb(
            spec["species"], int(spec.get("annotation_taxid", spec["taxid"])),
            species_latin=spec["latin_name"].replace(" ", "_"),
            gene_info_path=str(gene_info),
        )
        builder.generate_gmt_files(spec["species"], outdir)
    audit = validate_snapshot(database_root, matrix)
    versions = {
        "go": go_version,
        "reactome": latest_version(source_basic / "reactome", "Reactome"),
        "wikipathways": latest_version(source_basic / "wikipathways", "WP").removeprefix("WP"),
        "disgenet": "AllEnricher-v1-v20190612-free-snapshot",
        "trrust": "v2",
        "chea3": "v2024-frozen-local",
        "animaltfdb": "v4.0",
        "htftarget": "2020",
    }
    (database_root / "DATABASE_AUDIT.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    (database_root / "SOURCE_MANIFEST.json").write_text(
        json.dumps(source_manifest(database_root, source_basic, versions), ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return audit


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--database-root", type=Path, default=DEFAULT_DATABASE_ROOT)
    parser.add_argument("--source-basic", type=Path, default=DEFAULT_SOURCE_BASIC)
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--repair-animaltfdb", action="store_true")
    args = parser.parse_args()
    if args.repair_animaltfdb:
        audit = repair_animaltfdb(
            args.matrix.resolve(), args.database_root.resolve(), args.source_basic.resolve()
        )
    else:
        audit = prepare(
            args.matrix.resolve(), args.database_root.resolve(), args.source_basic.resolve(), args.offline
        )
    print(json.dumps(audit, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
