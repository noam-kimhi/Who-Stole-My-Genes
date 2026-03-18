#!/usr/bin/env python3
"""
refseq_fetch_proteins.py

End-to-end pipeline:
1) Read assembly_summary_refseq.txt (streaming)
2) Normalize organism_name -> species_binomial
3) Filter assemblies by:
   - species in a provided species list
   - (optional) version_status == latest
   - assembly_level in allowed set
4) For each species, pick up to K assemblies, prioritizing:
   - refseq_category: reference genome > representative genome > others
   - deterministic tie-breakers
5) Write manifest.tsv (all needed columns + constructed protein_faa_url)
6) Download *_protein.faa.gz for each chosen assembly
7) Parse FASTAs into proteins.tsv

Usage example:
  python refseq_fetch_proteins.py \
    --assembly_summary assembly_summary_refseq.txt \
    --species_list config/species.txt \
    --out_dir out_refseq \
    --max_assemblies_per_species 2 \
    --require_latest \
    --download

Notes:
- This script does NOT use /refseq/release/bacteria shards at all.
- It downloads assembly-specific protein FASTAs from ftp_path.
"""

import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Set, Tuple
import urllib.request

try:
    from graph_construction.fasta_parsing import iter_fasta_gz, protein_id_from_header
except ImportError:
    from fasta_parsing import iter_fasta_gz, protein_id_from_header


def assembly_dirname_from_ftp_path(ftp_path: str) -> str:
    return ftp_path.rstrip("/").split("/")[-1]


# ---------------- Species normalization ----------------

# Matches:
#   "Escherichia coli BL21(DE3)" -> ("Escherichia", "coli")
#   "Bacillus subtilis subsp. subtilis str. 168" -> ("Bacillus", "subtilis")
#   "Candidatus Pelagibacter ubique HTCC1062" -> ("Pelagibacter", "ubique")  [drops "Candidatus"]
SPECIES_RE = re.compile(r'^(?:Candidatus\s+)?([A-Z][a-z]+)\s+([a-z][a-z0-9_-]*)')


def normalize_species(organism_name: str) -> Optional[str]:
    organism_name = organism_name.strip()
    m = SPECIES_RE.match(organism_name)
    if not m:
        return None
    genus, species = m.group(1), m.group(2)
    if species == "sp.":  # drop ambiguous entries like "Bacillus sp. XYZ"
        return None
    return f"{genus} {species}"


# ---------------- Assembly row representation ----------------

@dataclass(frozen=True)
class AssemblyRow:
    assembly_accession: str
    organism_name: str
    species_binomial: str
    taxid: str
    species_taxid: str
    assembly_level: str
    version_status: str
    refseq_category: str
    ftp_path: str  # converted to https://


# RefSeq category priority (lower is better)
REFSEQ_PRIORITY = {
    "reference genome": 0,
    "representative genome": 1,
    # everything else -> 2
}


# ---------------- Read species list ----------------

def read_species_list(path: Path) -> Set[str]:
    """
    species_list file format:
      - one species per line, e.g. "Escherichia coli"
      - blank lines ignored
      - lines starting with # ignored
    """
    s = set()
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            s.add(line)
    return s


# ---------------- Stream assembly_summary_refseq.txt ----------------

def iter_assembly_summary_rows(assembly_summary_path: Path) -> Iterator[Dict[str, str]]:
    """
    Streams rows as dictionaries from assembly_summary_refseq.txt (tab-delimited).
    Skips comments and uses the header line that begins with "# assembly_accession".
    """
    with open(assembly_summary_path, "r", encoding="utf-8", newline="") as f:
        header = None
        for line in f:
            if line.startswith("#assembly_accession"):
                header = line.lstrip("# ").rstrip("\n").split("\t")
                break

        if header is None:
            raise RuntimeError("Header not found. Expected a line starting with '# assembly_accession'.")

        reader = csv.DictReader(f, fieldnames=header, delimiter="\t")
        for row in reader:
            yield row


# ---------------- Filter + select assemblies ----------------

def select_assemblies(
    assembly_summary_path: Path,
    wanted_species: Set[str],
    allowed_levels: Set[str],
    require_latest: bool,
    max_assemblies_per_species: int,
    prefer_refseq_categories: Optional[Set[str]] = None,
) -> List[AssemblyRow]:
    """
    Returns selected AssemblyRow records.

    Selection rule:
      - filter rows by species, allowed_levels, (optional) latest, ftp_path not 'na'
      - group by species_binomial
      - within each species, sort by:
          refseq_category priority (reference > representative > others),
          then assembly_level priority (Complete Genome > Chromosome > others),
          then deterministic tie-breakers (assembly_accession)
      - take up to K rows per species
    """
    by_species: Dict[str, List[AssemblyRow]] = defaultdict(list)

    # Assembly level priority (lower is better)
    level_priority = {
        "Complete Genome": 0,
        "Chromosome": 1,
        "Scaffold": 2,
        "Contig": 3,
    }

    for row in iter_assembly_summary_rows(assembly_summary_path):
        organism_name = (row.get("organism_name") or "").strip()
        ftp_path = (row.get("ftp_path") or "").strip()
        if not organism_name or not ftp_path or ftp_path == "na":
            continue

        species_binomial = normalize_species(organism_name)
        if not species_binomial:
            continue
        if wanted_species and species_binomial not in wanted_species:
            continue

        assembly_level = (row.get("assembly_level") or "").strip()
        if allowed_levels and assembly_level not in allowed_levels:
            continue

        version_status = (row.get("version_status") or "").strip()
        if require_latest and version_status != "latest":
            continue

        refseq_category = (row.get("refseq_category") or "").strip()
        if prefer_refseq_categories is not None and refseq_category not in prefer_refseq_categories:
            continue

        assembly_accession = (row.get("assembly_accession") or "").strip()
        taxid = (row.get("taxid") or "").strip()
        species_taxid = (row.get("species_taxid") or "").strip()

        by_species[species_binomial].append(
            AssemblyRow(
                assembly_accession=assembly_accession,
                organism_name=organism_name,
                species_binomial=species_binomial,
                taxid=taxid,
                species_taxid=species_taxid,
                assembly_level=assembly_level,
                version_status=version_status,
                refseq_category=refseq_category,
                ftp_path=ftp_path.replace("ftp://", "https://"),
            )
        )

    selected: List[AssemblyRow] = []

    for species, rows in by_species.items():
        rows.sort(
            key=lambda r: (
                REFSEQ_PRIORITY.get(r.refseq_category, 2),
                level_priority.get(r.assembly_level, 9),
                r.assembly_accession,  # deterministic tie-break
            )
        )
        selected.extend(rows[:max_assemblies_per_species])

    # Deterministic ordering of the full output too
    selected.sort(key=lambda r: (r.species_binomial, r.assembly_accession))
    return selected


# ---------------- Manifest writing ----------------

def write_manifest(rows: List[AssemblyRow], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "species_binomial",
            "organism_name",
            "assembly_accession",
            "taxid",
            "species_taxid",
            "assembly_level",
            "version_status",
            "refseq_category",
            "ftp_path",
            "protein_faa_url",
        ])
        for r in rows:
            dirname = assembly_dirname_from_ftp_path(r.ftp_path)
            url = f"{r.ftp_path}/{dirname}_protein.faa.gz"
            w.writerow([
                r.species_binomial,
                r.organism_name,
                r.assembly_accession,
                r.taxid,
                r.species_taxid,
                r.assembly_level,
                r.version_status,
                r.refseq_category,
                r.ftp_path,
                url,
            ])


# ---------------- Downloading ----------------

def download_file(url: str, out_path: Path, retries: int = 2) -> bool:
    """
    Download url -> out_path. Returns True on success.
    Skips download if file already exists and is non-empty.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists() and out_path.stat().st_size > 0:
        return True

    last_err = None
    for _ in range(retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                # Some urllib responses do not expose .status; but NCBI will raise on HTTPError for non-200
                with open(out_path, "wb") as f:
                    while True:
                        chunk = resp.read(1 << 20)
                        if not chunk:
                            break
                        f.write(chunk)
            return True
        except Exception as e:
            last_err = e

    print(f"[WARN] Failed download: {url} ({last_err})", file=sys.stderr)
    return False


def write_proteins_tsv(
    proteins_out: Path,
    rows: List[AssemblyRow],
    downloads_dir: Path,
    min_len: int,
    include_sequence: bool,
) -> None:
    """
    Write a table with proteins from all downloaded FASTAs.

    Columns always include:
      protein_id, assembly_accession, species_binomial, organism_name, taxid, species_taxid, length
    Optionally include: sequence
    """
    proteins_out.parent.mkdir(parents=True, exist_ok=True)

    with open(proteins_out, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        cols = [
            "protein_id",
            "assembly_accession",
            "species_binomial",
            "organism_name",
            "taxid",
            "species_taxid",
            "length",
        ]
        if include_sequence:
            cols.append("sequence")
        w.writerow(cols)

        for r in rows:
            dirname = assembly_dirname_from_ftp_path(r.ftp_path)
            faa_path = downloads_dir / f"{dirname}_protein.faa.gz"
            if not faa_path.exists():
                continue

            for header, seq in iter_fasta_gz(faa_path):
                if len(seq) < min_len:
                    continue
                pid = protein_id_from_header(header)
                out_row = [
                    pid,
                    r.assembly_accession,
                    r.species_binomial,
                    r.organism_name,
                    r.taxid,
                    r.species_taxid,
                    str(len(seq)),
                ]
                if include_sequence:
                    out_row.append(seq)
                w.writerow(out_row)


# ---------------- CLI ----------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Filter assembly_summary_refseq.txt, select K assemblies/species, download *_protein.faa.gz, parse proteins."
    )
    ap.add_argument("--assembly_summary", required=True, type=Path,
                    help="Path to assembly_summary_refseq.txt")
    ap.add_argument("--species_list", required=True, type=Path,
                    help="Text file: one species per line (Genus species)")
    ap.add_argument("--out_dir", required=True, type=Path,
                    help="Output directory for manifest, downloads, proteins.tsv")
    ap.add_argument("--max_assemblies_per_species", type=int, default=2,
                    help="Cap number of assemblies per species (default: 2)")
    ap.add_argument("--assembly_levels", type=str, default="Complete Genome,Chromosome",
                    help="Comma-separated allowed assembly levels (default: Complete Genome,Chromosome)")
    ap.add_argument("--require_latest", action="store_true",
                    help="Require version_status == latest")
    ap.add_argument("--prefer_refseq_categories", type=str, default="",
                    help="Optional hard filter: comma-separated refseq_category values. "
                         "Example: 'reference genome,representative genome'")
    ap.add_argument("--min_protein_len", type=int, default=50,
                    help="Minimum protein length to keep (default: 50 aa)")
    ap.add_argument("--include_sequence", action="store_true",
                    help="Include full sequences in proteins.tsv (will be huge)")
    ap.add_argument("--download", action="store_true",
                    help="Actually download FASTAs and parse them (otherwise only write manifest)")
    args = ap.parse_args()

    wanted_species = read_species_list(args.species_list)
    allowed_levels = {s.strip() for s in args.assembly_levels.split(",") if s.strip()}

    prefer_cats = None
    if args.prefer_refseq_categories.strip():
        prefer_cats = {s.strip() for s in args.prefer_refseq_categories.split(",") if s.strip()}

    out_dir = args.out_dir
    downloads_dir = out_dir / "downloads"
    manifest_path = out_dir / "manifest.tsv"
    proteins_path = out_dir / "proteins.tsv"

    rows = select_assemblies(
        assembly_summary_path=args.assembly_summary,
        wanted_species=wanted_species,
        allowed_levels=allowed_levels,
        require_latest=args.require_latest,
        max_assemblies_per_species=args.max_assemblies_per_species,
        prefer_refseq_categories=prefer_cats,
    )

    if not rows:
        print("[ERROR] No assemblies matched your filters/species list.", file=sys.stderr)
        sys.exit(2)

    write_manifest(rows, manifest_path)
    print(f"[OK] Wrote manifest: {manifest_path} ({len(rows)} assemblies)")

    if not args.download:
        print("[INFO] Skipping downloads (run again with --download).")
        return

    ok = 0
    for r in rows:
        dirname = assembly_dirname_from_ftp_path(r.ftp_path)
        url = f"{r.ftp_path}/{dirname}_protein.faa.gz"
        out_faa = downloads_dir / f"{dirname}_protein.faa.gz"
        if download_file(url, out_faa, retries=2):
            ok += 1
    print(f"[OK] Downloaded: {ok}/{len(rows)} protein FASTAs to {downloads_dir}")

    write_proteins_tsv(
        proteins_out=proteins_path,
        rows=rows,
        downloads_dir=downloads_dir,
        min_len=args.min_protein_len,
        include_sequence=args.include_sequence,
    )
    print(f"[OK] Wrote proteins table: {proteins_path}")


if __name__ == "__main__":
    main()
