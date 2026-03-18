import csv
import gzip
from pathlib import Path
from typing import Dict, Iterator, List, Tuple


def iter_fasta_gz(path: Path) -> Iterator[Tuple[str, str]]:
    """Yield (header_without_>, sequence) from a gzipped FASTA."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        header = None
        seq_parts: List[str] = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            yield header, "".join(seq_parts)


def protein_id_from_header(header: str) -> str:
    """First token of header is usually accession (NP_/WP_/...)."""
    return header.split()[0]


# ---------------- Species mapping from manifest ----------------

def basename_no_slash(p: str) -> str:
    return p.rstrip("/").split("/")[-1]


def load_manifest_species_map(manifest_path: Path) -> Dict[str, str]:
    """
    Return map: assembly_dirname -> species_binomial

    In the fixed pipeline, downloaded FASTAs are named:
      {assembly_dirname}_protein.faa.gz
    where assembly_dirname is basename of ftp_path, e.g.:
      ftp_path .../GCF_000005845.2_ASM584v2
      dirname = GCF_000005845.2_ASM584v2
      file = GCF_000005845.2_ASM584v2_protein.faa.gz
    """
    m: Dict[str, str] = {}
    with open(manifest_path, "r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            ftp_path = row["ftp_path"]
            species = row["species_binomial"]
            dirname = basename_no_slash(ftp_path)
            m[dirname] = species
    return m

