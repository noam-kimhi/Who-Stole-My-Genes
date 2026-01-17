import time
import json
from pathlib import Path
from typing import Dict, Union, Optional, Mapping, Any
from ..constants import ORGANISM_KEY, NCBI_SLEEP, PathLike
from .normalize import normalize_species_name
from .ncbi import species_to_tax_id, tax_id_to_lineage


def build_taxonomy_map(node_attrs: Dict[str, Dict[str, Union[str, int]]],
                       save_path: str, sleep_seconds: float) -> Dict[str, Mapping[str, Any]]:
    """
    Build a taxonomy map from species names to their taxonomic lineage.
    :param node_attrs: A dictionary of node attributes.
    :param save_path: Path to save the taxonomy cache JSON file.
    :param sleep_seconds: Seconds to sleep between NCBI queries.
    :return: A dictionary mapping species names to their taxonomic lineage.
    """
    unique_species = set()
    for attrs in node_attrs.values():
        org = attrs.get(ORGANISM_KEY, '')
        unique_species.add(normalize_species_name(str(org)))

    taxonomy_cache: Dict[str, Mapping[str, Any]] = {}

    for sp in sorted(unique_species):
        if not sp:
            continue

        tax_id = species_to_tax_id(sp)
        if tax_id is None:
            taxonomy_cache[sp] = {}  # keep empty to avoid re-query
            continue

        lineage = tax_id_to_lineage(tax_id)
        taxonomy_cache[sp] = lineage

        time.sleep(sleep_seconds)

    # Save the taxonomy cache to a JSON file
    with open(save_path, 'w') as f:
        json.dump(taxonomy_cache, f, indent=4)

    return taxonomy_cache


def load_taxonomy_cache(path: PathLike) -> Dict[str, Mapping[str, Any]]:
    """
    Load taxonomy cache JSON from disk.
    :param path: Path to JSON cache file.
    :return: Taxonomy cache dict.
    :raises FileNotFoundError: if the file does not exist.
    :raises json.JSONDecodeError: if the file is not valid JSON.
    """
    p = Path(path)
    with p.open('r') as f:
        return json.load(f)


def get_or_build_taxonomy_cache(cache_path: Optional[PathLike], *, node_attrs: Dict[str, Dict[str, Union[str, int]]],
                                save_path: PathLike, sleep_seconds: float = NCBI_SLEEP) -> Dict[str, Mapping[str, Any]]:
    """
    Get taxonomy cache from disk or build it if not available.
    If cache_path is provided, load the cache from that path.
    Otherwise, build the taxonomy map from node attributes and save it to save_path.
    :param cache_path: The path to the cache file (if any).
    :param node_attrs: The node attributes dictionary.
    :param save_path: The path to save the cache file if built.
    :param sleep_seconds: Seconds to sleep between NCBI queries.
    :return: The taxonomy cache dictionary.
    """
    if cache_path:
        return load_taxonomy_cache(cache_path)

    return build_taxonomy_map(node_attrs, save_path=str(save_path), sleep_seconds=sleep_seconds)
