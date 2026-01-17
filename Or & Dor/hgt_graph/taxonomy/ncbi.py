from Bio import Entrez
from typing import Optional, Dict
from ..constants import ENTREZ_EMAIL, ENTREZ_TOOL, TAXONOMY_DB

Entrez.email = ENTREZ_EMAIL
Entrez.tool = ENTREZ_TOOL


def species_to_tax_id(species_name: str) -> Optional[str]:
    """
    Given a species name, return its taxonomic ID using NCBI Entrez.
    :param species_name: Name of the species to look up.
    :return: Taxonomic ID as a string, or None if not found.
    """
    try:
        handle = Entrez.esearch(
            db=TAXONOMY_DB,
            term=species_name,
            retmode='xml'
        )
        record = Entrez.read(handle)
        handle.close()
        return record['IdList'][0] if record['IdList'] else None
    except Exception:
        return None


def tax_id_to_lineage(tax_id: str) -> Dict[str, str]:
    """
    Return the taxonomic lineage for a given taxonomic ID.
    :param tax_id: Taxonomic ID to look up.
    :return: A dictionary mapping taxonomic ranks to names.
    """
    try:
        handle = Entrez.efetch(
            db='taxonomy',
            id=tax_id,
            retmode='xml'
        )
        record = Entrez.read(handle)
        handle.close()
    except Exception:
        return {}

    lineage = {}
    for item in record[0]['LineageEx']:
        rank = item['Rank']
        name = item['ScientificName']
        lineage[rank] = name

    lineage['species'] = record[0]['ScientificName']
    return lineage