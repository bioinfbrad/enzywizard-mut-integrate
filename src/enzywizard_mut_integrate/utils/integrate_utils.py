from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Tuple
import json
import re

from ..utils.logging_utils import Logger
from ..utils.IO_utils import write_json_from_dict_inline_leaf_lists
from ..utils.common_utils import get_optimized_filename
from ..resources.aa_resources import AA_20NAME_INDEX, AA_8CLASSES, DSSP_8STATE_INDEX
from ..utils.sequence_utils import normalize_aa_name_to_one_letter


SUPPORTED_OUTPUT_TYPES = {
    "enzywizard_clean",
    "enzywizard_aaprops",
    "enzywizard_hydrocluster",
    "enzywizard_energy",
    "enzywizard_flexibility",
    "enzywizard_disorder",
    "enzywizard_conservation",
    "enzywizard_embedding",
    "enzywizard_pocket",
    "enzywizard_substrate",
    "enzywizard_dock",
    "enzywizard_interaction",
}

INTERACTION_ORDER = ["HBOND", "IONIC", "VDW", "PIPISTACK", "PICATION", "SSBOND"]
NODE_TYPE_ORDER = ["residue", "substrate"]
AA_CLASS_ORDER = [item[0] for item in AA_8CLASSES]
AA_NAME_ORDER = [aa for aa, _ in sorted(AA_20NAME_INDEX.items(), key=lambda x: x[1])]
AA_SS_ORDER = [ss for ss, _ in sorted(DSSP_8STATE_INDEX.items(), key=lambda x: x[1])]
AA_NAME_COUNT_FIELD_ORDER = [
    "alanine_count",
    "cysteine_count",
    "aspartic_acid_count",
    "glutamic_acid_count",
    "phenylalanine_count",
    "glycine_count",
    "histidine_count",
    "isoleucine_count",
    "lysine_count",
    "leucine_count",
    "methionine_count",
    "asparagine_count",
    "proline_count",
    "glutamine_count",
    "arginine_count",
    "serine_count",
    "threonine_count",
    "valine_count",
    "tryptophan_count",
    "tyrosine_count",
]

AA_CLASS_COUNT_FIELD_ORDER = [
    "uncharged_polar_count",
    "positively_charged_count",
    "negatively_charged_count",
    "hydrophobic_count",
    "aromatic_count",
    "aliphatic_count",
    "heterocyclic_count",
    "sulfur_containing_count",
]

AA_SS_COUNT_FIELD_ORDER = [
    "unassigned_count",
    "alpha_helix_count",
    "beta_bridge_count",
    "extended_strand_count",
    "three_ten_helix_count",
    "pi_helix_count",
    "turn_count",
    "bend_count",
]


def load_json_file(path: str | Path, logger: Logger) -> Dict[str, Any] | None:
    try:
        p = Path(path)
        if not p.exists() or not p.is_file():
            logger.print(f"[ERROR] JSON file not found: {p}")
            return None

        with p.open("r", encoding="utf-8") as f:
            data = json.load(f)

        if not isinstance(data, dict):
            logger.print(f"[ERROR] JSON top-level object must be a dict: {p}")
            return None

        return data

    except Exception as e:
        logger.print(f"[ERROR] Failed to load JSON file {path}: {e}")
        return None


def list_json_files(input_dir: str | Path, logger: Logger) -> List[Path] | None:
    try:
        p = Path(input_dir)
        if not p.exists() or not p.is_dir():
            logger.print(f"[ERROR] Invalid input_dir: {p}")
            return None

        json_files = sorted([x for x in p.iterdir() if x.is_file() and x.suffix.lower() == ".json"])
        return json_files

    except Exception as e:
        logger.print(f"[ERROR] Failed to list JSON files from {input_dir}: {e}")
        return None

def find_unique_clean_report_path(json_path_list: List[Path], logger: Logger) -> Path | None:
    if not isinstance(json_path_list, list):
        logger.print("[ERROR] json_path_list must be a list.")
        return None

    clean_report_path_list: List[Path] = []

    for json_path in json_path_list:
        data = load_json_file(json_path, logger)
        if data is None:
            return None

        report_type = get_supported_output_type(data, logger)
        if report_type is None:
            return None

        if report_type == "enzywizard_clean":
            clean_report_path_list.append(json_path)

    if len(clean_report_path_list) == 0:
        logger.print("[ERROR] No enzywizard_clean report found in input_dir.")
        return None

    if len(clean_report_path_list) > 1:
        logger.print(f"[ERROR] Multiple enzywizard_clean reports found in input_dir: {len(clean_report_path_list)}")
        return None

    return clean_report_path_list[0]


def extract_protein_name_from_clean_report_path(clean_report_path: str | Path, logger: Logger) -> str | None:
    try:
        name = Path(clean_report_path).name
        m = re.fullmatch(r"clean_report_(.+)\.json", name)
        if m is None:
            logger.print(f"[ERROR] clean_report file name must match clean_report_{{protein_name}}.json: {name}")
            return None

        protein_name = m.group(1).strip()
        if protein_name == "":
            logger.print(f"[ERROR] Invalid protein name in clean_report file name: {name}")
            return None

        return protein_name

    except Exception as e:
        logger.print(f"[ERROR] Failed to parse clean_report file name: {e}")
        return None


def save_integrate_json(report: Dict[str, Any] | List[Any], output_path: str | Path, logger: Logger) -> bool:
    try:
        write_json_from_dict_inline_leaf_lists(report, output_path)
        return True
    except Exception as e:
        logger.print(f"[ERROR] Failed to save integrate JSON: {e}")
        return False

def split_integrated_graph_entries(
    integrated_graph: List[Dict[str, Any]],
    logger: Logger,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]] | None:
    if not isinstance(integrated_graph, list):
        logger.print("[ERROR] integrated_graph must be a list.")
        return None

    node_by_index: Dict[int, Dict[str, Any]] = {}
    edge_list: List[Dict[str, Any]] = []

    def add_node(node: Any) -> bool:
        if not isinstance(node, dict):
            logger.print("[ERROR] Invalid node in integrated_graph.")
            return False

        node_index = node.get("node_index")
        if not isinstance(node_index, int):
            logger.print("[ERROR] Invalid or missing node_index in integrated_graph node.")
            return False

        if node_index in node_by_index:
            return True

        node_by_index[node_index] = node
        return True

    for item in integrated_graph:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid integrated_graph item.")
            return None

        if "isolated_node" in item:
            if not add_node(item["isolated_node"]):
                return None

        elif (
            "molecular_interaction" in item
            and "source_node" in item
            and "target_node" in item
        ):
            source_node = item["source_node"]
            target_node = item["target_node"]

            if not add_node(source_node):
                return None
            if not add_node(target_node):
                return None

            source_node_index = source_node.get("node_index")
            target_node_index = target_node.get("node_index")

            if not isinstance(source_node_index, int) or not isinstance(target_node_index, int):
                logger.print("[ERROR] Invalid node_index in edge source_node or target_node.")
                return None

            edge_list.append({
                "molecular_interaction": item["molecular_interaction"],
                "source_node": {
                    "node_index": source_node_index,
                },
                "target_node": {
                    "node_index": target_node_index,
                },
            })

        else:
            logger.print("[ERROR] Invalid integrated_graph entry structure.")
            return None

    node_list = [
        node_by_index[node_index]
        for node_index in sorted(node_by_index.keys())
    ]

    return node_list, edge_list



def interaction_one_hot(interaction: str) -> List[int]:
    vec = [0] * len(INTERACTION_ORDER)
    if interaction in INTERACTION_ORDER:
        vec[INTERACTION_ORDER.index(interaction)] = 1
    return vec


def node_type_one_hot(node_type: str) -> List[int]:
    vec = [0] * len(NODE_TYPE_ORDER)
    if node_type in NODE_TYPE_ORDER:
        vec[NODE_TYPE_ORDER.index(node_type)] = 1
    return vec


def aa_name_statistics_to_list(statistics: Dict[str, int], logger: Logger) -> List[int] | None:
    if not isinstance(statistics, dict):
        logger.print("[ERROR] residue_name_statistics must be a dict.")
        return None

    out: List[int] = []
    for field in AA_NAME_COUNT_FIELD_ORDER:
        value = statistics.get(field)
        if not isinstance(value, int):
            logger.print(f"[ERROR] Missing or invalid residue_name_statistics field: {field}.")
            return None
        out.append(value)

    return out


def aa_class_statistics_to_list(statistics: Dict[str, int], logger: Logger) -> List[int] | None:
    if not isinstance(statistics, dict):
        logger.print("[ERROR] residue_chemical_classification_statistics must be a dict.")
        return None

    out: List[int] = []
    for field in AA_CLASS_COUNT_FIELD_ORDER:
        value = statistics.get(field)
        if not isinstance(value, int):
            logger.print(
                f"[ERROR] Missing or invalid residue_chemical_classification_statistics field: {field}."
            )
            return None
        out.append(value)

    return out


def aa_ss_statistics_to_list(statistics: Dict[str, int], logger: Logger) -> List[int] | None:
    if not isinstance(statistics, dict):
        logger.print("[ERROR] residue_secondary_structure_statistics must be a dict.")
        return None

    out: List[int] = []
    for field in AA_SS_COUNT_FIELD_ORDER:
        value = statistics.get(field)
        if not isinstance(value, int):
            logger.print(
                f"[ERROR] Missing or invalid residue_secondary_structure_statistics field: {field}."
            )
            return None
        out.append(value)

    return out



def residue_key(aa_id: int, aa_name: str) -> Tuple[int, str]:
    aa_name_norm = normalize_aa_name_to_one_letter(aa_name)
    return int(aa_id), aa_name_norm


def residue_key_from_dict(data: Dict[str, Any], id_key: str, name_key: str) -> Tuple[int, str] | None:
    aa_id = data.get(id_key)
    aa_name = data.get(name_key)

    if not isinstance(aa_id, int):
        return None
    if not isinstance(aa_name, str) or aa_name.strip() == "":
        return None

    return residue_key(aa_id, aa_name)


def substrate_key(substrate_index: int, substrate_name: str) -> Tuple[int, str]:
    return int(substrate_index), str(substrate_name).strip()


def build_residue_membership_set(residue_list: List[Dict[str, Any]], logger: Logger) -> set[Tuple[int, str]] | None:
    if not isinstance(residue_list, list):
        logger.print("[ERROR] residue_list must be a list.")
        return None

    out: set[Tuple[int, str]] = set()

    for item in residue_list:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid residue item.")
            return None

        key = residue_key_from_dict(item, "residue_index", "residue_name")
        if key is None:
            logger.print("[ERROR] Invalid residue key in residue list.")
            return None
        out.add(key)

    return out


def build_lookup_by_residue(
    item_list: List[Dict[str, Any]],
    id_key: str,
    name_key: str,
    logger: Logger
) -> Dict[Tuple[int, str], Dict[str, Any]] | None:
    if not isinstance(item_list, list):
        logger.print("[ERROR] item_list must be a list.")
        return None

    lookup: Dict[Tuple[int, str], Dict[str, Any]] = {}

    for item in item_list:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid item in lookup list.")
            return None

        key = residue_key_from_dict(item, id_key, name_key)
        if key is None:
            logger.print("[ERROR] Invalid residue identifier in lookup list.")
            return None

        if key in lookup:
            logger.print(f"[ERROR] Duplicate residue entry found: {key}")
            return None

        lookup[key] = item

    return lookup


def get_hydrophobic_cluster_membership_set(hydro_report: Dict[str, Any], logger: Logger) -> set[Tuple[int, str]] | None:
    clusters = hydro_report.get("hydrophobic_clusters")
    if not isinstance(clusters, list):
        logger.print("[ERROR] Invalid hydrophobic_clusters.")
        return None

    out: set[Tuple[int, str]] = set()

    for cluster in clusters:
        if not isinstance(cluster, dict):
            logger.print("[ERROR] Invalid hydrophobic cluster entry.")
            return None
        residues = cluster.get("residues")
        membership = build_residue_membership_set(residues, logger)
        if membership is None:
            return None
        out |= membership

    return out


def get_disorder_membership_set(disorder_report: Dict[str, Any], logger: Logger) -> set[Tuple[int, str]] | None:
    regions = disorder_report.get("disordered_regions")
    if not isinstance(regions, list):
        logger.print("[ERROR] Invalid disordered_regions.")
        return None

    out: set[Tuple[int, str]] = set()

    for region in regions:
        if not isinstance(region, dict):
            logger.print("[ERROR] Invalid disordered region entry.")
            return None
        residues = region.get("residues")
        membership = build_residue_membership_set(residues, logger)
        if membership is None:
            return None
        out |= membership

    return out


def get_pocket_membership_set(pocket_report: Dict[str, Any], logger: Logger) -> set[Tuple[int, str]] | None:
    regions = pocket_report.get("binding_pockets")
    if not isinstance(regions, list):
        logger.print("[ERROR] Invalid binding_pockets.")
        return None

    out: set[Tuple[int, str]] = set()

    for region in regions:
        if not isinstance(region, dict):
            logger.print("[ERROR] Invalid binding pocket entry.")
            return None
        residues = region.get("residues")
        membership = build_residue_membership_set(residues, logger)
        if membership is None:
            return None
        out |= membership

    return out


def get_clean_new_residue_list(clean_report: Dict[str, Any], logger: Logger) -> List[Dict[str, Any]] | None:
    mapping = clean_report.get("residue_mapping_old_to_new")
    if not isinstance(mapping, list):
        logger.print("[ERROR] Invalid residue_mapping_old_to_new in clean report.")
        return None

    if len(mapping) == 0:
        logger.print("[ERROR] residue_mapping_old_to_new is empty.")
        return None

    new_residue_list: List[Dict[str, Any]] = []

    for item in mapping:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid mapping entry in clean report.")
            return None

        new_residue = item.get("new_residue")
        if not isinstance(new_residue, dict):
            logger.print("[ERROR] Missing new_residue in clean report.")
            return None

        residue_index = new_residue.get("residue_index")
        residue_name = new_residue.get("residue_name")

        if not isinstance(residue_index, int):
            logger.print("[ERROR] Invalid new_residue residue_index in clean report.")
            return None
        if not isinstance(residue_name, str) or residue_name.strip() == "":
            logger.print("[ERROR] Invalid new_residue residue_name in clean report.")
            return None

        new_residue_list.append({
            "residue_index": residue_index,
            "residue_name": residue_name.strip().upper(),
        })

    new_residue_list.sort(key=lambda x: x["residue_index"])
    return new_residue_list


def get_supported_output_type(data: Dict[str, Any], logger: Logger) -> str | None:
    report_type = data.get("report_type")
    if not isinstance(report_type, str):
        logger.print("[ERROR] Missing or invalid report_type.")
        return None
    if report_type not in SUPPORTED_OUTPUT_TYPES:
        logger.print(f"[ERROR] Unsupported report_type for integrate: {report_type}")
        return None
    return report_type


def validate_report_by_type(data: Dict[str, Any], logger: Logger) -> bool:
    report_type = get_supported_output_type(data, logger)
    if report_type is None:
        return False

    if report_type == "enzywizard_clean":
        return validate_clean_report(data, logger)
    if report_type == "enzywizard_aaprops":
        return validate_aaprops_report(data, logger)
    if report_type == "enzywizard_hydrocluster":
        return validate_hydrocluster_report(data, logger)
    if report_type == "enzywizard_energy":
        return validate_energy_report(data, logger)
    if report_type == "enzywizard_flexibility":
        return validate_flexibility_report(data, logger)
    if report_type == "enzywizard_disorder":
        return validate_disorder_report(data, logger)
    if report_type == "enzywizard_conservation":
        return validate_conservation_report(data, logger)
    if report_type == "enzywizard_embedding":
        return validate_embedding_report(data, logger)
    if report_type == "enzywizard_pocket":
        return validate_pocket_report(data, logger)
    if report_type == "enzywizard_substrate":
        return validate_substrate_report(data, logger)
    if report_type == "enzywizard_dock":
        return validate_dock_report(data, logger)
    if report_type == "enzywizard_interaction":
        return validate_interaction_report(data, logger)

    logger.print("[ERROR] Unsupported report type.")
    return False


def _is_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def _is_number_list(x: Any, min_len: int | None = None) -> bool:
    if not isinstance(x, list):
        return False
    if min_len is not None and len(x) < min_len:
        return False
    return all(_is_number(v) for v in x)


def _is_int01_list(x: Any) -> bool:
    if not isinstance(x, list):
        return False
    return all(isinstance(v, int) and v in {0, 1} for v in x)


def _validate_residue_index_name(item: Dict[str, Any], logger: Logger) -> bool:
    if not isinstance(item.get("residue_index"), int):
        logger.print("[ERROR] Invalid residue_index.")
        return False
    if not isinstance(item.get("residue_name"), str) or item["residue_name"].strip() == "":
        logger.print("[ERROR] Invalid residue_name.")
        return False
    return True


def validate_clean_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("report_type") != "enzywizard_clean":
        logger.print("[ERROR] clean report report_type mismatch.")
        return False

    mapping = data.get("residue_mapping_old_to_new")
    stats = data.get("clean_statistics")

    if not isinstance(mapping, list):
        logger.print("[ERROR] Invalid residue_mapping_old_to_new.")
        return False
    if not isinstance(stats, dict):
        logger.print("[ERROR] Invalid clean_statistics.")
        return False

    required_stat_keys = [
        "removed_heterogen_count",
        "standardized_residue_name_count",
        "repaired_residue_count",
        "added_heavy_atom_count",
        "added_hydrogen_atom_count",
        "retained_residue_count",
    ]
    for key in required_stat_keys:
        if not isinstance(stats.get(key), int):
            logger.print(f"[ERROR] Invalid clean_statistics field: {key}")
            return False

    for item in mapping:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid residue mapping entry.")
            return False

        old_residue = item.get("old_residue")
        new_residue = item.get("new_residue")

        if not isinstance(old_residue, dict) or not isinstance(new_residue, dict):
            logger.print("[ERROR] Invalid residue mapping structure.")
            return False

        for residue_item in [old_residue, new_residue]:
            if not _validate_residue_index_name(residue_item, logger):
                return False
            if not isinstance(residue_item.get("hydrogen_atom_count"), int):
                logger.print("[ERROR] Invalid hydrogen_atom_count.")
                return False

    return True


def validate_aaprops_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("report_type") != "enzywizard_aaprops":
        logger.print("[ERROR] aaprops report report_type mismatch.")
        return False

    residue_properties = data.get("amino_acid_residue_properties")
    stats = data.get("residue_properties_statistics")

    if not isinstance(residue_properties, list):
        logger.print("[ERROR] Invalid amino_acid_residue_properties.")
        return False
    if not isinstance(stats, dict):
        logger.print("[ERROR] Invalid residue_properties_statistics.")
        return False

    for key in [
        "residue_name_statistics",
        "residue_chemical_classification_statistics",
        "residue_secondary_structure_statistics",
    ]:
        if not isinstance(stats.get(key), dict):
            logger.print(f"[ERROR] Invalid residue_properties_statistics field: {key}")
            return False

    for item in residue_properties:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid amino_acid_residue_properties entry.")
            return False

        if not isinstance(item.get("residue_index"), int):
            logger.print("[ERROR] Invalid aaprops field: residue_index")
            return False

        for key in [
            "residue_name",
            "residue_chemical_classification",
            "residue_secondary_structure",
        ]:
            if not isinstance(item.get(key), str):
                logger.print(f"[ERROR] Invalid aaprops field: {key}")
                return False

        for key in [
            "residue_name_one_hot_encoding",
            "residue_chemical_classification_one_hot_encoding",
            "residue_secondary_structure_one_hot_encoding",
        ]:
            if not _is_int01_list(item.get(key)):
                logger.print(f"[ERROR] Invalid one-hot field: {key}")
                return False

        for key in [
            "residue_relative_solvent_accessibility",
            "residue_backbone_phi_angle",
            "residue_backbone_psi_angle",
            "residue_net_charge",
            "residue_pka",
            "residue_volume",
            "residue_hydrophobicity",
            "residue_molecular_weight",
            "residue_isoelectric_point",
        ]:
            if not _is_number(item.get(key)):
                logger.print(f"[ERROR] Invalid numeric field: {key}")
                return False

        if not _is_number_list(item.get("residue_alpha_carbon_coordinate"), min_len=3):
            logger.print("[ERROR] Invalid residue_alpha_carbon_coordinate.")
            return False

    return True


def validate_hydrocluster_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("report_type") != "enzywizard_hydrocluster":
        logger.print("[ERROR] hydrocluster report report_type mismatch.")
        return False

    stats = data.get("hydrophobic_cluster_statistics")
    clusters = data.get("hydrophobic_clusters")

    if not isinstance(stats, dict):
        logger.print("[ERROR] Invalid hydrophobic_cluster_statistics.")
        return False

    if not isinstance(stats.get("hydrophobic_cluster_count"), int):
        logger.print("[ERROR] Invalid hydrophobic_cluster_statistics.hydrophobic_cluster_count.")
        return False
    if not _is_number(stats.get("max_hydrophobic_cluster_area")):
        logger.print("[ERROR] Invalid hydrophobic_cluster_statistics.max_hydrophobic_cluster_area.")
        return False
    if not _is_number(stats.get("total_hydrophobic_cluster_area")):
        logger.print("[ERROR] Invalid hydrophobic_cluster_statistics.total_hydrophobic_cluster_area.")
        return False

    if not isinstance(clusters, list):
        logger.print("[ERROR] Invalid hydrophobic_clusters.")
        return False

    for item in clusters:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid hydrophobic_clusters entry.")
            return False

        if not _is_number(item.get("hydrophobic_cluster_area")):
            logger.print("[ERROR] Invalid hydrophobic_cluster_area.")
            return False

        residues = item.get("residues")
        if not isinstance(residues, list):
            logger.print("[ERROR] Invalid hydrophobic cluster residues.")
            return False

        for residue in residues:
            if not isinstance(residue, dict) or not _validate_residue_index_name(residue, logger):
                return False

    return True


def validate_energy_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("report_type") != "enzywizard_energy":
        logger.print("[ERROR] energy report report_type mismatch.")
        return False

    terms = data.get("energy_terms")
    if not isinstance(terms, dict):
        logger.print("[ERROR] Invalid energy_terms.")
        return False

    keys = [
        "total_potential_energy",
        "harmonic_bond_potential_energy",
        "harmonic_angle_potential_energy",
        "custom_bond_potential_energy",
        "custom_torsion_potential_energy",
        "custom_nonbonded_potential_energy",
        "nonbonded_potential_energy",
        "periodic_torsion_potential_energy",
        "cmap_torsion_potential_energy",
    ]

    for key in keys:
        if not _is_number(terms.get(key)):
            logger.print(f"[ERROR] Invalid energy term: {key}")
            return False

    return True


def validate_flexibility_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("report_type") != "enzywizard_flexibility":
        logger.print("[ERROR] flexibility report report_type mismatch.")
        return False

    entries = data.get("protein_flexibility")
    if not isinstance(entries, list):
        logger.print("[ERROR] Invalid protein_flexibility.")
        return False

    for item in entries:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid protein_flexibility entry.")
            return False

        if not _validate_residue_index_name(item, logger):
            return False

        if not _is_number(item.get("residue_root_mean_square_fluctuation")):
            logger.print("[ERROR] Invalid residue_root_mean_square_fluctuation.")
            return False

    return True


def validate_disorder_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("report_type") != "enzywizard_disorder":
        logger.print("[ERROR] disorder report report_type mismatch.")
        return False

    stats = data.get("disordered_region_statistics")
    regions = data.get("disordered_regions")

    if not isinstance(stats, dict):
        logger.print("[ERROR] Invalid disordered_region_statistics.")
        return False

    if not isinstance(stats.get("disordered_region_count"), int):
        logger.print("[ERROR] Invalid disordered_region_statistics.disordered_region_count.")
        return False
    if not isinstance(stats.get("max_disordered_region_length"), int):
        logger.print("[ERROR] Invalid disordered_region_statistics.max_disordered_region_length.")
        return False
    if not isinstance(stats.get("total_disordered_region_length"), int):
        logger.print("[ERROR] Invalid disordered_region_statistics.total_disordered_region_length.")
        return False

    if not isinstance(regions, list):
        logger.print("[ERROR] Invalid disordered_regions.")
        return False

    for item in regions:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid disordered region entry.")
            return False

        if not isinstance(item.get("disordered_region_length"), int):
            logger.print("[ERROR] Invalid disordered_region_length.")
            return False

        residues = item.get("residues")
        if not isinstance(residues, list):
            logger.print("[ERROR] Invalid disordered region residues.")
            return False

        for residue in residues:
            if not isinstance(residue, dict) or not _validate_residue_index_name(residue, logger):
                return False

    return True


def validate_conservation_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("report_type") != "enzywizard_conservation":
        logger.print("[ERROR] conservation report report_type mismatch.")
        return False

    entries = data.get("sequence_conservation_scores")
    if not isinstance(entries, list):
        logger.print("[ERROR] Invalid sequence_conservation_scores.")
        return False

    for item in entries:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid sequence_conservation_scores entry.")
            return False

        if not _validate_residue_index_name(item, logger):
            return False

        for key in [
            "hmm_profile_raw_score",
            "normalized_emission_probability",
            "normalized_shannon_information_content",
        ]:
            if not _is_number(item.get(key)):
                logger.print(f"[ERROR] Invalid conservation field: {key}")
                return False

    return True


def validate_embedding_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("report_type") != "enzywizard_embedding":
        logger.print("[ERROR] embedding report report_type mismatch.")
        return False

    entries = data.get("sequence_embeddings")
    if not isinstance(entries, list):
        logger.print("[ERROR] Invalid sequence_embeddings.")
        return False

    for item in entries:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid sequence_embeddings entry.")
            return False

        if not _validate_residue_index_name(item, logger):
            return False

        if not _is_number_list(item.get("residue_embedding"), min_len=1):
            logger.print("[ERROR] Invalid residue_embedding vector.")
            return False

    return True


def validate_pocket_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("report_type") != "enzywizard_pocket":
        logger.print("[ERROR] pocket report report_type mismatch.")
        return False

    stats = data.get("binding_pocket_statistics")
    pockets = data.get("binding_pockets")

    if not isinstance(stats, dict):
        logger.print("[ERROR] Invalid binding_pocket_statistics.")
        return False

    if not isinstance(stats.get("binding_pocket_count"), int):
        logger.print("[ERROR] Invalid binding_pocket_statistics.binding_pocket_count.")
        return False
    if not _is_number(stats.get("max_binding_pocket_volume")):
        logger.print("[ERROR] Invalid binding_pocket_statistics.max_binding_pocket_volume.")
        return False
    if not _is_number(stats.get("total_binding_pocket_volume")):
        logger.print("[ERROR] Invalid binding_pocket_statistics.total_binding_pocket_volume.")
        return False

    if not isinstance(pockets, list):
        logger.print("[ERROR] Invalid binding_pockets.")
        return False

    for item in pockets:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid binding_pockets entry.")
            return False

        if not _is_number(item.get("binding_pocket_volume")):
            logger.print("[ERROR] Invalid binding_pocket_volume.")
            return False

        if not isinstance(item.get("binding_pocket_sphere_count"), int):
            logger.print("[ERROR] Invalid binding_pocket_sphere_count.")
            return False

        if not _is_number_list(item.get("binding_pocket_center_coordinate"), min_len=3):
            logger.print("[ERROR] Invalid binding_pocket_center_coordinate.")
            return False

        if not _is_number_list(item.get("binding_pocket_box_size"), min_len=3):
            logger.print("[ERROR] Invalid binding_pocket_box_size.")
            return False

        residues = item.get("residues")
        if not isinstance(residues, list):
            logger.print("[ERROR] Invalid binding pocket residues.")
            return False

        for residue in residues:
            if not isinstance(residue, dict) or not _validate_residue_index_name(residue, logger):
                return False

    return True


def validate_substrate_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("report_type") != "enzywizard_substrate":
        logger.print("[ERROR] substrate report report_type mismatch.")
        return False

    substrates = data.get("substrates")
    if not isinstance(substrates, list):
        logger.print("[ERROR] Invalid substrates.")
        return False

    for item in substrates:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid substrate entry.")
            return False

        for key in ["substrate_name", "substrate_smiles"]:
            if not isinstance(item.get(key), str):
                logger.print(f"[ERROR] Invalid substrate field: {key}")
                return False

        if not _is_int01_list(item.get("substrate_fingerprint_encoding")):
            logger.print("[ERROR] Invalid substrate_fingerprint_encoding.")
            return False

        if not isinstance(item.get("substrate_atom_count"), int):
            logger.print("[ERROR] Invalid substrate_atom_count.")
            return False

        for key in ["substrate_molecular_weight", "substrate_logp"]:
            if not _is_number(item.get(key)):
                logger.print(f"[ERROR] Invalid substrate field: {key}")
                return False

        structures = item.get("substrate_possible_structures")
        if not isinstance(structures, list):
            logger.print("[ERROR] Invalid substrate_possible_structures.")
            return False

        for structure in structures:
            if not isinstance(structure, dict):
                logger.print("[ERROR] Invalid substrate_possible_structures entry.")
                return False

            if not isinstance(structure.get("substrate_structure_name"), str):
                logger.print("[ERROR] Invalid substrate_structure_name.")
                return False

            if not _is_number(structure.get("substrate_structure_energy")):
                logger.print("[ERROR] Invalid substrate_structure_energy.")
                return False

    return True


def validate_dock_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("report_type") != "enzywizard_dock":
        logger.print("[ERROR] dock report report_type mismatch.")
        return False

    result = data.get("enzyme_substrate_docking_result")
    if not isinstance(result, dict):
        logger.print("[ERROR] Invalid enzyme_substrate_docking_result.")
        return False

    for key in ["enzyme_substrate_complex_name", "docked_substrate_names"]:
        if not isinstance(result.get(key), str):
            logger.print(f"[ERROR] Invalid dock field: {key}")
            return False

    if not _is_number(result.get("enzyme_substrate_binding_affinity")):
        logger.print("[ERROR] Invalid enzyme_substrate_binding_affinity.")
        return False

    if not _is_number_list(result.get("docking_box_center_coordinate"), min_len=3):
        logger.print("[ERROR] Invalid docking_box_center_coordinate.")
        return False

    if not _is_number_list(result.get("docking_box_size"), min_len=3):
        logger.print("[ERROR] Invalid docking_box_size.")
        return False

    docked_substrates = result.get("docked_substrates")
    if not isinstance(docked_substrates, list):
        logger.print("[ERROR] Invalid docked_substrates.")
        return False

    for item in docked_substrates:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid docked_substrates entry.")
            return False

        for key in ["docked_substrate_name", "docked_substrate_structure_name"]:
            if not isinstance(item.get(key), str):
                logger.print(f"[ERROR] Invalid docked substrate field: {key}")
                return False

        if not _is_number_list(item.get("docked_substrate_center_coordinate"), min_len=3):
            logger.print("[ERROR] Invalid docked_substrate_center_coordinate.")
            return False

    return True


def validate_interaction_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("report_type") != "enzywizard_interaction":
        logger.print("[ERROR] interaction report report_type mismatch.")
        return False

    interactions = data.get("molecular_interactions")
    stats = data.get("molecular_interaction_statistics")

    if not isinstance(interactions, list):
        logger.print("[ERROR] Invalid molecular_interactions.")
        return False
    if not isinstance(stats, dict):
        logger.print("[ERROR] Invalid molecular_interaction_statistics.")
        return False

    block_name_list = [
        "overall_molecular_interaction_statistics",
        "intra_enzyme_interaction_statistics",
        "enzyme_substrate_interaction_statistics",
    ]

    count_field_list = [
        "hydrogen_bond_count",
        "ionic_bond_count",
        "van_der_waals_contact_count",
        "pi_pi_stacking_count",
        "pi_cation_interaction_count",
        "disulfide_bond_count",
    ]

    for block_name in block_name_list:
        block = stats.get(block_name)
        if not isinstance(block, dict):
            logger.print(f"[ERROR] Invalid interaction statistics block: {block_name}")
            return False

        for sub_key in ["interaction_count", "unique_pair_interaction_count"]:
            sub_block = block.get(sub_key)
            if not isinstance(sub_block, dict):
                logger.print(f"[ERROR] Invalid interaction statistics sub-block: {block_name}.{sub_key}")
                return False

            for count_field in count_field_list:
                if not isinstance(sub_block.get(count_field), int):
                    logger.print(f"[ERROR] Invalid interaction statistic: {block_name}.{sub_key}.{count_field}")
                    return False

    for item in interactions:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid molecular_interactions entry.")
            return False

        interaction_type = item.get("molecular_interaction_type")
        source_node = item.get("source_node")
        target_node = item.get("target_node")

        if interaction_type not in INTERACTION_ORDER:
            logger.print("[ERROR] Invalid molecular_interaction_type.")
            return False

        if not isinstance(source_node, dict) or not isinstance(target_node, dict):
            logger.print("[ERROR] Invalid molecular interaction node.")
            return False

        if not validate_interaction_node(source_node, logger):
            return False
        if not validate_interaction_node(target_node, logger):
            return False

    return True


def validate_interaction_node(node: Dict[str, Any], logger: Logger) -> bool:
    node_type = node.get("node_type")

    if node_type == "residue":
        if not isinstance(node.get("residue_index"), int):
            logger.print("[ERROR] Invalid residue node residue_index.")
            return False
        if not isinstance(node.get("residue_name"), str):
            logger.print("[ERROR] Invalid residue node residue_name.")
            return False
        return True

    if node_type == "substrate":
        if not isinstance(node.get("substrate_index"), int):
            logger.print("[ERROR] Invalid substrate node substrate_index.")
            return False
        if not isinstance(node.get("substrate_name"), str):
            logger.print("[ERROR] Invalid substrate node substrate_name.")
            return False
        return True

    logger.print("[ERROR] Invalid interaction node_type.")
    return False