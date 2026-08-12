from __future__ import annotations

from pathlib import Path
import math
from typing import Any, Dict, List, Tuple
import json
import re

from ..utils.logging_utils import Logger
from ..utils.IO_utils import write_json_from_dict_inline_leaf_lists

from ..utils.integrate_utils import (
    SUPPORTED_OUTPUT_TYPES,
    get_disorder_membership_set,
    get_hydrophobic_cluster_membership_set,
    get_pocket_membership_set,
    load_json_file,
    list_json_files,
    residue_key,
    split_integrated_graph_entries,
    validate_clean_report,
    validate_aaprops_report,
    validate_hydrocluster_report,
    validate_energy_report,
    validate_flexibility_report,
    validate_disorder_report,
    validate_conservation_report,
    validate_embedding_report,
    validate_pocket_report,
    validate_substrate_report,
    validate_dock_report,
    validate_interaction_report,
)
from ..utils.sequence_utils import normalize_aa_name_to_one_letter


MUT_INTEGRATE_SUPPORTED_OUTPUT_TYPES = set(SUPPORTED_OUTPUT_TYPES) | {"enzywizard_mut_clean"}

MUT_INTEGRATE_SIDE_OUTPUT_TYPES = {
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


def save_mut_integrate_json(report: Dict[str, Any] | List[Any], output_path: str | Path, logger: Logger) -> bool:
    try:
        write_json_from_dict_inline_leaf_lists(report, output_path)
        return True
    except Exception as e:
        logger.print(f"[ERROR] Failed to save mut-integrate JSON: {e}")
        return False


def extract_wt_mut_protein_names_from_mutclean_report_path(
    mutclean_report_path: str | Path,
    logger: Logger,
) -> Tuple[str, str] | None:

    try:
        name = Path(mutclean_report_path).name
        m = re.fullmatch(r"mut_clean_report_(.+)\.json", name)
        if m is None:
            logger.print(
                "[ERROR] mut_clean_report file name must match "
                "mut_clean_report_{wt_protein_name}_to_{mut_protein_name}.json"
            )
            return None

        body = m.group(1).strip()
        if body == "":
            logger.print(f"[ERROR] Invalid mut_clean_report file name: {name}")
            return None

        split_token = "_to_"
        if split_token not in body:
            logger.print(
                f"[ERROR] Invalid mut_clean_report file name. Cannot find '_to_' separator: {name}"
            )
            return None

        wt_protein_name, mut_protein_name = body.split(split_token, 1)
        wt_protein_name = wt_protein_name.strip()
        mut_protein_name = mut_protein_name.strip()

        if wt_protein_name == "" or mut_protein_name == "":
            logger.print(
                f"[ERROR] Invalid WT or MUT protein name in mut_clean_report file name: {name}"
            )
            return None

        return wt_protein_name, mut_protein_name
    except Exception as e:
        logger.print(f"[ERROR] Failed to parse mut_clean_report file name: {e}")
        return None


def get_mut_integrate_supported_output_type(data: Dict[str, Any], logger: Logger) -> str | None:
    report_type = data.get("report_type")
    if not isinstance(report_type, str):
        logger.print("[ERROR] Missing or invalid report_type.")
        return None
    if report_type not in MUT_INTEGRATE_SUPPORTED_OUTPUT_TYPES:
        logger.print(f"[ERROR] Unsupported report_type for mut_integrate: {report_type}")
        return None
    return report_type


def validate_mut_integrate_report_by_type(data: Dict[str, Any], logger: Logger) -> bool:
    report_type = get_mut_integrate_supported_output_type(data, logger)
    if report_type is None:
        return False

    if report_type == "enzywizard_mut_clean":
        return validate_mutclean_report(data, logger)
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


def validate_mutclean_report(data: Dict[str, Any], logger: Logger) -> bool:
    if data.get("report_type") != "enzywizard_mut_clean":
        logger.print("[ERROR] mut_clean report report_type mismatch.")
        return False

    amino_acid_substitution = data.get("amino_acid_substitution")
    cleaned_amino_acid_substitution = data.get("cleaned_amino_acid_substitution")
    wt_mapping = data.get("wild_type_residue_mapping_old_to_new")
    wt_stats = data.get("wild_type_clean_statistics")
    mut_mapping = data.get("mutant_residue_mapping_old_to_new")
    mut_stats = data.get("mutant_clean_statistics")

    if not isinstance(amino_acid_substitution, str) or amino_acid_substitution.strip() == "":
        logger.print("[ERROR] Invalid amino_acid_substitution in mut_clean report.")
        return False

    if not isinstance(cleaned_amino_acid_substitution, str) or cleaned_amino_acid_substitution.strip() == "":
        logger.print("[ERROR] Invalid cleaned_amino_acid_substitution in mut_clean report.")
        return False

    if not isinstance(wt_mapping, list):
        logger.print("[ERROR] Invalid wild_type_residue_mapping_old_to_new.")
        return False

    if not isinstance(mut_mapping, list):
        logger.print("[ERROR] Invalid mutant_residue_mapping_old_to_new.")
        return False

    if not isinstance(wt_stats, dict):
        logger.print("[ERROR] Invalid wild_type_clean_statistics.")
        return False

    if not isinstance(mut_stats, dict):
        logger.print("[ERROR] Invalid mutant_clean_statistics.")
        return False

    required_stat_keys = [
        "removed_heterogen_count",
        "standardized_residue_name_count",
        "repaired_residue_count",
        "added_heavy_atom_count",
        "added_hydrogen_atom_count",
        "retained_residue_count",
    ]

    for stats, stats_name in [(wt_stats, "wild_type_clean_statistics"), (mut_stats, "mutant_clean_statistics")]:
        for k in required_stat_keys:
            if not isinstance(stats.get(k), int):
                logger.print(f"[ERROR] Invalid {stats_name} field: {k}")
                return False

    for mapping, mapping_name in [
        (wt_mapping, "wild_type_residue_mapping_old_to_new"),
        (mut_mapping, "mutant_residue_mapping_old_to_new"),
    ]:
        for item in mapping:
            if not isinstance(item, dict):
                logger.print(f"[ERROR] Invalid residue mapping entry in {mapping_name}.")
                return False

            old_residue = item.get("old_residue")
            new_residue = item.get("new_residue")

            if not isinstance(old_residue, dict) or not isinstance(new_residue, dict):
                logger.print(f"[ERROR] Invalid residue mapping structure in {mapping_name}.")
                return False

            for residue_item in [old_residue, new_residue]:
                if not isinstance(residue_item.get("residue_index"), int):
                    logger.print(f"[ERROR] Invalid residue_index in {mapping_name}.")
                    return False
                if not isinstance(residue_item.get("residue_name"), str) or residue_item["residue_name"].strip() == "":
                    logger.print(f"[ERROR] Invalid residue_name in {mapping_name}.")
                    return False
                if not isinstance(residue_item.get("hydrogen_atom_count"), int):
                    logger.print(f"[ERROR] Invalid hydrogen_atom_count in {mapping_name}.")
                    return False

    return True


def synthesize_clean_report_from_mutclean(
    mutclean_report: Dict[str, Any],
    side: str,
    logger: Logger,
) -> Dict[str, Any] | None:
    if not isinstance(mutclean_report, dict):
        logger.print("[ERROR] mut_clean_report must be a dict.")
        return None

    if side == "wt":
        mapping = mutclean_report.get("wild_type_residue_mapping_old_to_new")
        stats = mutclean_report.get("wild_type_clean_statistics")
    elif side == "mut":
        mapping = mutclean_report.get("mutant_residue_mapping_old_to_new")
        stats = mutclean_report.get("mutant_clean_statistics")
    else:
        logger.print("[ERROR] side must be 'wt' or 'mut'.")
        return None

    if not isinstance(mapping, list):
        logger.print(f"[ERROR] Missing mapping for side={side} in mut_clean_report.")
        return None

    if not isinstance(stats, dict):
        logger.print(f"[ERROR] Missing clean statistics for side={side} in mut_clean_report.")
        return None

    clean_report = {
        "report_type": "enzywizard_clean",
        "residue_mapping_old_to_new": mapping,
        "clean_statistics": stats,
    }

    if not validate_clean_report(clean_report, logger):
        logger.print(f"[ERROR] Failed to synthesize clean report from mut_clean for side={side}.")
        return None

    return clean_report


def build_mutation_site_distance_features(
    mut_list: List[Tuple[str, int, str]],
    wt_residue_lookup: Dict[int, str],
    mut_residue_lookup: Dict[int, str],
    wt_report_dict: Dict[str, Dict[str, Any]],
    mut_report_dict: Dict[str, Dict[str, Any]],
    logger: Logger,
) -> Dict[str, Any] | None:
    wt_features = build_single_side_mutation_site_distance_features(
        mut_list=mut_list,
        residue_lookup=wt_residue_lookup,
        report_dict=wt_report_dict,
        logger=logger,
    )
    if wt_features is None:
        return None

    mut_features = build_single_side_mutation_site_distance_features(
        mut_list=mut_list,
        residue_lookup=mut_residue_lookup,
        report_dict=mut_report_dict,
        logger=logger,
    )
    if mut_features is None:
        return None

    result: Dict[str, Any] = {}
    field_order = [
        "mutation_site_distance_to_centroid",
        "mutation_site_distance_to_nearest_binding_pocket",
        "mutation_site_distance_to_nearest_hydrophobic_cluster",
        "mutation_site_distance_to_nearest_disordered_region",
        "mutation_site_distance_to_nearest_substrate",
    ]

    for field_name in field_order:
        wt_has = field_name in wt_features
        mut_has = field_name in mut_features

        if wt_has:
            result[f"wild_type_{field_name}"] = wt_features[field_name]
        if mut_has:
            result[f"mutant_{field_name}"] = mut_features[field_name]
        if wt_has and mut_has:
            result[f"difference_{field_name}"] = mut_features[field_name] - wt_features[field_name]

    return result


def build_single_side_mutation_site_distance_features(
    mut_list: List[Tuple[str, int, str]],
    residue_lookup: Dict[int, str],
    report_dict: Dict[str, Dict[str, Any]],
    logger: Logger,
) -> Dict[str, float] | None:
    aaprops_report = report_dict.get("enzywizard_aaprops")
    if aaprops_report is None:
        return {}

    residue_properties = aaprops_report.get("amino_acid_residue_properties")
    if not isinstance(residue_properties, list):
        logger.print("[ERROR] amino_acid_residue_properties must be a list.")
        return None

    coordinate_lookup: Dict[Tuple[int, str], List[float]] = {}
    all_coordinates: List[List[float]] = []

    for item in residue_properties:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid amino_acid_residue_properties entry.")
            return None

        residue_index = item.get("residue_index")
        residue_name = item.get("residue_name")
        coordinate = item.get("residue_alpha_carbon_coordinate")
        if not isinstance(residue_index, int) or not isinstance(residue_name, str):
            continue
        clean_coordinate = _normalize_coordinate(coordinate)
        if clean_coordinate is None:
            continue

        normalized_name = normalize_aa_name_to_one_letter(residue_name)
        coordinate_lookup[residue_key(residue_index, normalized_name)] = clean_coordinate
        all_coordinates.append(clean_coordinate)

    if len(all_coordinates) == 0:
        return {}

    mutation_coordinates: List[List[float]] = []
    for _, pos, _ in mut_list:
        residue_name = residue_lookup.get(pos)
        if residue_name is None:
            logger.print(f"[ERROR] Mutation position missing in cleaned residue list: {pos}")
            return None
        coordinate = coordinate_lookup.get(residue_key(pos, residue_name))
        if coordinate is not None:
            mutation_coordinates.append(coordinate)

    if len(mutation_coordinates) == 0:
        return {}

    result: Dict[str, float] = {}
    centroid = _coordinate_centroid(all_coordinates)
    result["mutation_site_distance_to_centroid"] = _mean_distance_to_targets(mutation_coordinates, [centroid])

    pocket_report = report_dict.get("enzywizard_pocket")
    if pocket_report is not None:
        pocket_membership = get_pocket_membership_set(pocket_report, logger)
        if pocket_membership is None:
            return None
        pocket_coordinates = _coordinates_for_membership(pocket_membership, coordinate_lookup)
        nearest_distance = _mean_nearest_distance(mutation_coordinates, pocket_coordinates)
        if nearest_distance is not None:
            result["mutation_site_distance_to_nearest_binding_pocket"] = nearest_distance

    hydro_report = report_dict.get("enzywizard_hydrocluster")
    if hydro_report is not None:
        hydro_membership = get_hydrophobic_cluster_membership_set(hydro_report, logger)
        if hydro_membership is None:
            return None
        hydro_coordinates = _coordinates_for_membership(hydro_membership, coordinate_lookup)
        nearest_distance = _mean_nearest_distance(mutation_coordinates, hydro_coordinates)
        if nearest_distance is not None:
            result["mutation_site_distance_to_nearest_hydrophobic_cluster"] = nearest_distance

    disorder_report = report_dict.get("enzywizard_disorder")
    if disorder_report is not None:
        disorder_membership = get_disorder_membership_set(disorder_report, logger)
        if disorder_membership is None:
            return None
        disorder_coordinates = _coordinates_for_membership(disorder_membership, coordinate_lookup)
        nearest_distance = _mean_nearest_distance(mutation_coordinates, disorder_coordinates)
        if nearest_distance is not None:
            result["mutation_site_distance_to_nearest_disordered_region"] = nearest_distance

    substrate_coordinates = _get_docked_substrate_center_coordinates(report_dict.get("enzywizard_dock"))
    nearest_distance = _mean_nearest_distance(mutation_coordinates, substrate_coordinates)
    if nearest_distance is not None:
        result["mutation_site_distance_to_nearest_substrate"] = nearest_distance

    return result


def _coordinates_for_membership(
    membership_set: set[Tuple[int, str]],
    coordinate_lookup: Dict[Tuple[int, str], List[float]],
) -> List[List[float]]:
    out: List[List[float]] = []
    for key in membership_set:
        coordinate = coordinate_lookup.get(key)
        if coordinate is not None:
            out.append(coordinate)
    return out


def _get_docked_substrate_center_coordinates(dock_report: Dict[str, Any] | None) -> List[List[float]]:
    if not isinstance(dock_report, dict):
        return []
    docking_result = dock_report.get("enzyme_substrate_docking_result")
    if not isinstance(docking_result, dict):
        return []
    docked_substrates = docking_result.get("docked_substrates")
    if not isinstance(docked_substrates, list):
        return []

    out: List[List[float]] = []
    for substrate in docked_substrates:
        if not isinstance(substrate, dict):
            continue
        coordinate = substrate.get("docked_substrate_center_coordinate")
        clean_coordinate = _normalize_coordinate(coordinate)
        if clean_coordinate is not None:
            out.append(clean_coordinate)
    return out


def _mean_distance_to_targets(
    source_coordinates: List[List[float]],
    target_coordinates: List[List[float]],
) -> float:
    distance_list: List[float] = []
    for source_coordinate in source_coordinates:
        for target_coordinate in target_coordinates:
            distance_list.append(_distance(source_coordinate, target_coordinate))
    return sum(distance_list) / float(len(distance_list))


def _mean_nearest_distance(
    source_coordinates: List[List[float]],
    target_coordinates: List[List[float]],
) -> float | None:
    if len(source_coordinates) == 0 or len(target_coordinates) == 0:
        return None

    nearest_distance_list: List[float] = []
    for source_coordinate in source_coordinates:
        nearest_distance_list.append(min(_distance(source_coordinate, target_coordinate) for target_coordinate in target_coordinates))

    return sum(nearest_distance_list) / float(len(nearest_distance_list))


def _normalize_coordinate(value: Any) -> List[float] | None:
    if not isinstance(value, (list, tuple)) or len(value) < 3:
        return None
    out: List[float] = []
    for item in value[:3]:
        if isinstance(item, bool):
            return None
        try:
            out.append(float(item))
        except (TypeError, ValueError):
            return None
    return out


def _coordinate_centroid(coordinate_list: List[List[float]]) -> List[float]:
    return [
        sum(coord[i] for coord in coordinate_list) / float(len(coordinate_list))
        for i in range(3)
    ]


def _distance(coord_1: List[float], coord_2: List[float]) -> float:
    return math.sqrt(sum((coord_1[i] - coord_2[i]) ** 2 for i in range(3)))