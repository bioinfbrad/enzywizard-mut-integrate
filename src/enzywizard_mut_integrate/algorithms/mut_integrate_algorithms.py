from __future__ import annotations

import math
from typing import Any, Dict, List, Tuple

from ..utils.logging_utils import Logger
from ..utils.integrate_utils import build_lookup_by_residue,get_clean_new_residue_list

from ..utils.mut_clean_utils import check_amino_acid_substitution, get_muts_from_aas
from ..utils.mut_integrate_utils import synthesize_clean_report_from_mutclean

from ..algorithms.integrate_algorithms import build_overall_statistics,build_integrated_graph
from ..utils.sequence_utils import normalize_aa_name_to_one_letter


def generate_mut_integrate_report(
    cleaned_amino_acid_substitution: str,
    overall_statistics: Dict[str, Any],
    amino_acid_substitution_properties: Dict[str, Any],
    wt_integrated_graph: List[Dict[str, Any]],
    mut_integrated_graph: List[Dict[str, Any]],
) -> Dict[str, Any]:
    return {
        "report_type": "enzywizard_mut_integrate",
        "cleaned_amino_acid_substitution": cleaned_amino_acid_substitution,
        "overall_statistics": overall_statistics,
        "amino_acid_substitution_properties": amino_acid_substitution_properties,
        "wild_type_integrated_graph": wt_integrated_graph,
        "mutant_integrated_graph": mut_integrated_graph,
    }


def integrate_mut_reports(
    mutclean_report: Dict[str, Any],
    wt_report_dict: Dict[str, Dict[str, Any]],
    mut_report_dict: Dict[str, Dict[str, Any]],
    strict: bool,
    logger: Logger,
) -> Dict[str, Any] | None:
    if not isinstance(mutclean_report, dict):
        logger.print("[ERROR] mut_clean_report must be a dict.")
        return None

    if not isinstance(wt_report_dict, dict):
        logger.print("[ERROR] wt_report_dict must be a dict.")
        return None

    if not isinstance(mut_report_dict, dict):
        logger.print("[ERROR] mut_report_dict must be a dict.")
        return None

    wt_clean_report = synthesize_clean_report_from_mutclean(mutclean_report, "wt", logger)
    if wt_clean_report is None:
        return None

    mut_clean_report = synthesize_clean_report_from_mutclean(mutclean_report, "mut", logger)
    if mut_clean_report is None:
        return None

    wt_report_dict_full = dict(wt_report_dict)
    wt_report_dict_full["enzywizard_clean"] = wt_clean_report

    mut_report_dict_full = dict(mut_report_dict)
    mut_report_dict_full["enzywizard_clean"] = mut_clean_report

    cleaned_amino_acid_substitution = mutclean_report.get("cleaned_amino_acid_substitution")
    if not isinstance(cleaned_amino_acid_substitution, str) or cleaned_amino_acid_substitution.strip() == "":
        logger.print("[ERROR] Invalid cleaned_amino_acid_substitution in mut_clean report.")
        return None
    cleaned_amino_acid_substitution = cleaned_amino_acid_substitution.strip()

    wt_new_residue_list = get_clean_new_residue_list(wt_clean_report, logger)
    if wt_new_residue_list is None:
        return None

    mut_new_residue_list = get_clean_new_residue_list(mut_clean_report, logger)
    if mut_new_residue_list is None:
        return None

    if len(wt_new_residue_list) == 0:
        logger.print("[ERROR] wild_type_residue_mapping_old_to_new cannot be empty.")
        return None

    if len(mut_new_residue_list) == 0:
        logger.print("[ERROR] mutant_residue_mapping_old_to_new cannot be empty.")
        return None

    if not check_amino_acid_substitution(
        mutation=cleaned_amino_acid_substitution,
        wt_length=len(wt_new_residue_list),
        mut_length=len(mut_new_residue_list),
        logger=logger,
    ):
        return None

    overall_statistics = build_mut_overall_statistics(
        wt_report_dict=wt_report_dict_full,
        mut_report_dict=mut_report_dict_full,
        strict=strict,
        logger=logger,
    )
    if overall_statistics is None:
        return None

    amino_acid_substitution_properties = build_mutation_site_features(
        mutclean_report=mutclean_report,
        wt_report_dict=wt_report_dict_full,
        mut_report_dict=mut_report_dict_full,
        logger=logger,
    )
    if amino_acid_substitution_properties is None:
        return None

    wt_integrated_graph, mut_integrated_graph = build_mut_integrated_graphs(
        wt_report_dict=wt_report_dict_full,
        mut_report_dict=mut_report_dict_full,
        strict=strict,
        logger=logger,
    )
    if wt_integrated_graph is None or mut_integrated_graph is None:
        return None

    return generate_mut_integrate_report(
        cleaned_amino_acid_substitution=cleaned_amino_acid_substitution,
        overall_statistics=overall_statistics,
        amino_acid_substitution_properties=amino_acid_substitution_properties,
        wt_integrated_graph=wt_integrated_graph,
        mut_integrated_graph=mut_integrated_graph,
    )


def build_mut_overall_statistics(
    wt_report_dict: Dict[str, Dict[str, Any]],
    mut_report_dict: Dict[str, Dict[str, Any]],
    strict: bool,
    logger: Logger,
) -> Dict[str, Any] | None:
    wt_has_any = any(
        key in wt_report_dict
        for key in [
            "enzywizard_aaprops",
            "enzywizard_hydrocluster",
            "enzywizard_disorder",
            "enzywizard_pocket",
            "enzywizard_energy",
            "enzywizard_dock",
            "enzywizard_interaction",
        ]
    )
    mut_has_any = any(
        key in mut_report_dict
        for key in [
            "enzywizard_aaprops",
            "enzywizard_hydrocluster",
            "enzywizard_disorder",
            "enzywizard_pocket",
            "enzywizard_energy",
            "enzywizard_dock",
            "enzywizard_interaction",
        ]
    )

    if not wt_has_any and not mut_has_any:
        if strict:
            logger.print("[ERROR] overall_statistics cannot be empty in strict mode.")
            return None
        return {}

    wt_stats = build_overall_statistics(wt_report_dict, strict=False, logger=logger)
    if wt_stats is None:
        return None

    mut_stats = build_overall_statistics(mut_report_dict, strict=False, logger=logger)
    if mut_stats is None:
        return None

    overall_statistics = reorder_mut_overall_statistics(wt_stats, mut_stats)
    return overall_statistics


def build_mut_integrated_graphs(
    wt_report_dict: Dict[str, Dict[str, Any]],
    mut_report_dict: Dict[str, Dict[str, Any]],
    strict: bool,
    logger: Logger,
) -> Tuple[List[Dict[str, Any]] | None, List[Dict[str, Any]] | None]:
    wt_integrated_graph = build_integrated_graph(wt_report_dict, strict=strict, logger=logger)
    if wt_integrated_graph is None:
        return None, None

    mut_integrated_graph = build_integrated_graph(mut_report_dict, strict=strict, logger=logger)
    if mut_integrated_graph is None:
        return None, None

    return wt_integrated_graph, mut_integrated_graph

def reorder_mutation_site_features(data: Dict[str, Any]) -> Dict[str, Any]:
    ordered: Dict[str, Any] = {}

    field_order = [
        "residue_name",
        "residue_name_one_hot_encoding",
        "residue_chemical_classification",
        "residue_chemical_classification_one_hot_encoding",
        "residue_secondary_structure",
        "residue_secondary_structure_one_hot_encoding",
        "residue_relative_solvent_accessibility",
        "residue_backbone_phi_angle",
        "residue_backbone_psi_angle",
        "residue_net_charge",
        "residue_pka",
        "residue_volume",
        "residue_hydrophobicity",
        "residue_molecular_weight",
        "residue_isoelectric_point",
        "residue_root_mean_square_fluctuation",
        "residue_sequence_conservation_score",
    ]

    for field_name in field_order:
        wt_key = f"wild_type_{field_name}"
        mut_key = f"mutant_{field_name}"
        diff_key = f"difference_{field_name}"

        if wt_key in data:
            ordered[wt_key] = data[wt_key]
        if mut_key in data:
            ordered[mut_key] = data[mut_key]
        if diff_key in data:
            ordered[diff_key] = data[diff_key]

    for key, value in data.items():
        if key not in ordered:
            ordered[key] = value

    return ordered

def build_mutation_site_features(
    mutclean_report: Dict[str, Any],
    wt_report_dict: Dict[str, Dict[str, Any]],
    mut_report_dict: Dict[str, Dict[str, Any]],
    logger: Logger,
) -> Dict[str, Any] | None:
    result: Dict[str, Any] = {}

    wt_clean_report = wt_report_dict.get("enzywizard_clean")
    mut_clean_report = mut_report_dict.get("enzywizard_clean")
    if not isinstance(wt_clean_report, dict) or not isinstance(mut_clean_report, dict):
        logger.print("[ERROR] Missing synthesized clean report while building amino_acid_substitution_properties.")
        return None

    cleaned_amino_acid_substitution = mutclean_report.get("cleaned_amino_acid_substitution")
    if not isinstance(cleaned_amino_acid_substitution, str) or cleaned_amino_acid_substitution.strip() == "":
        logger.print("[ERROR] Invalid cleaned_amino_acid_substitution.")
        return None

    wt_new_residue_list = get_clean_new_residue_list(wt_clean_report, logger)
    if wt_new_residue_list is None:
        return None

    mut_new_residue_list = get_clean_new_residue_list(mut_clean_report, logger)
    if mut_new_residue_list is None:
        return None

    wt_residue_lookup = {
        int(item["residue_index"]): normalize_aa_name_to_one_letter(item["residue_name"])
        for item in wt_new_residue_list
    }
    mut_residue_lookup = {
        int(item["residue_index"]): normalize_aa_name_to_one_letter(item["residue_name"])
        for item in mut_new_residue_list
    }

    mut_list = get_muts_from_aas(cleaned_amino_acid_substitution)

    wt_aaprops_lookup = None
    mut_aaprops_lookup = None
    wt_flexibility_lookup = None
    mut_flexibility_lookup = None
    wt_conservation_lookup = None
    mut_conservation_lookup = None

    if "enzywizard_aaprops" in wt_report_dict:
        wt_aaprops_lookup = build_lookup_by_residue(
            wt_report_dict["enzywizard_aaprops"]["amino_acid_residue_properties"],
            "residue_index",
            "residue_name",
            logger,
        )
        if wt_aaprops_lookup is None:
            return None

    if "enzywizard_aaprops" in mut_report_dict:
        mut_aaprops_lookup = build_lookup_by_residue(
            mut_report_dict["enzywizard_aaprops"]["amino_acid_residue_properties"],
            "residue_index",
            "residue_name",
            logger,
        )
        if mut_aaprops_lookup is None:
            return None

    if "enzywizard_flexibility" in wt_report_dict:
        wt_flexibility_lookup = build_lookup_by_residue(
            wt_report_dict["enzywizard_flexibility"]["protein_flexibility"],
            "residue_index",
            "residue_name",
            logger,
        )
        if wt_flexibility_lookup is None:
            return None

    if "enzywizard_flexibility" in mut_report_dict:
        mut_flexibility_lookup = build_lookup_by_residue(
            mut_report_dict["enzywizard_flexibility"]["protein_flexibility"],
            "residue_index",
            "residue_name",
            logger,
        )
        if mut_flexibility_lookup is None:
            return None

    if "enzywizard_conservation" in wt_report_dict:
        wt_conservation_lookup = build_lookup_by_residue(
            wt_report_dict["enzywizard_conservation"]["sequence_conservation_scores"],
            "residue_index",
            "residue_name",
            logger,
        )
        if wt_conservation_lookup is None:
            return None

    if "enzywizard_conservation" in mut_report_dict:
        mut_conservation_lookup = build_lookup_by_residue(
            mut_report_dict["enzywizard_conservation"]["sequence_conservation_scores"],
            "residue_index",
            "residue_name",
            logger,
        )
        if mut_conservation_lookup is None:
            return None

    wt_name_list: List[str] = []
    mut_name_list: List[str] = []

    wt_class_list: List[str] = []
    mut_class_list: List[str] = []

    wt_ss_list: List[str] = []
    mut_ss_list: List[str] = []

    wt_rsa_list: List[float] = []
    mut_rsa_list: List[float] = []

    wt_phi_list: List[float] = []
    mut_phi_list: List[float] = []

    wt_psi_list: List[float] = []
    mut_psi_list: List[float] = []

    wt_net_charge_list: List[float] = []
    mut_net_charge_list: List[float] = []

    wt_pka_list: List[float] = []
    mut_pka_list: List[float] = []

    wt_volume_list: List[float] = []
    mut_volume_list: List[float] = []

    wt_hydrophobicity_list: List[float] = []
    mut_hydrophobicity_list: List[float] = []

    wt_molecular_weight_list: List[float] = []
    mut_molecular_weight_list: List[float] = []

    wt_pi_list: List[float] = []
    mut_pi_list: List[float] = []

    wt_rmsf_list: List[float] = []
    mut_rmsf_list: List[float] = []

    wt_conservation_score_list: List[float] = []
    mut_conservation_score_list: List[float] = []

    wt_residue_name_one_hot_list: List[List[float]] = []
    mut_residue_name_one_hot_list: List[List[float]] = []

    wt_residue_class_one_hot_list: List[List[float]] = []
    mut_residue_class_one_hot_list: List[List[float]] = []

    wt_residue_ss_one_hot_list: List[List[float]] = []
    mut_residue_ss_one_hot_list: List[List[float]] = []

    for wt_aa_expected, pos, mut_aa_expected in mut_list:
        wt_aa_actual = wt_residue_lookup.get(pos)
        mut_aa_actual = mut_residue_lookup.get(pos)

        if wt_aa_actual is None:
            logger.print(f"[ERROR] Mutation position missing in WT cleaned residue list: {pos}")
            return None

        if mut_aa_actual is None:
            logger.print(f"[ERROR] Mutation position missing in MUT cleaned residue list: {pos}")
            return None

        if wt_aa_actual != wt_aa_expected:
            logger.print(
                f"[ERROR] WT amino acid mismatch at mutation position {pos}: "
                f"expected {wt_aa_expected}, got {wt_aa_actual}"
            )
            return None

        if mut_aa_actual != mut_aa_expected:
            logger.print(
                f"[ERROR] MUT amino acid mismatch at mutation position {pos}: "
                f"expected {mut_aa_expected}, got {mut_aa_actual}"
            )
            return None

        wt_name_list.append(wt_aa_actual)
        mut_name_list.append(mut_aa_actual)

        wt_key = (pos, wt_aa_actual)
        mut_key = (pos, mut_aa_actual)

        if wt_aaprops_lookup is not None:
            wt_item = wt_aaprops_lookup.get(wt_key)
            if wt_item is None:
                logger.print(f"[ERROR] Missing WT amino_acid_residue_properties entry for mutation site: {pos} {wt_aa_actual}")
                return None
            wt_class_list.append(wt_item["residue_chemical_classification"])
            wt_ss_list.append(wt_item["residue_secondary_structure"])
            wt_rsa_list.append(float(wt_item["residue_relative_solvent_accessibility"]))
            wt_phi_list.append(float(wt_item["residue_backbone_phi_angle"]))
            wt_psi_list.append(float(wt_item["residue_backbone_psi_angle"]))
            wt_net_charge_list.append(float(wt_item["residue_net_charge"]))
            wt_pka_list.append(float(wt_item["residue_pka"]))
            wt_volume_list.append(float(wt_item["residue_volume"]))
            wt_hydrophobicity_list.append(float(wt_item["residue_hydrophobicity"]))
            wt_molecular_weight_list.append(float(wt_item["residue_molecular_weight"]))
            wt_pi_list.append(float(wt_item["residue_isoelectric_point"]))
            wt_residue_name_one_hot_list.append([float(x) for x in wt_item["residue_name_one_hot_encoding"]])
            wt_residue_class_one_hot_list.append([float(x) for x in wt_item["residue_chemical_classification_one_hot_encoding"]])
            wt_residue_ss_one_hot_list.append([float(x) for x in wt_item["residue_secondary_structure_one_hot_encoding"]])

        if mut_aaprops_lookup is not None:
            mut_item = mut_aaprops_lookup.get(mut_key)
            if mut_item is None:
                logger.print(f"[ERROR] Missing MUT amino_acid_residue_properties entry for mutation site: {pos} {mut_aa_actual}")
                return None
            mut_class_list.append(mut_item["residue_chemical_classification"])
            mut_ss_list.append(mut_item["residue_secondary_structure"])
            mut_rsa_list.append(float(mut_item["residue_relative_solvent_accessibility"]))
            mut_phi_list.append(float(mut_item["residue_backbone_phi_angle"]))
            mut_psi_list.append(float(mut_item["residue_backbone_psi_angle"]))
            mut_net_charge_list.append(float(mut_item["residue_net_charge"]))
            mut_pka_list.append(float(mut_item["residue_pka"]))
            mut_volume_list.append(float(mut_item["residue_volume"]))
            mut_hydrophobicity_list.append(float(mut_item["residue_hydrophobicity"]))
            mut_molecular_weight_list.append(float(mut_item["residue_molecular_weight"]))
            mut_pi_list.append(float(mut_item["residue_isoelectric_point"]))
            mut_residue_name_one_hot_list.append([float(x) for x in mut_item["residue_name_one_hot_encoding"]])
            mut_residue_class_one_hot_list.append([float(x) for x in mut_item["residue_chemical_classification_one_hot_encoding"]])
            mut_residue_ss_one_hot_list.append([float(x) for x in mut_item["residue_secondary_structure_one_hot_encoding"]])

        if wt_flexibility_lookup is not None:
            wt_rmsf_item = wt_flexibility_lookup.get(wt_key)
            if wt_rmsf_item is None:
                logger.print(f"[ERROR] Missing WT RMSF entry for mutation site: {pos} {wt_aa_actual}")
                return None
            wt_rmsf_list.append(float(wt_rmsf_item["residue_root_mean_square_fluctuation"]))

        if mut_flexibility_lookup is not None:
            mut_rmsf_item = mut_flexibility_lookup.get(mut_key)
            if mut_rmsf_item is None:
                logger.print(f"[ERROR] Missing MUT RMSF entry for mutation site: {pos} {mut_aa_actual}")
                return None
            mut_rmsf_list.append(float(mut_rmsf_item["residue_root_mean_square_fluctuation"]))

        if wt_conservation_lookup is not None:
            wt_cons_item = wt_conservation_lookup.get(wt_key)
            if wt_cons_item is None:
                logger.print(f"[ERROR] Missing WT conservation entry for mutation site: {pos} {wt_aa_actual}")
                return None
            wt_conservation_score_list.append(float(wt_cons_item["normalized_shannon_information_content"]))

        if mut_conservation_lookup is not None:
            mut_cons_item = mut_conservation_lookup.get(mut_key)
            if mut_cons_item is None:
                logger.print(f"[ERROR] Missing MUT conservation entry for mutation site: {pos} {mut_aa_actual}")
                return None
            mut_conservation_score_list.append(float(mut_cons_item["normalized_shannon_information_content"]))

    if len(wt_name_list) > 0:
        result["wild_type_residue_name"] = ",".join(wt_name_list)
    if len(mut_name_list) > 0:
        result["mutant_residue_name"] = ",".join(mut_name_list)

    if len(wt_class_list) > 0:
        result["wild_type_residue_chemical_classification"] = ",".join(wt_class_list)
    if len(mut_class_list) > 0:
        result["mutant_residue_chemical_classification"] = ",".join(mut_class_list)

    if len(wt_ss_list) > 0:
        result["wild_type_residue_secondary_structure"] = ",".join(wt_ss_list)
    if len(mut_ss_list) > 0:
        result["mutant_residue_secondary_structure"] = ",".join(mut_ss_list)

    _write_mutation_numeric_triplet(result, "residue_relative_solvent_accessibility", wt_rsa_list, mut_rsa_list)
    _write_mutation_angle_triplet(result, "residue_backbone_phi_angle", wt_phi_list, mut_phi_list)
    _write_mutation_angle_triplet(result, "residue_backbone_psi_angle", wt_psi_list, mut_psi_list)
    _write_mutation_numeric_triplet(result, "residue_net_charge", wt_net_charge_list, mut_net_charge_list)
    _write_mutation_numeric_triplet(result, "residue_pka", wt_pka_list, mut_pka_list)
    _write_mutation_numeric_triplet(result, "residue_volume", wt_volume_list, mut_volume_list)
    _write_mutation_numeric_triplet(result, "residue_hydrophobicity", wt_hydrophobicity_list, mut_hydrophobicity_list)
    _write_mutation_numeric_triplet(result, "residue_molecular_weight", wt_molecular_weight_list, mut_molecular_weight_list)
    _write_mutation_numeric_triplet(result, "residue_isoelectric_point", wt_pi_list, mut_pi_list)
    _write_mutation_numeric_triplet(result, "residue_root_mean_square_fluctuation", wt_rmsf_list, mut_rmsf_list)
    _write_mutation_numeric_triplet(result, "residue_sequence_conservation_score", wt_conservation_score_list, mut_conservation_score_list)

    _write_mutation_vector_triplet(result, "residue_name_one_hot_encoding", wt_residue_name_one_hot_list, mut_residue_name_one_hot_list)
    _write_mutation_vector_triplet(result, "residue_chemical_classification_one_hot_encoding", wt_residue_class_one_hot_list, mut_residue_class_one_hot_list)
    _write_mutation_vector_triplet(result, "residue_secondary_structure_one_hot_encoding", wt_residue_ss_one_hot_list, mut_residue_ss_one_hot_list)

    result = reorder_mutation_site_features(result)

    return result



def reorder_mut_overall_statistics(
    wt_stats: Dict[str, Any],
    mut_stats: Dict[str, Any],
) -> Dict[str, Any]:
    ordered: Dict[str, Any] = {}

    field_order = [
        "residue_name_count",
        "residue_chemical_classification_count",
        "residue_secondary_structure_count",
        "hydrophobic_cluster_count",
        "max_hydrophobic_cluster_area",
        "total_hydrophobic_cluster_area",
        "disordered_region_count",
        "max_disordered_region_length",
        "total_disordered_region_length",
        "binding_pocket_count",
        "max_binding_pocket_volume",
        "total_binding_pocket_volume",
        "total_potential_energy",
        "harmonic_bond_potential_energy",
        "harmonic_angle_potential_energy",
        "custom_bond_potential_energy",
        "custom_torsion_potential_energy",
        "custom_nonbonded_potential_energy",
        "nonbonded_potential_energy",
        "periodic_torsion_potential_energy",
        "cmap_torsion_potential_energy",
        "enzyme_substrate_binding_affinity",
        "hydrogen_bond_count",
        "ionic_bond_count",
        "van_der_waals_contact_count",
        "pi_pi_stacking_count",
        "pi_cation_interaction_count",
        "disulfide_bond_count",
    ]

    for field_name in field_order:
        wt_has = field_name in wt_stats
        mut_has = field_name in mut_stats

        if wt_has:
            ordered[f"wild_type_{field_name}"] = wt_stats[field_name]
        if mut_has:
            ordered[f"mutant_{field_name}"] = mut_stats[field_name]

        if wt_has and mut_has:
            diff_value = _diff_scalar_or_list(wt_stats[field_name], mut_stats[field_name])
            if diff_value is not None:
                ordered[f"difference_{field_name}"] = diff_value

    return ordered




def _diff_scalar_or_list(wt_value: Any, mut_value: Any) -> Any | None:
    if isinstance(wt_value, (int, float)) and not isinstance(wt_value, bool):
        if isinstance(mut_value, (int, float)) and not isinstance(mut_value, bool):
            return mut_value - wt_value

    if isinstance(wt_value, list) and isinstance(mut_value, list):
        if len(wt_value) != len(mut_value):
            return None
        out: List[float] = []
        for x, y in zip(wt_value, mut_value):
            if not isinstance(x, (int, float)) or isinstance(x, bool):
                return None
            if not isinstance(y, (int, float)) or isinstance(y, bool):
                return None
            out.append(float(y) - float(x))
        return out

    return None


def _mean_number_list(value_list: List[float]) -> float | None:
    if len(value_list) == 0:
        return None
    return sum(value_list) / float(len(value_list))


def _mean_vector_list(vector_list: List[List[float]]) -> List[float] | None:
    if len(vector_list) == 0:
        return None
    dim = len(vector_list[0])
    if dim == 0:
        return None
    for vec in vector_list:
        if len(vec) != dim:
            return None
    out: List[float] = []
    for i in range(dim):
        out.append(sum(vec[i] for vec in vector_list) / float(len(vector_list)))
    return out


def _circular_mean_deg(value_list: List[float]) -> float | None:
    if len(value_list) == 0:
        return None
    sin_sum = 0.0
    cos_sum = 0.0
    for value in value_list:
        rad = math.radians(value)
        sin_sum += math.sin(rad)
        cos_sum += math.cos(rad)
    if abs(sin_sum) < 1e-12 and abs(cos_sum) < 1e-12:
        return 0.0
    mean_rad = math.atan2(sin_sum, cos_sum)
    mean_deg = math.degrees(mean_rad)
    return _normalize_angle_deg(mean_deg)


def _normalize_angle_deg(value: float) -> float:
    out = ((value + 180.0) % 360.0) - 180.0
    if out == -180.0:
        return 180.0
    return out


def _circular_diff_deg(wt_value: float, mut_value: float) -> float:
    return _normalize_angle_deg(mut_value - wt_value)


def _write_mutation_numeric_triplet(
    out: Dict[str, Any],
    field_name: str,
    wt_list: List[float],
    mut_list: List[float],
) -> None:
    wt_mean = _mean_number_list(wt_list)
    mut_mean = _mean_number_list(mut_list)

    if wt_mean is not None:
        out[f"wild_type_{field_name}"] = wt_mean
    if mut_mean is not None:
        out[f"mutant_{field_name}"] = mut_mean
    if wt_mean is not None and mut_mean is not None:
        out[f"difference_{field_name}"] = mut_mean - wt_mean


def _write_mutation_angle_triplet(
    out: Dict[str, Any],
    field_name: str,
    wt_list: List[float],
    mut_list: List[float],
) -> None:
    wt_mean = _circular_mean_deg(wt_list)
    mut_mean = _circular_mean_deg(mut_list)

    if wt_mean is not None:
        out[f"wild_type_{field_name}"] = wt_mean
    if mut_mean is not None:
        out[f"mutant_{field_name}"] = mut_mean
    if wt_mean is not None and mut_mean is not None:
        out[f"difference_{field_name}"] = _circular_diff_deg(wt_mean, mut_mean)


def _write_mutation_vector_triplet(
    out: Dict[str, Any],
    field_name: str,
    wt_vector_list: List[List[float]],
    mut_vector_list: List[List[float]],
) -> None:
    wt_mean = _mean_vector_list(wt_vector_list)
    mut_mean = _mean_vector_list(mut_vector_list)

    if wt_mean is not None:
        out[f"wild_type_{field_name}"] = wt_mean
    if mut_mean is not None:
        out[f"mutant_{field_name}"] = mut_mean
    if wt_mean is not None and mut_mean is not None:
        out[f"difference_{field_name}"] = [y - x for x, y in zip(wt_mean, mut_mean)]