from __future__ import annotations
from argparse import Namespace, ArgumentParser
from ..services.mut_integrate_service import run_mut_integrate_service


def add_mut_integrate_parser(parser: ArgumentParser) -> None:
    parser.add_argument("-i","--mut_clean_report_path",required=True,help="Path to the required mut_clean report JSON file.",)
    parser.add_argument("-w","--wt_input_dir",required=True,help="Path to a directory containing wild-type JSON reports to integrate.",)
    parser.add_argument("-m","--mut_input_dir",required=True,help="Path to a directory containing mutant JSON reports to integrate.",)
    parser.add_argument("-wo","--wt_output_dir", required=True,help="Path to wild-type output directory for mut_integrated JSON files.", )
    parser.add_argument("-mo","--mut_output_dir", required=True,help="Path to mutant output directory for mut_integrated JSON files.", )
    parser.add_argument("--strict", dest="strict", action="store_true",help="Enable strict mode requiring all 12 report types and all node fields (default: Disabled).")
    parser.set_defaults(strict=False)
    parser.set_defaults(func=run_mut_integrate)


def run_mut_integrate(args: Namespace) -> None:
    run_mut_integrate_service(
        mutclean_report_path=args.mut_clean_report_path,
        wt_input_dir=args.wt_input_dir,
        mut_input_dir=args.mut_input_dir,
        wt_output_dir=args.wt_output_dir,
        mut_output_dir=args.mut_output_dir,
        strict=args.strict,
    )

# ==============================
# Command: enzywizard-mut_integrate
# ==============================

# brief introduction:
'''
EnzyWizard-Mut-Integrate is a command-line tool for integrating wild-type and
mutant EnzyWizard JSON reports and constructing matched wild-type / mutant
protein graph representations.
It takes a enzywizard_mut_clean report, a wild-type report
directory and a mutant report directory as input, and merges information
from supported report types into structured paired graph datasets, where nodes
represent amino acids or substrates, and edges represent interactions.
The input enzywizard_mut_clean report is used as the anchor for cleaned amino acid 
substitution, residue indexing, and wild-type / mutant residue mapping.
Additional reports on each side provide complementary features.
The tool integrates the paired graph data, enabling direct downstream use in
mutation effect analysis, wild-type / mutant graph comparison, graph-based
machine learning, and enzyme engineering studies.
'''

# example usage:
'''
Example command:

enzywizard-mut_integrate -i examples/wt_input/mut_clean_report_1ZG4_WT_to_1ZG6_S70G.json -w examples/wt_input/ -m examples/mut_input/ -wo examples/wt_output/ -mo examples/mut_output/
'''

# input parameters:
'''
-i, --mut_clean_report_path
Required.
Path to the required enzywizard_mut_clean JSON report.

This report provides:
- amino acid substitution information
- cleaned amino acid substitution information
- wild-type residue mapping after cleaning
- mutant residue mapping after cleaning
- wild-type clean statistics
- mutant clean statistics

-w, --wt_input_dir
Required.
Path to a directory containing wild-type JSON reports to integrate.

-m, --mut_input_dir
Required.
Path to a directory containing mutant JSON reports to integrate.

Supported side report types include:
- enzywizard_aaprops
- enzywizard_hydrocluster
- enzywizard_energy
- enzywizard_flexibility
- enzywizard_disorder
- enzywizard_conservation
- enzywizard_embedding
- enzywizard_pocket
- enzywizard_substrate
- enzywizard_dock
- enzywizard_interaction

Duplicate report types are not allowed on either side.

-wo, --wt_output_dir
Required.
Directory to save wild-type mut-integrated JSON outputs.

-mo, --mut_output_dir
Required.
Directory to save mutant mut-integrated JSON outputs.

wt_output_dir and mut_output_dir must be different directories.

--strict
Optional.
Enable strict mode requiring all 12 report types and all node fields.

'''

# output content:
'''
The program outputs the following files into the output directories:

1. A shared mut-integrate JSON report
   - mut_integrate_report_{wt_protein_name}_to_{mut_protein_name}.json

   This file is saved into both wt_output_dir and mut_output_dir.

   The JSON report contains:

   - "output_type"
     A string identifying the report type:
     "enzywizard_mut_integrate"

   - "cleaned_amino_acid_substitution"
     Mutation definition after cleaning.

   - "overall_statistics"
     A dictionary summarizing paired wild-type / mutant overall features.

     It may contain:
     - wild-type statistics with prefix "wt_"
     - mutant statistics with prefix "mut_"
     - difference values with prefix "diff_"

   - "mutation_site_features"
     A dictionary summarizing mutation-site features extracted from
     wild-type and mutant reports.

     It may contain:
     - wild-type mutation-site features with prefix "wt_"
     - mutant mutation-site features with prefix "mut_"
     - mutation differences with prefix "diff_"

   - "wt_integrated_graph"
     Integrated graph entries for the wild-type protein.

   - "mut_integrated_graph"
     Integrated graph entries for the mutant protein.

2. Wild-type node-only JSON file
   - wt_integrate_nodes_{wt_protein_name}.json

   Contains all wild-type node records extracted from wt_integrated_graph.

3. Wild-type edge-only JSON file
   - wt_integrate_edges_{wt_protein_name}.json

   Contains all wild-type edge records extracted from wt_integrated_graph.

4. Mutant node-only JSON file
   - mut_integrate_nodes_{mut_protein_name}.json

   Contains all mutant node records extracted from mut_integrated_graph.

5. Mutant edge-only JSON file
   - mut_integrate_edges_{mut_protein_name}.json

   Contains all mutant edge records extracted from mut_integrated_graph.
'''

# Process:
'''
This command processes the mutation-paired JSON reports as follows:

1. Validate input paths
   - Check that mut_clean_report_path exists.
   - Check that wt_input_dir and mut_input_dir are valid directories.
   - Check that wt_output_dir and mut_output_dir are different directories.
   - Create output directories if needed.

2. Resolve protein names
   - Parse wild-type and mutant protein names from the mut_clean report filename.

3. Load and validate mut_clean report
   - Read the enzywizard_mut_clean JSON report.
   - Validate its schema and required mutation fields.
   - Confirm that output_type is enzywizard_mut_clean.

4. List side JSON files
   - Search wt_input_dir for wild-type JSON files.
   - Search mut_input_dir for mutant JSON files.

5. Validate report counts
   - Ensure the number of JSON files on each side does not exceed the maximum number of supported side report types.
   - In strict mode, require exactly all supported side report types on both sides.

6. Load and validate wild-type reports
   - Read each JSON report from wt_input_dir.
   - Validate each report using output_type-specific schema checks.
   - Reject unsupported report types for mut_integrate.
   - Reject duplicated report types on the wild-type side.

7. Load and validate mutant reports
   - Read each JSON report from mut_input_dir.
   - Validate each report using output_type-specific schema checks.
   - Reject unsupported report types for mut_integrate.
   - Reject duplicated report types on the mutant side.

8. Synthesize clean reports
   - Generate a synthetic enzywizard_clean report for the wild-type side from mut_clean_report.
   - Generate a synthetic enzywizard_clean report for the mutant side from mut_clean_report.

9. Build full report dictionaries
   - Add synthesized clean reports into the wild-type and mutant report dictionaries.
   - Use these synthesized clean reports as anchor reports for residue indexing.

10. Validate cleaned mutation definition
   - Read cleaned_amino_acid_substitution from mut_clean_report.
   - Extract cleaned residue lists for wild-type and mutant proteins.
   - Validate that the cleaned mutation definition is compatible with the cleaned residue mappings.

11. Build overall statistics
   - Construct paired wild-type / mutant overall statistics.
   - Reorder fields into stable WT / MUT / DIFF output structure.

12. Build mutation-site features
   - Parse mutation positions from cleaned_amino_acid_substitution.
   - Match mutation residues in both cleaned residue lists.
   - Extract mutation-site features from available reports, including:
     - amino acid properties
     - flexibility
     - conservation
   - Compute WT / MUT / DIFF mutation-site summaries.

13. Build paired integrated graphs
   - Construct wild-type integrated graph from wild-type reports.
   - Construct mutant integrated graph from mutant reports.
   - If interaction reports are missing on either side, return empty integrated graphs.

14. Generate mut-integrate report
   - Assemble cleaned mutation definition, overall statistics, mutation-site features,
     wild-type integrated graph, and mutant integrated graph into one report.

15. Save mut-integrate report
   - Write the same mut-integrate report into wt_output_dir.
   - Write the same mut-integrate report into mut_output_dir.

16. Split integrated graphs
   - Parse wt_integrated_graph into:
     - wild-type node list
     - wild-type edge list
   - Parse mut_integrated_graph into:
     - mutant node list
     - mutant edge list

17. Save node and edge outputs
   - Save wild-type node and edge JSON files.
   - Save mutant node and edge JSON files.

18. Finalize outputs
   - Finish mut_integrate processing.
   - Copy the generated log file from wt_output_dir to mut_output_dir.
'''

# dependencies:
'''
- Biopython
- NumPy
- JSON
'''

# references:
'''
- Biopython:
  https://biopython.org/

- JSON:
  https://www.json.org/
'''
