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