from __future__ import annotations

import argparse

from .commands.mut_integrate import add_mut_integrate_parser


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="enzywizard-mut-integrate",
        description="EnzyWizard-Mut-Integrate: Integrate wild-type and mutant EnzyWizard JSON reports and constructing matched wild-type / mutant protein graph representations."
    )
    add_mut_integrate_parser(parser)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)