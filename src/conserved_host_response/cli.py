from __future__ import annotations

import argparse

from .audit import audit_manifest
from .construct import construct_framework
from .freeze import freeze_framework
from .score import score_framework
from .validate import validate_scores


def build_parser():
    parser = argparse.ArgumentParser(
        prog="chrf",
        description="Construct, freeze, score, audit, and validate conserved host-response frameworks.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    audit = sub.add_parser("audit", help="Audit dataset roles, dependence, and leakage risks")
    audit.add_argument("--manifest", required=True)
    audit.add_argument("--output-dir", required=True)

    construct = sub.add_parser("construct", help="Apply configuration-driven cross-system construction rules")
    construct.add_argument("--evidence", required=True)
    construct.add_argument("--config", required=True)
    construct.add_argument("--output-dir", required=True)

    freeze = sub.add_parser("freeze", help="Freeze gene, module, direction, and scoring definitions")
    freeze.add_argument("--gene-table", required=True)
    freeze.add_argument("--config", required=True)
    freeze.add_argument("--output-dir", required=True)

    score = sub.add_parser("score", help="Apply a frozen framework to an expression matrix")
    score.add_argument("--expression", required=True)
    score.add_argument("--gene-manifest", required=True)
    score.add_argument("--scoring-spec", required=True)
    score.add_argument("--output-dir", required=True)
    score.add_argument("--expression-gene-column", default="gene")

    validate = sub.add_parser("validate", help="Estimate post-freeze case-control or paired effects")
    validate.add_argument("--scores", required=True)
    validate.add_argument("--metadata", required=True)
    validate.add_argument("--output-dir", required=True)
    validate.add_argument("--group-column", required=True)
    validate.add_argument("--case", required=True)
    validate.add_argument("--control", required=True)
    validate.add_argument("--pair-column")
    validate.add_argument("--score-columns", help="Comma-separated score columns; default is all numeric scores")
    validate.add_argument("--filter-column")
    validate.add_argument("--filter-value")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.command == "audit":
        audit_manifest(args.manifest, args.output_dir)
    elif args.command == "construct":
        construct_framework(args.evidence, args.config, args.output_dir)
    elif args.command == "freeze":
        freeze_framework(args.gene_table, args.config, args.output_dir)
    elif args.command == "score":
        score_framework(
            args.expression,
            args.gene_manifest,
            args.scoring_spec,
            args.output_dir,
            args.expression_gene_column,
        )
    elif args.command == "validate":
        columns = args.score_columns.split(",") if args.score_columns else None
        validate_scores(
            args.scores,
            args.metadata,
            args.output_dir,
            args.group_column,
            args.case,
            args.control,
            args.pair_column,
            columns,
            args.filter_column,
            args.filter_value,
        )


if __name__ == "__main__":
    main()
