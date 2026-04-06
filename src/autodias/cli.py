from __future__ import annotations

import argparse

from autodias.workflow import analyze, prepare, run


def main() -> int:
    parser = argparse.ArgumentParser(prog="autodias")
    subparsers = parser.add_subparsers(dest="command", required=True)

    for name in ("prepare", "run", "analyze"):
        command_parser = subparsers.add_parser(name)
        command_parser.add_argument("config", help="Path to a TOML configuration file")

    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args.config)
    elif args.command == "run":
        run(args.config)
    elif args.command == "analyze":
        analyze(args.config)
    return 0
