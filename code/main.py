from typing import Sequence
from pathlib import Path
import argparse


def add_preprocess_lv_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-o",
        "--mesh-folder",
        type=Path,
        default=Path.cwd().parent / "meshes-lv",
    )
    parser.add_argument("--r-short-endo", type=float, default=4.0)
    parser.add_argument("--r-long-endo", type=float, default=8.0)
    parser.add_argument("--r-short-epi", type=float, default=5.5)
    parser.add_argument("--r-long-epi", type=float, default=9.5)
    parser.add_argument("--psize-ref", type=float, default=0.5)


def add_preprocess_cylinder_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-o",
        "--mesh-folder",
        type=Path,
        default=Path.cwd().parent / "meshes-cylinder",
    )
    parser.add_argument("-c", "--char_length", type=float, default=50.0)
    parser.add_argument("-l", "--length", type=float, default=4000.0)
    parser.add_argument("-r", "--radius", type=float, default=300.0)


def add_run_lv_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-i",
        "--mesh-folder",
        type=Path,
        default=Path.cwd().parent / "meshes-lv",
    )
    parser.add_argument(
        "-o",
        "--output-folder",
        type=Path,
        default=Path.cwd().parent / "results-lv",
    )


def add_run_cylinder_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-i",
        "--mesh-folder",
        type=Path,
        default=Path.cwd().parent / "meshes-cylinder",
    )
    parser.add_argument(
        "-o",
        "--output-folder",
        type=Path,
        default=Path.cwd().parent / "results-cylinder",
    )
    parser.add_argument(
        "--varying-gamma",
        action="store_true",
        help="Vary gamma linearly from 0 in the center to amp in the border",
    )
    parser.add_argument(
        "--amp",
        type=float,
        default=0.2,
        help="Amplitude of gamma",
    )


def add_run_cylinder_twitch_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-i",
        "--mesh-folder",
        type=Path,
        default=Path.cwd().parent / "meshes-cylinder",
    )
    parser.add_argument(
        "-o",
        "--output-folder",
        type=Path,
        default=Path.cwd().parent / "results-cylinder",
    )
    parser.add_argument("--spring", type=float, default=3000.0, help="Spring constant")
    parser.add_argument(
        "--varying-gamma",
        action="store_true",
        help="Vary gamma linearly from 0 in the center to amp in the border",
    )
    parser.add_argument(
        "--amp",
        type=float,
        default=0.2,
        help="Amplitude of gamma",
    )


def add_postprocess_lv_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-r",
        "--resultsdir",
        type=Path,
        default=Path.cwd().parent / "results-lv",
    )
    parser.add_argument(
        "-i",
        "--mesh-folder",
        type=Path,
        default=Path.cwd().parent / "meshes-lv",
    )
    parser.add_argument(
        "-o",
        "--figdir",
        type=Path,
        default=Path.cwd().parent / "figures-lv",
    )
    parser.add_argument(
        "--print-stats",
        action="store_true",
        help="Print statistics",
    )
    parser.add_argument(
        "--create-paraview",
        action="store_true",
        help="Create Paraview files",
    )


def add_postprocess_cylinder_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-r",
        "--resultsdir",
        type=Path,
        default=Path.cwd().parent / "results-cylinder",
    )
    parser.add_argument(
        "-i",
        "--mesh-folder",
        type=Path,
        default=Path.cwd().parent / "meshes-cylinder",
    )
    parser.add_argument(
        "-o",
        "--figdir",
        type=Path,
        default=Path.cwd().parent / "figures-cylinder",
    )


def add_postprocess_cylinder_twitch_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-r",
        "--resultsdir",
        type=Path,
        default=Path.cwd().parent / "results-cylinder",
    )
    parser.add_argument(
        "-i",
        "--mesh-folder",
        type=Path,
        default=Path.cwd().parent / "meshes-cylinder",
    )
    parser.add_argument(
        "-o",
        "--figdir",
        type=Path,
        default=Path.cwd().parent / "figures-cylinder",
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    # Root parser
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Just print the command and do not run it",
    )

    subparsers = parser.add_subparsers(dest="command")

    # Subparsers

    # Preprocess LV
    preprocess_lv_parser = subparsers.add_parser("preprocess-lv", help="Create left ventricle mesh")
    add_preprocess_lv_arguments(preprocess_lv_parser)

    # Preprocess cylinder
    preprocess_cylinder_parser = subparsers.add_parser(
        "preprocess-cylinder", help="Create cylinder mesh"
    )
    add_preprocess_cylinder_arguments(preprocess_cylinder_parser)

    # Run LV
    lv_parser = subparsers.add_parser("run-lv", help="Run simulations with left ventricle model")
    add_run_lv_arguments(lv_parser)

    # Run cylinder
    cylinder_parser = subparsers.add_parser(
        "run-cylinder", help="Run simulations with cylinder model"
    )
    add_run_cylinder_arguments(cylinder_parser)

    # Run cylinder twitch
    cylinder_twitch_parser = subparsers.add_parser(
        "run-cylinder-twitch", help="Run simulations with cylinder model"
    )
    add_run_cylinder_twitch_arguments(cylinder_twitch_parser)

    # Postprocess LV
    postprocess_lv_parser = subparsers.add_parser("postprocess-lv", help="Postprocess LV results")
    add_postprocess_lv_arguments(postprocess_lv_parser)

    # Postprocess cylinder
    postprocess_cylinder_parser = subparsers.add_parser(
        "postprocess-cylinder", help="Postprocess cylinder results"
    )
    add_postprocess_cylinder_arguments(postprocess_cylinder_parser)

    # Postprocess cylinder twitch
    postprocess_cylinder_twitch_parser = subparsers.add_parser(
        "postprocess-cylinder-twitch", help="Postprocess cylinder results"
    )
    add_postprocess_cylinder_twitch_arguments(postprocess_cylinder_twitch_parser)

    args = vars(parser.parse_args(argv))

    if args.pop("dry_run"):
        import pprint

        pprint.pprint(args)
        return 0

    cmd = args.pop("command")

    if cmd == "preprocess-lv":
        from geometry import preprocess_lv

        preprocess_lv(**args)

    elif cmd == "preprocess-cylinder":
        from geometry import preprocess_cylinder

        preprocess_cylinder(**args)

    elif cmd == "run-lv":
        import lv

        lv.main(**args)

    elif cmd == "run-cylinder":
        import cylinder

        cylinder.main(**args)

    elif cmd == "run-cylinder-twitch":
        import cylinder

        cylinder.main_twitch(**args)

    elif cmd == "postprocess-lv":
        import postprocess_lv

        postprocess_lv.postprocess_lv(**args)

    elif cmd == "postprocess-cylinder":
        import postprocess_cylinder

        postprocess_cylinder.postprocess_cylinder(**args)

    elif cmd == "postprocess-cylinder-twitch":
        import postprocess_cylinder

        postprocess_cylinder.postprocess_cylinder_twitch(**args)

    else:
        parser.print_help()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
