"""Command-line interface for headless SimSpace simulations."""

from __future__ import annotations

import argparse
import copy
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from . import __version__, spatial, util
from .optimize import spatial_fit


OUTPUT_FILES = (
    "metadata.csv",
    "molecular_profiles.csv.gz",
    "spatial_scatter.png",
    "parameters.json",
)


class _HelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    """Show defaults while preserving readable multi-line help text."""

    def _get_help_string(self, action: argparse.Action) -> str:
        help_text = action.help
        if (
            "%(default)" not in help_text
            and action.default is not argparse.SUPPRESS
            and action.default is not None
            and action.option_strings
        ):
            help_text += " (default: %(default)s)"
        return help_text


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def _unit_interval(value: str) -> float:
    parsed = float(value)
    if not 0 <= parsed <= 1:
        raise argparse.ArgumentTypeError("must be between 0 and 1")
    return parsed


def _nonnegative_float(value: str) -> float:
    parsed = float(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be non-negative")
    return parsed


def _add_common_spatial_arguments(parser: argparse.ArgumentParser) -> None:
    spatial_options = parser.add_argument_group("spatial simulation options")
    spatial_options.add_argument(
        "--shape",
        nargs=2,
        type=_positive_int,
        metavar=("ROWS", "COLS"),
        default=(100, 100),
        help="Size of the two-dimensional lattice before density sampling.",
    )
    spatial_options.add_argument(
        "--neighbor-distance",
        type=_positive_int,
        default=3,
        help="Grid radius used to define neighboring locations.",
    )
    spatial_options.add_argument(
        "--neighbor-method",
        choices=("manhattan", "euclidean"),
        default="manhattan",
        help="Distance rule used to construct the spatial neighborhood.",
    )
    spatial_options.add_argument(
        "--cell-sweeps",
        type=_positive_int,
        default=4,
        help="Number of Gibbs-sampling sweeps used to organize cell types.",
    )
    spatial_options.add_argument(
        "--niche-sweeps",
        type=_positive_int,
        default=6,
        help="Number of MRF sweeps used to organize spatial niches.",
    )
    spatial_options.add_argument(
        "--jitter",
        type=_nonnegative_float,
        default=0.2,
        help="Standard deviation of Gaussian noise added to cell coordinates; use 0 for an exact lattice.",
    )

    run_options = parser.add_argument_group("run and output options")
    run_options.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Non-negative random seed for reproducible simulation.",
    )
    run_options.add_argument(
        "--output-dir",
        "-o",
        type=Path,
        required=True,
        help="Directory for metadata, molecular profiles, parameters, and the scatterplot.",
    )
    run_options.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace the four SimSpace output files if they already exist.",
    )


def _add_omics_arguments(parser: argparse.ArgumentParser) -> None:
    omics_options = parser.add_argument_group("native molecular-profile options")
    omics_options.add_argument(
        "--n-genes",
        type=_positive_int,
        default=1000,
        help="Number of molecular features to generate with the native Gamma-Poisson model.",
    )
    omics_options.add_argument(
        "--bg-ratio",
        type=_unit_interval,
        default=0.2,
        help="Fraction of all features assigned as non-marker background features.",
    )
    omics_options.add_argument(
        "--lr-ratio",
        type=_unit_interval,
        default=0.5,
        help="Fraction of marker features assigned to ligand-receptor pairs.",
    )
    omics_options.add_argument(
        "--bg-param",
        nargs=2,
        type=float,
        metavar=("SHAPE", "SCALE"),
        default=(1.0, 1.0),
        help="Gamma shape and scale used for background-feature mean expression.",
    )
    omics_options.add_argument(
        "--marker-param",
        nargs=2,
        type=float,
        metavar=("SHAPE", "SCALE"),
        default=(5.0, 2.0),
        help="Gamma shape and scale used for marker-feature mean expression.",
    )
    omics_options.add_argument(
        "--paired-feature-effects",
        action="store_true",
        help="Apply the existing neighborhood-dependent effect to ligand-receptor feature pairs.",
    )
    omics_options.add_argument(
        "--omics-neighbors",
        type=_positive_int,
        default=20,
        help="Number of nearest cells used by paired-feature effects.",
    )
    omics_options.add_argument(
        "--spatial-effect",
        type=float,
        default=3.0,
        help="Expression multiplier used by paired-feature effects.",
    )
    omics_options.add_argument(
        "--spatial-effect-threshold",
        type=float,
        default=1.5,
        help="Neighbor-ratio threshold above which a paired-feature effect is increased.",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="simspace",
        description="Run reference-free or reference-based SimSpace simulations without a notebook.",
        formatter_class=_HelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"SimSpace {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    reference_free = subparsers.add_parser(
        "reference-free",
        help="Generate spatial organization and molecular profiles without a reference dataset.",
        description=(
            "Generate a two-dimensional spatial layout and native molecular profiles without a reference dataset. "
            "Supply --params to reuse a design, or use the random-parameter options below."
        ),
        formatter_class=_HelpFormatter,
    )
    _add_common_spatial_arguments(reference_free)
    _add_omics_arguments(reference_free)
    design_options = reference_free.add_argument_group("reference-free design options")
    design_options.add_argument(
        "--params",
        type=Path,
        help="SimSpace parameters.json file to reuse instead of generating random spatial parameters.",
    )
    design_options.add_argument(
        "--n-niches",
        type=_positive_int,
        default=3,
        help="Number of spatial niches when generating random parameters.",
    )
    design_options.add_argument(
        "--n-cell-types",
        type=_positive_int,
        default=8,
        help="Number of cell types when generating random parameters.",
    )
    design_options.add_argument(
        "--theta-max",
        type=_nonnegative_float,
        default=0.8,
        help="Maximum absolute random cell-type affinity within each niche.",
    )
    design_options.add_argument(
        "--niche-theta-max",
        type=_nonnegative_float,
        default=0.5,
        help="Maximum absolute random affinity between spatial niches.",
    )
    design_options.add_argument(
        "--density-min",
        type=_unit_interval,
        default=0.01,
        help="Minimum random cell-retention probability across cell types.",
    )
    design_options.add_argument(
        "--density-max",
        type=_unit_interval,
        default=0.4,
        help="Maximum random cell-retention probability across cell types.",
    )
    design_options.add_argument(
        "--phi-min",
        type=float,
        default=4.4,
        help="Minimum random MRF spatial-smoothness scale.",
    )
    design_options.add_argument(
        "--phi-max",
        type=float,
        default=5.0,
        help="Maximum random MRF spatial-smoothness scale.",
    )
    reference_free.set_defaults(handler=_run_reference_free)

    reference_based = subparsers.add_parser(
        "reference-based",
        help="Calibrate or reuse spatial parameters from a reference and generate a new realization.",
        description=(
            "Generate a new two-dimensional realization guided by reference cell metadata. "
            "Choose either --params for a fitted design or --fit-spatial to calibrate one."
        ),
        formatter_class=_HelpFormatter,
    )
    _add_common_spatial_arguments(reference_based)
    _add_omics_arguments(reference_based)
    reference_options = reference_based.add_argument_group("reference data and parameter source")
    reference_options.add_argument(
        "--reference-metadata",
        type=Path,
        required=True,
        help="CSV with cell IDs in the first column plus cell-type and numeric coordinate columns.",
    )
    reference_options.add_argument(
        "--reference-counts",
        type=Path,
        help="Reference count-matrix CSV; required only for scdesign3 or srtsim profiles.",
    )
    reference_options.add_argument(
        "--cell-type-column",
        required=True,
        help="Name of the cell-type label column in --reference-metadata.",
    )
    reference_options.add_argument(
        "--x-column",
        required=True,
        help="Name of the numeric horizontal-coordinate column in --reference-metadata.",
    )
    reference_options.add_argument(
        "--y-column",
        required=True,
        help="Name of the numeric vertical-coordinate column in --reference-metadata.",
    )
    parameter_source = reference_options.add_mutually_exclusive_group(required=True)
    parameter_source.add_argument(
        "--params",
        type=Path,
        help="Previously fitted SimSpace parameters.json file to use directly.",
    )
    parameter_source.add_argument(
        "--fit-spatial",
        action="store_true",
        help="Fit new spatial parameters from the reference metadata.",
    )

    fitting_options = reference_based.add_argument_group("spatial fitting options (used with --fit-spatial)")
    fitting_options.add_argument(
        "--n-niches",
        type=_positive_int,
        default=2,
        help="Number of niches to fit in the simulated spatial design.",
    )
    fitting_options.add_argument(
        "--population-size",
        type=_positive_int,
        default=50,
        help="Number of candidate parameter sets per evolutionary-algorithm generation.",
    )
    fitting_options.add_argument(
        "--generations",
        type=_positive_int,
        default=20,
        help="Number of evolutionary-algorithm generations.",
    )
    fitting_options.add_argument(
        "--mutation-rate",
        type=_unit_interval,
        default=0.2,
        help="Probability of mutating each candidate parameter set.",
    )
    fitting_options.add_argument(
        "--crossover-rate",
        type=_unit_interval,
        default=0.6,
        help="Probability of crossover between selected candidate parameter sets.",
    )
    fitting_options.add_argument(
        "--replicates",
        type=_positive_int,
        default=1,
        help="Simulation replicates used to evaluate each candidate's fitness.",
    )
    fitting_options.add_argument(
        "--parallel-fit",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Evaluate spatial-calibration candidates in parallel; use --no-parallel-fit to disable.",
    )

    profile_options = reference_based.add_argument_group("reference-based molecular-profile options")
    profile_options.add_argument(
        "--profile-model",
        choices=("native", "scdesign3", "srtsim"),
        default="native",
        help="Molecular-profile generator; scdesign3 and srtsim use packages installed in the local R environment.",
    )
    profile_options.add_argument(
        "--scdesign-family",
        choices=("auto", "nb", "gaussian"),
        default="auto",
        help=(
            "Marginal family for scDesign3 profiles. Auto uses negative binomial for "
            "integer-valued data and log-Gaussian for continuous nonnegative data."
        ),
    )
    reference_based.set_defaults(handler=_run_reference_based)
    return parser


def _validate_common_arguments(args: argparse.Namespace) -> None:
    if args.seed < 0:
        raise ValueError("--seed must be non-negative")
    if args.bg_param[0] <= 0 or args.bg_param[1] <= 0:
        raise ValueError("--bg-param values must be positive")
    if args.marker_param[0] <= 0 or args.marker_param[1] <= 0:
        raise ValueError("--marker-param values must be positive")


def _prepare_output_directory(path: Path, overwrite: bool) -> None:
    conflicts = [path / name for name in OUTPUT_FILES if (path / name).exists()]
    if conflicts and not overwrite:
        names = ", ".join(item.name for item in conflicts)
        raise FileExistsError(f"Output files already exist in {path}: {names}. Use --overwrite to replace them.")
    path.mkdir(parents=True, exist_ok=True)


def _simulate_from_parameters(args: argparse.Namespace, parameters: dict):
    neighbors = spatial.generate_offsets(args.neighbor_distance, args.neighbor_method)
    return util.sim_from_params(
        parameters=parameters,
        shape=tuple(args.shape),
        num_iteration=args.cell_sweeps,
        n_iter=args.niche_sweeps,
        custom_neighbor=neighbors,
        step=args.jitter,
        seed=args.seed,
    )


def _create_native_profiles(simulation, args: argparse.Namespace) -> None:
    simulation.create_omics(
        n_genes=args.n_genes,
        bg_ratio=args.bg_ratio,
        lr_ratio=args.lr_ratio,
        bg_param=tuple(args.bg_param),
        marker_param=tuple(args.marker_param),
        spatial=args.paired_feature_effects,
        k_neighors=args.omics_neighbors,
        spatial_effect=args.spatial_effect,
        se_threshold=args.spatial_effect_threshold,
    )


def _load_reference_metadata(args: argparse.Namespace) -> pd.DataFrame:
    if not args.reference_metadata.exists():
        raise FileNotFoundError(f"Reference metadata not found: {args.reference_metadata}")
    metadata = pd.read_csv(args.reference_metadata, index_col=0)
    required = (args.cell_type_column, args.x_column, args.y_column)
    missing = [column for column in required if column not in metadata.columns]
    if missing:
        raise ValueError(f"Reference metadata is missing required columns: {', '.join(missing)}")
    if metadata.empty:
        raise ValueError("Reference metadata is empty")
    if metadata[list(required)].isna().any().any():
        raise ValueError("Reference coordinate and cell-type columns cannot contain missing values")
    coordinate_columns = (args.x_column, args.y_column)
    nonnumeric = [column for column in coordinate_columns if not pd.api.types.is_numeric_dtype(metadata[column])]
    if nonnumeric:
        raise ValueError(f"Reference coordinate columns must be numeric: {', '.join(nonnumeric)}")
    return metadata


def _reference_target(metadata: pd.DataFrame, args: argparse.Namespace) -> np.ndarray:
    cell_types = metadata[args.cell_type_column].unique()
    morans_i = spatial.integrate_morans_I(
        metadata[args.cell_type_column],
        metadata[[args.x_column, args.y_column]],
        cell_types,
    )
    counts = metadata[args.cell_type_column].value_counts().sort_index()
    ordered_morans_i = np.asarray([morans_i[index] for index in np.argsort(-counts)])
    local_entropy = spatial.calculate_local_entropy(
        metadata[args.cell_type_column],
        metadata[[args.x_column, args.y_column]],
    )
    entropy_density, _ = np.histogram(local_entropy, bins=[value / 4 for value in range(12)], density=True)
    return np.hstack((entropy_density, ordered_morans_i))


def _fit_reference_parameters(metadata: pd.DataFrame, args: argparse.Namespace) -> dict:
    target = _reference_target(metadata, args)
    return spatial_fit(
        target=target,
        population_size=args.population_size,
        generations=args.generations,
        mutation_rate=args.mutation_rate,
        crossover_rate=args.crossover_rate,
        shape=tuple(args.shape),
        n_group=args.n_niches,
        n_state=metadata[args.cell_type_column].nunique(),
        custom_neighbor=spatial.generate_offsets(args.neighbor_distance, args.neighbor_method),
        num_iterations=args.cell_sweeps,
        n_iter=args.niche_sweeps,
        replicate=args.replicates,
        seed=args.seed,
        parallel=args.parallel_fit,
        verbose=True,
    )


def _map_reference_cell_types(simulation, metadata: pd.DataFrame, cell_type_column: str) -> None:
    ranked = {
        cell_type: rank
        for rank, cell_type in enumerate(metadata[cell_type_column].value_counts().index, start=1)
    }
    rank_to_cell_type = {rank: cell_type for cell_type, rank in ranked.items()}
    simulation.meta["fitted_celltype"] = simulation.meta["state_rank"].map(rank_to_cell_type)
    if simulation.meta["fitted_celltype"].isna().any():
        raise ValueError("The fitted parameters and reference cell-type labels are incompatible")


def _create_reference_profiles(simulation, metadata: pd.DataFrame, args: argparse.Namespace) -> None:
    _map_reference_cell_types(simulation, metadata, args.cell_type_column)
    if args.profile_model == "native":
        _create_native_profiles(simulation, args)
        return
    if args.reference_counts is None:
        raise ValueError(f"--reference-counts is required for --profile-model {args.profile_model}")
    if not args.reference_counts.exists():
        raise FileNotFoundError(f"Reference counts not found: {args.reference_counts}")

    if args.profile_model == "scdesign3":
        simulation.fit_scdesign(
            str(args.reference_counts),
            str(args.reference_metadata),
            args.cell_type_column,
            args.x_column,
            args.y_column,
            seed=args.seed,
            family=args.scdesign_family,
        )
    else:
        simulation.fit_srtsim(
            str(args.reference_counts),
            str(args.reference_metadata),
            args.cell_type_column,
            args.x_column,
            args.y_column,
            seed=args.seed,
        )


def _metadata_output(simulation) -> pd.DataFrame:
    source = simulation.meta.reset_index(drop=True)
    width = max(6, len(str(len(source))))
    cell_ids = [f"cell_{index:0{width}d}" for index in range(1, len(source) + 1)]
    if "fitted_celltype" in source.columns:
        cell_type = source["fitted_celltype"].astype(str)
    else:
        cell_type = source["state"].astype(int).map(lambda value: f"CellType_{value + 1}")

    return pd.DataFrame(
        {
            "cell_id": cell_ids,
            "x": source["col"].to_numpy(),
            "y": source["row"].to_numpy(),
            "niche": source["niche"].astype(int).to_numpy(),
            "cell_type": cell_type.to_numpy(),
        }
    )


def _write_outputs(simulation, parameters: dict, args: argparse.Namespace) -> None:
    metadata = _metadata_output(simulation)
    profiles = simulation.omics.reset_index(drop=True).copy()
    if len(metadata) != len(profiles):
        raise ValueError("Metadata and molecular profiles contain different numbers of cells")
    profiles.insert(0, "cell_id", metadata["cell_id"])

    metadata.to_csv(args.output_dir / "metadata.csv", index=False)
    profiles.to_csv(args.output_dir / "molecular_profiles.csv.gz", index=False, compression="gzip")
    util.save_params(copy.deepcopy(parameters), args.output_dir / "parameters.json")

    plt.switch_backend("Agg")
    figure, axis = plt.subplots(figsize=(7, 6), dpi=150)
    sns.scatterplot(data=metadata, x="x", y="y", hue="cell_type", s=12, linewidth=0, ax=axis)
    axis.set_aspect("equal")
    axis.set_title("SimSpace cell-type distribution")
    axis.legend(bbox_to_anchor=(1.02, 1), loc="upper left", title="Cell type")
    figure.tight_layout()
    figure.savefig(args.output_dir / "spatial_scatter.png", bbox_inches="tight")
    plt.close(figure)

    print(f"SimSpace {__version__} completed successfully.")
    print(f"Cells: {len(metadata):,}; features: {simulation.omics.shape[1]:,}")
    for name in OUTPUT_FILES:
        print(args.output_dir / name)


def _run_reference_free(args: argparse.Namespace) -> None:
    _validate_common_arguments(args)
    if args.density_min > args.density_max:
        raise ValueError("--density-min cannot exceed --density-max")
    if args.phi_min > args.phi_max:
        raise ValueError("--phi-min cannot exceed --phi-max")
    _prepare_output_directory(args.output_dir, args.overwrite)

    if args.params is not None:
        parameters = util.load_params(args.params)
    else:
        parameters = util.generate_random_parameters(
            n_group=args.n_niches,
            n_state=args.n_cell_types,
            theta=args.theta_max,
            niche_theta=args.niche_theta_max,
            density_min=args.density_min,
            density_max=args.density_max,
            phi_min=args.phi_min,
            phi_max=args.phi_max,
            seed=args.seed,
        )
    simulation = _simulate_from_parameters(args, parameters)
    _create_native_profiles(simulation, args)
    _write_outputs(simulation, parameters, args)


def _run_reference_based(args: argparse.Namespace) -> None:
    _validate_common_arguments(args)
    _prepare_output_directory(args.output_dir, args.overwrite)
    metadata = _load_reference_metadata(args)
    parameters = _fit_reference_parameters(metadata, args) if args.fit_spatial else util.load_params(args.params)
    if int(parameters["n_state"]) != metadata[args.cell_type_column].nunique():
        raise ValueError("The parameter file's cell-type count does not match the reference metadata")

    simulation = _simulate_from_parameters(args, parameters)
    _create_reference_profiles(simulation, metadata, args)
    _write_outputs(simulation, parameters, args)


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        args.handler(args)
    except Exception as error:
        print(f"simspace: error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
