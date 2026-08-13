import gzip

import pandas as pd
import pytest

from simspace import SimSpace, cli, util


def _small_parameters():
    return util.generate_random_parameters(
        n_group=2,
        n_state=3,
        density_min=1.0,
        density_max=1.0,
        seed=7,
    )


def _common_small_arguments(output_dir):
    return [
        "--output-dir", str(output_dir),
        "--shape", "8", "8",
        "--neighbor-distance", "1",
        "--cell-sweeps", "1",
        "--niche-sweeps", "1",
        "--n-genes", "6",
        "--bg-ratio", "0.5",
        "--lr-ratio", "0",
        "--seed", "7",
    ]


def _assert_output_contract(output_dir):
    for name in cli.OUTPUT_FILES:
        assert (output_dir / name).is_file()

    metadata = pd.read_csv(output_dir / "metadata.csv")
    with gzip.open(output_dir / "molecular_profiles.csv.gz", "rt") as stream:
        profiles = pd.read_csv(stream)

    assert list(metadata.columns[:5]) == ["cell_id", "x", "y", "niche", "cell_type"]
    assert metadata["cell_id"].tolist() == profiles["cell_id"].tolist()
    assert (output_dir / "spatial_scatter.png").stat().st_size > 0
    return metadata, profiles


@pytest.mark.parametrize(
    ("command", "expected_sections", "expected_phrases"),
    (
        (
            "reference-free",
            ("spatial simulation options:", "reference-free design options:", "native molecular-profile options:"),
            ("Number of spatial niches", "cell-retention probability", "Gamma shape and scale"),
        ),
        (
            "reference-based",
            (
                "reference data and parameter source:",
                "spatial fitting options (used with --fit-spatial):",
                "reference-based molecular-profile options:",
            ),
            ("cell IDs in the first column", "evolutionary-algorithm generations", "local R environment"),
        ),
    ),
)
def test_subcommand_help_explains_parameters(command, expected_sections, expected_phrases, capsys):
    with pytest.raises(SystemExit) as exit_info:
        cli.main([command, "--help"])

    assert exit_info.value.code == 0
    help_text = capsys.readouterr().out
    for text in (*expected_sections, *expected_phrases):
        assert text in help_text
    assert "(default:" in help_text


def test_reference_free_cli_writes_documented_outputs(tmp_path):
    output_dir = tmp_path / "reference_free"
    arguments = [
        "reference-free",
        *_common_small_arguments(output_dir),
        "--n-niches", "2",
        "--n-cell-types", "3",
        "--density-min", "1",
        "--density-max", "1",
    ]

    assert cli.main(arguments) == 0
    metadata, profiles = _assert_output_contract(output_dir)
    assert len(metadata) == 64
    assert profiles.shape == (64, 7)
    assert metadata["cell_type"].str.startswith("CellType_").all()


def test_reference_based_cli_maps_reference_labels(tmp_path):
    reference_metadata = pd.DataFrame(
        {
            "cell_type": ["A"] * 6 + ["B"] * 4 + ["C"] * 2,
            "x": list(range(12)),
            "y": list(reversed(range(12))),
        },
        index=[f"reference_{index}" for index in range(12)],
    )
    reference_path = tmp_path / "reference_metadata.csv"
    reference_metadata.to_csv(reference_path)
    parameter_path = tmp_path / "parameters.json"
    util.save_params(_small_parameters(), parameter_path)
    output_dir = tmp_path / "reference_based"

    arguments = [
        "reference-based",
        *_common_small_arguments(output_dir),
        "--reference-metadata", str(reference_path),
        "--cell-type-column", "cell_type",
        "--x-column", "x",
        "--y-column", "y",
        "--params", str(parameter_path),
        "--profile-model", "native",
    ]

    assert cli.main(arguments) == 0
    metadata, profiles = _assert_output_contract(output_dir)
    assert set(metadata["cell_type"]) <= {"A", "B", "C"}
    assert profiles.shape[1] == 7


def test_reference_based_fit_option_uses_existing_optimizer(tmp_path, monkeypatch):
    reference_metadata = pd.DataFrame(
        {
            "cell_type": ["A"] * 6 + ["B"] * 4 + ["C"] * 2,
            "x": list(range(12)),
            "y": list(reversed(range(12))),
        }
    )
    reference_path = tmp_path / "reference_metadata.csv"
    reference_metadata.to_csv(reference_path)
    output_dir = tmp_path / "fitted"
    calls = []

    def fitted_parameters(metadata, args):
        calls.append((len(metadata), args.n_niches, args.population_size, args.generations))
        return _small_parameters()

    monkeypatch.setattr(cli, "_fit_reference_parameters", fitted_parameters)
    arguments = [
        "reference-based",
        *_common_small_arguments(output_dir),
        "--reference-metadata", str(reference_path),
        "--cell-type-column", "cell_type",
        "--x-column", "x",
        "--y-column", "y",
        "--fit-spatial",
        "--n-niches", "2",
        "--population-size", "4",
        "--generations", "1",
    ]

    assert cli.main(arguments) == 0
    _assert_output_contract(output_dir)
    assert calls == [(12, 2, 4, 1)]


@pytest.mark.parametrize(
    ("profile_model", "method_name"),
    (("scdesign3", "fit_scdesign"), ("srtsim", "fit_srtsim")),
)
def test_optional_profile_models_dispatch_without_running_r(tmp_path, monkeypatch, profile_model, method_name):
    reference_metadata = pd.DataFrame(
        {
            "cell_type": ["A"] * 6 + ["B"] * 4 + ["C"] * 2,
            "x": list(range(12)),
            "y": list(reversed(range(12))),
        },
        index=[f"reference_{index}" for index in range(12)],
    )
    reference_path = tmp_path / "reference_metadata.csv"
    reference_metadata.to_csv(reference_path)
    counts_path = tmp_path / "reference_counts.csv"
    pd.DataFrame({"Gene_A": [1] * 12}, index=reference_metadata.index).to_csv(counts_path)
    parameter_path = tmp_path / "parameters.json"
    util.save_params(_small_parameters(), parameter_path)
    output_dir = tmp_path / profile_model
    calls = []

    def fake_adapter(self, count_path, meta_path, group_col, spatial_x, spatial_y, seed=0, **kwargs):
        calls.append((count_path, meta_path, group_col, spatial_x, spatial_y, seed, kwargs))
        self.omics = pd.DataFrame({"ReferenceGene": [1] * len(self.meta)})

    monkeypatch.setattr(SimSpace, method_name, fake_adapter)
    arguments = [
        "reference-based",
        *_common_small_arguments(output_dir),
        "--reference-metadata", str(reference_path),
        "--reference-counts", str(counts_path),
        "--cell-type-column", "cell_type",
        "--x-column", "x",
        "--y-column", "y",
        "--params", str(parameter_path),
        "--profile-model", profile_model,
    ]

    assert cli.main(arguments) == 0
    _, profiles = _assert_output_contract(output_dir)
    assert list(profiles.columns) == ["cell_id", "ReferenceGene"]
    expected_options = {"family": "auto"} if profile_model == "scdesign3" else {}
    assert calls == [(str(counts_path), str(reference_path), "cell_type", "x", "y", 7, expected_options)]


def test_cli_refuses_to_replace_outputs_without_overwrite(tmp_path):
    output_dir = tmp_path / "existing"
    arguments = [
        "reference-free",
        *_common_small_arguments(output_dir),
        "--n-niches", "2",
        "--n-cell-types", "3",
        "--density-min", "1",
        "--density-max", "1",
    ]

    assert cli.main(arguments) == 0
    assert cli.main(arguments) == 1
    assert cli.main([*arguments, "--overwrite"]) == 0
