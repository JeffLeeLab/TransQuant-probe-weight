"""Cross-check against the earlier R translation, and a CLI smoke test.

tests/data/golden_* is a synthetic 12,345 bp locus with 48 unique 20-nt probes.
Expected W was produced by archive/jeff_r_implementation/calculate_probe_weight_factor.R
on the same files (R_W = 0.529840). Tolerance absorbs the +-1 base threshold difference.
"""

from pathlib import Path

import pytest

from transquant_w import cli, core

DATA = Path(__file__).parent / "data"
R_W = 0.529840


def test_golden_matches_r_translation():
    target, _ = core.parse_target((DATA / "golden_target.fa").read_text())
    probes, _ = core.parse_probes((DATA / "golden_probes.txt").read_text())
    res = core.compute(target, probes)
    assert res.L == 12345
    assert res.n_sites == 48
    assert res.W == pytest.approx(R_W, abs=1e-3)


def test_cli_prints_results_and_writes_plot(tmp_path, capsys):
    out = tmp_path / "gene"
    rc = cli.main([
        "--target", str(DATA / "golden_target.fa"),
        "--probes", str(DATA / "golden_probes.txt"),
        "--out", str(out), "--svg", "--profile-tsv",
    ])
    assert rc == 0
    stdout = capsys.readouterr().out
    assert "W\t0.5298" in stdout
    assert "L\t12345" in stdout
    assert (tmp_path / "gene.png").read_bytes()[:4] == b"\x89PNG"
    assert (tmp_path / "gene.svg").exists()
    assert sum(1 for _ in open(tmp_path / "gene_profile.tsv")) == 12346


def test_cli_input_error_exit_code(tmp_path, capsys):
    bad = tmp_path / "bad.txt"
    bad.write_text("not a sequence\n")
    rc = cli.main(["--target", str(DATA / "golden_target.fa"), "--probes", str(bad)])
    assert rc == 1
    assert "error:" in capsys.readouterr().err
