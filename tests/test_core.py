"""Analytic and behavioural tests for transquant_w.core (no MATLAB needed)."""

import random

import numpy as np
import pytest

from transquant_w import core

# ---------------------------------------------------------------- helpers ---

def random_dna(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def probes_at(target, starts, length=20):
    """Antisense probes whose binding sites start at the given 1-based positions."""
    return [core.revcomp(target[s - 1 : s - 1 + length]) for s in starts]


def assert_unique(target, probes):
    for p in probes:
        assert target.count(core.revcomp(p)) == 1


TARGET = random_dna(10_000, seed=1)

# ------------------------------------------------------- sequence handling ---

def test_revcomp_pure_python():
    assert core.revcomp("ACGTN") == "NACGT"
    assert core.revcomp("AAAC") == "GTTT"


def test_parse_target_plain_text_normalises():
    seq, warnings = core.parse_target("  acg\nuUg 12 34\r\n ")
    assert seq == "ACGTTG"


def test_parse_target_single_fasta():
    seq, warnings = core.parse_target(">chr1:1-10 strand=+\nACGT\nacgt\n")
    assert seq == "ACGTACGT"
    assert warnings == []


def test_parse_target_multi_fasta_concatenates_and_warns():
    text = ">exon1\nAAAA\n>intron1\nCCCC\n>exon2\nGGGG\n"
    seq, warnings = core.parse_target(text)
    assert seq == "AAAACCCCGGGG"
    assert any("3 entries" in w for w in warnings)


def test_parse_target_multi_fasta_no_concat_uses_first():
    text = ">a\nAAAA\n>b\nCCCC\n"
    seq, warnings = core.parse_target(text, concatenate=False)
    assert seq == "AAAA"
    assert any("first" in w.lower() for w in warnings)


def test_parse_target_allows_N_rejects_others():
    seq, _ = core.parse_target("ACGNNNT")
    assert seq == "ACGNNNT"
    with pytest.raises(core.InputError, match="position 4"):
        core.parse_target("ACGXT")


def test_parse_target_empty():
    with pytest.raises(core.InputError):
        core.parse_target(">only a header\n")


# ---------------------------------------------------------- probe parsing ---

def test_parse_probes_one_per_line():
    probes, warnings = core.parse_probes("acgtacgtacgtacgtacgt\nUUUUAAAACCCCGGGGTTTT\n")
    assert probes == ["ACGTACGTACGTACGTACGT", "TTTTAAAACCCCGGGGTTTT"]
    assert warnings == []


def test_parse_probes_excel_paste_with_header_and_extra_columns():
    text = (
        "Probe\tSequence\tGC%\tPosition\n"
        "1\tACGTACGTACGTACGTACGT\t50\t12\n"
        "2\tTTTTAAAACCCCGGGGTTTT\t45\t80\n"
    )
    probes, warnings = core.parse_probes(text)
    assert probes == ["ACGTACGTACGTACGTACGT", "TTTTAAAACCCCGGGGTTTT"]
    assert any("1 non-sequence row" in w for w in warnings)


def test_parse_probes_collapses_duplicates():
    probes, warnings = core.parse_probes("ACGTACGTACGTACGTACGT\nACGTACGTACGTACGTACGT\n")
    assert probes == ["ACGTACGTACGTACGTACGT"]
    assert any("duplicate" in w.lower() for w in warnings)


def test_parse_probes_ambiguous_columns_error():
    with pytest.raises(core.InputError, match="2 columns"):
        core.parse_probes("ACGTACGTACGTACGTACGT\tTTTTAAAACCCCGGGGTTTT\n")


def test_parse_probes_none_found():
    with pytest.raises(core.InputError):
        core.parse_probes("name\tgc\nfoo\t50\n")


# --------------------------------------------------------------- locating ---

def test_locate_probes_positions_and_midpoint():
    target = "T" * 100 + "ACGTACGTAC" + "T" * 100
    probe = core.revcomp("ACGTACGTAC")
    hits = core.locate_probes(target, [probe])
    assert len(hits) == 1
    h = hits[0]
    assert (h.start, h.end, h.mid) == (101, 110, 105.5)


def test_locate_probes_counts_every_site():
    site = "ACGTACGTAC"
    target = "T" * 50 + site + "G" * 50 + site + "T" * 50
    hits = core.locate_probes(target, [core.revcomp(site)])
    assert [h.start for h in hits] == [51, 111]


def test_locate_probes_never_spans_N():
    # Site would be GTAAC if the N were an A; an N in the target blocks the match.
    target = "ACGTNACGT"
    assert core.locate_probes(target, [core.revcomp("GTAAC")]) == []
    assert core.locate_probes("ACGTAACGT", [core.revcomp("GTAAC")]) != []


# ------------------------------------------------------ profile and W ------

def test_profile_and_W_hand_calculated():
    # L = 10, single site with mid 3.5: N(i) = 1 for i >= 4 -> sum = 7 -> W = 7/10
    N = core.probe_profile(10, [3.5])
    assert N.tolist() == [0, 0, 0, 1, 1, 1, 1, 1, 1, 1]
    assert core.weight_factor(N, 10) == pytest.approx(0.7)


def test_W_uniform_probes_is_half():
    starts = [491 + 1000 * k for k in range(10)]  # mids at 500.5, 1500.5, ...
    probes = probes_at(TARGET, starts)
    assert_unique(TARGET, probes)
    res = core.compute(TARGET, probes)
    assert res.W == pytest.approx(0.5, abs=0.01)
    assert res.n_sites == 10


def test_W_probes_at_5prime_near_one():
    probes = probes_at(TARGET, [10 + 40 * k for k in range(10)])
    assert_unique(TARGET, probes)
    assert core.compute(TARGET, probes).W > 0.9


def test_W_probes_at_3prime_near_zero():
    probes = probes_at(TARGET, [9_500 + 40 * k for k in range(10)])
    assert_unique(TARGET, probes)
    assert core.compute(TARGET, probes).W < 0.1


def test_compute_multi_site_probe_counted_twice():
    site = random_dna(20, seed=7)
    target = random_dna(2000, seed=8) + site + random_dna(2000, seed=9) + site + random_dna(2000, seed=10)
    assert target.count(site) == 2
    res = core.compute(target, [core.revcomp(site)])
    assert res.n_probes == 1
    assert res.n_sites == 2
    assert res.N[-1] == 2
    assert res.multi_site == {core.revcomp(site): 2}


def test_compute_unmatched_probe_reported_and_excluded():
    good = probes_at(TARGET, [1000, 5000])
    bad = "GGGGGGGGGGGGGGGGGGGG"
    res = core.compute(TARGET, good + [bad])
    assert res.unmatched == [bad]
    assert res.n_sites == 2
    assert any("1 probe" in w and "not found" in w for w in res.warnings)


def test_compute_zero_matches_gives_orientation_hint():
    sense = [TARGET[999:1019], TARGET[4999:5019]]  # target-sense, not antisense
    with pytest.raises(core.InputError, match="reverse complement"):
        core.compute(TARGET, sense)


def test_compute_revcomp_probes_false():
    sense = [TARGET[999:1019], TARGET[4999:5019]]
    res = core.compute(TARGET, sense, revcomp_probes=False)
    assert res.n_sites == 2


def test_compute_zero_matches_no_hint_when_nothing_matches():
    with pytest.raises(core.InputError, match="No probe"):
        core.compute(TARGET, ["GGGGGGGGGGGGGGGGGGGG"])


def test_profile_is_vectorised_for_long_genes():
    L = 2_000_000
    N = core.probe_profile(L, [float(x) for x in range(1000, L, 40_000)])
    assert N.shape == (L,)
    assert N[-1] == 50


def test_make_plot_returns_figure():
    probes = probes_at(TARGET, [1000, 5000, 9000])
    res = core.compute(TARGET, probes)
    fig = core.make_plot(res)
    assert fig.axes
    png = core.figure_bytes(fig, "png")
    assert png[:4] == b"\x89PNG"
