"""Core computation of the TransQuant probe weight factor W.

Reimplemented from the MATLAB functions ``compute_correction_factor.m``,
``output_UCSC_probe_track.m``, ``find_probe_fit.m`` and
``compute_weight_factor.m`` shipped with Bahar Halpern & Itzkovitz,
Methods 2016 (doi:10.1016/j.ymeth.2015.11.015), Eq. 7 (there called eta).

Definitions
-----------
L      length of the target (full genomic locus, in transcription direction)
mid    midpoint of a probe binding site on the target (1-based)
N(i)   number of binding sites with mid < i, for i = 1..L
W      (1/L) * sum_i N(i) / N(L)

No Streamlit imports here: this module is shared by the CLI and the app.
"""

from __future__ import annotations

import io
import re
from collections import Counter
from dataclasses import dataclass, field

import numpy as np

_RC_TABLE = str.maketrans("ACGTN", "TGCAN")
_TARGET_OK = re.compile(r"^[ACGTN]*$")
_PROBE_OK = re.compile(r"^[ACGTU]{8,}$")
_SPLIT = re.compile(r"[\t,;]+|\s{2,}|\s+")


class InputError(ValueError):
    """User input cannot be interpreted; message is safe to show verbatim."""


@dataclass(frozen=True)
class Hit:
    probe: str      # probe sequence as supplied (antisense, normalised)
    start: int      # 1-based, inclusive, on the target
    end: int        # 1-based, inclusive
    mid: float


@dataclass
class Result:
    W: float
    L: int
    N: np.ndarray
    hits: list[Hit]
    n_probes: int
    unmatched: list[str]
    multi_site: dict[str, int]
    warnings: list[str] = field(default_factory=list)

    @property
    def n_sites(self) -> int:
        return len(self.hits)

    @property
    def n_matched_probes(self) -> int:
        return self.n_probes - len(self.unmatched)


# ---------------------------------------------------------------- sequences

def revcomp(seq: str) -> str:
    """Reverse complement of an A/C/G/T/N string (upper case)."""
    return seq.translate(_RC_TABLE)[::-1]


def _normalise(seq: str) -> str:
    seq = re.sub(r"[\s\d]", "", seq).upper()
    return seq.replace("U", "T")


def _extract(text: str, concatenate: bool) -> tuple[str, list[str], list[str]]:
    """Return (sequence with original letter case, FASTA headers, warnings)."""
    warnings: list[str] = []
    headers: list[str] = []
    text = text.strip()
    if text.startswith(">"):
        entries = [e for e in re.split(r"^>", text, flags=re.M) if e.strip()]
        seqs = []
        for e in entries:
            lines = e.splitlines()
            headers.append(lines[0])
            seqs.append(re.sub(r"[\s\d]", "", "".join(lines[1:])))
        seqs = [q for q in seqs if q]
        if len(seqs) > 1:
            if concatenate:
                warnings.append(
                    f"FASTA contained {len(seqs)} entries; they were concatenated in order."
                )
                seq = "".join(seqs)
            else:
                warnings.append(
                    f"FASTA contained {len(seqs)} entries; only the first was used."
                )
                seq = seqs[0]
        else:
            seq = seqs[0] if seqs else ""
    else:
        seq = re.sub(r"[\s\d]", "", text)
    return seq.replace("U", "T").replace("u", "t"), headers, warnings


def parse_target(text: str, concatenate: bool = True) -> tuple[str, list[str]]:
    """Return (sequence, warnings) from pasted text or FASTA content.

    Whitespace and digits are stripped, U -> T, upper case. Allowed: A C G T N.
    Multi-entry FASTA is concatenated in order (UCSC exon/intron download)
    unless ``concatenate`` is False, in which case only the first entry is used.
    """
    cased, _, warnings = _extract(text, concatenate)
    seq = cased.upper()
    if not seq:
        raise InputError("Target sequence is empty.")
    if not _TARGET_OK.match(seq):
        bad = re.search(r"[^ACGTN]", seq)
        raise InputError(
            f"Target contains an unexpected character {bad.group()!r} at position {bad.start() + 1}. "
            "Only A, C, G, T, U and N are allowed."
        )
    return seq, warnings


def case_blocks(text: str, concatenate: bool = True) -> list[tuple[int, int]] | None:
    """1-based inclusive (start, end) runs of UPPERCASE letters in the target, or None.

    UCSC's Genomic Sequence download writes exons in upper case and introns in
    lower case, so these runs usually mark exons (or CDS/UTR, depending on the
    user's download options). Returns None when the sequence is single-case or
    when a FASTA header says ``repeatMasking=lower`` (lower case then means
    repeats, not introns).
    """
    cased, headers, _ = _extract(text, concatenate)
    if any("repeatMasking=lower" in h for h in headers):
        return None
    letters = re.sub(r"[^A-Za-z]", "", cased)
    if not letters or letters.isupper() or letters.islower():
        return None
    return [(m.start() + 1, m.end()) for m in re.finditer(r"[A-Z]+", cased)]


def parse_probes(text: str) -> tuple[list[str], list[str]]:
    """Return (probes, warnings) from one-per-line text or a spreadsheet paste.

    Rows are split on tabs, commas, semicolons or whitespace. Rows with no
    sequence-like cell (headers, blanks, notes) are dropped. Exactly one column
    must be sequence-like in every remaining row. Duplicates are collapsed.
    """
    warnings: list[str] = []
    rows = [_SPLIT.split(line.strip()) for line in text.splitlines() if line.strip()]
    rows = [[_normalise_cell(c) for c in r] for r in rows]

    def is_seq(cell: str) -> bool:
        return bool(_PROBE_OK.match(cell))

    kept = [r for r in rows if any(is_seq(c) for c in r)]
    dropped = len(rows) - len(kept)
    if dropped:
        warnings.append(
            f"{dropped} non-sequence row{'s' if dropped != 1 else ''} (header/blank) removed."
        )
    if not kept:
        raise InputError("No probe sequences found. Expected one probe per row (A/C/G/T/U only).")

    width = max(len(r) for r in kept)
    seq_cols = [
        c for c in range(width)
        if all(c < len(r) and is_seq(r[c]) for r in kept)
    ]
    if len(seq_cols) != 1:
        if not seq_cols:
            raise InputError(
                "Could not find a single column where every row is a probe sequence. "
                "Check for mixed content or missing cells."
            )
        raise InputError(
            f"{len(seq_cols)} columns look like sequences; paste only the probe sequence column."
        )

    probes = [r[seq_cols[0]].replace("U", "T") for r in kept]
    counts = Counter(probes)
    dups = sum(v - 1 for v in counts.values())
    if dups:
        warnings.append(f"{dups} duplicate probe sequence{'s' if dups != 1 else ''} collapsed.")
        probes = list(dict.fromkeys(probes))
    return probes, warnings


def _normalise_cell(cell: str) -> str:
    return cell.strip().upper()


# ----------------------------------------------------------------- matching

def locate_probes(target: str, probes: list[str], revcomp_probes: bool = True) -> list[Hit]:
    """Exact-match every probe (reverse complemented by default) on the target.

    Every occurrence is returned, including overlapping ones. Probes are
    expected to be antisense oligos as ordered; set ``revcomp_probes=False``
    if they are already target-sense.
    """
    hits: list[Hit] = []
    for probe in probes:
        query = revcomp(probe) if revcomp_probes else probe
        pos = target.find(query)
        while pos != -1:
            start, end = pos + 1, pos + len(query)
            hits.append(Hit(probe, start, end, (start + end) / 2))
            pos = target.find(query, pos + 1)
    return hits


# ---------------------------------------------------------------- W profile

def probe_profile(L: int, mids: list[float]) -> np.ndarray:
    """N(i) for i = 1..L: number of binding-site midpoints strictly below i."""
    sorted_mids = np.sort(np.asarray(mids, dtype=float))
    positions = np.arange(1, L + 1, dtype=float)
    return np.searchsorted(sorted_mids, positions, side="left")


def weight_factor(N: np.ndarray, L: int) -> float:
    """W = (1/L) * sum(N) / N(L)."""
    if N[-1] == 0:
        raise InputError("W is undefined: no probe binding sites on the target.")
    return float(N.sum() / (L * N[-1]))


def compute(
    target: str,
    probes: list[str],
    revcomp_probes: bool = True,
    warnings: list[str] | None = None,
) -> Result:
    """Full pipeline on already-parsed inputs. Warnings from parsing can be passed through."""
    warnings = list(warnings or [])
    hits = locate_probes(target, probes, revcomp_probes)
    if not hits:
        other = locate_probes(target, probes, not revcomp_probes)
        if other:
            hint = (
                "No probe matched, but probes do match without reverse complementing. "
                "They look like target-sense sequences: untick 'reverse complement probes' "
                "(CLI: --no-revcomp)."
                if revcomp_probes else
                "No probe matched, but probes do match after reverse complementing. "
                "They look like antisense probes: tick 'reverse complement probes' "
                "(CLI: drop --no-revcomp)."
            )
            raise InputError(hint)
        raise InputError(
            "No probe matched the target. Check that the target is the genomic locus in the "
            "direction of transcription and that the probes belong to this gene."
        )

    matched = {h.probe for h in hits}
    unmatched = [p for p in probes if p not in matched]
    site_counts = Counter(h.probe for h in hits)
    multi_site = {p: n for p, n in site_counts.items() if n > 1}

    if unmatched:
        warnings.append(
            f"{len(unmatched)} probe{'s' if len(unmatched) != 1 else ''} not found in the target "
            "and excluded from W: " + ", ".join(unmatched)
        )
    if multi_site:
        desc = "; ".join(
            f"{p} x{n} at " + ", ".join(str(h.start) for h in hits if h.probe == p)
            for p, n in multi_site.items()
        )
        warnings.append(
            f"{len(multi_site)} probe{'s' if len(multi_site) != 1 else ''} bind at several sites; "
            "every site is counted (" + desc + ")."
        )

    L = len(target)
    N = probe_profile(L, [h.mid for h in hits])
    W = weight_factor(N, L)
    return Result(W, L, N, hits, len(probes), unmatched, multi_site, warnings)


# ------------------------------------------------------------------- plot

def make_plot(res: Result, blocks: list[tuple[int, int]] | None = None, max_points: int = 50_000):
    """Probe map (top) and localisation profile (bottom, as in Fig. 3C of the paper).

    ``blocks`` are (start, end) runs to shade on the RNA bar (see ``case_blocks``).
    Returns a matplotlib Figure.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    main, accent = "#bb5a38", "#3d3a2a"
    multi_colour = "#2a6f97"

    fig, (ax_map, ax) = plt.subplots(
        2, 1, figsize=(8, 5.4), dpi=120, sharex=True,
        height_ratios=[1, 4], layout="constrained",
    )
    fig.get_layout_engine().set(hspace=0.02)
    fig.suptitle(f"Gene length L = {res.L:,} bp    W = {res.W:.4f}")

    # -- top: target RNA bar with probe barcode --------------------------------
    bar_y, bar_h = 0.0, 0.28
    ax_map.add_patch(Rectangle((1, bar_y - bar_h / 2), res.L, bar_h, color="#d3d2ca", lw=0))
    for b0, b1 in blocks or []:
        ax_map.add_patch(Rectangle((b0, bar_y - bar_h / 2), b1 - b0 + 1, bar_h, color="#8f8b7a", lw=0))
    single = [h.mid for h in res.hits if h.probe not in res.multi_site]
    multi = [h.mid for h in res.hits if h.probe in res.multi_site]
    ax_map.vlines(single, 0.28, 0.95, color=main, lw=1.0, alpha=0.85)
    if multi:
        ax_map.vlines(multi, 0.28, 0.95, color=multi_colour, lw=1.0, alpha=0.85)
    ax_map.text(1, -0.55, "5'", ha="left", va="top", fontsize=10, color=accent)
    ax_map.text(res.L, -0.55, "3'", ha="right", va="top", fontsize=10, color=accent)
    ax_map.set_ylim(-1.0, 1.05)
    ax_map.set_yticks([])
    for sp in ax_map.spines.values():
        sp.set_visible(False)
    ax_map.tick_params(axis="x", length=0)

    # -- bottom: N(i) step plot with a marker per binding site ----------------
    x = np.arange(1, res.L + 1)
    y = res.N
    if res.L > max_points:
        idx = np.unique(np.concatenate([
            np.linspace(0, res.L - 1, max_points).astype(int),
            np.clip(np.array([int(h.mid) for h in res.hits]) - 1, 0, res.L - 1),
            np.clip(np.array([int(h.mid) for h in res.hits]), 0, res.L - 1),
        ]))
        x, y = x[idx], y[idx]
    ax.step(x, y, where="post", color=main, linewidth=1.4)

    order = sorted(res.hits, key=lambda h: h.mid)
    mids = np.array([h.mid for h in order])
    ranks = np.arange(1, len(order) + 1)
    is_multi = np.array([h.probe in res.multi_site for h in order])
    ax.plot(mids[~is_multi], ranks[~is_multi], "o", ms=3.6, mfc="white", mec=main, mew=1.2, label="probe binding site")
    if is_multi.any():
        ax.plot(mids[is_multi], ranks[is_multi], "o", ms=3.6, mfc="white", mec=multi_colour, mew=1.2,
                label="probe binding at several sites")
        ax.legend(loc="upper left", frameon=False, fontsize=9)

    ax.set_xlabel("Position along the gene (bp)")
    ax.set_ylabel("Probes bound to nascent RNA")
    ax.set_xlim(1, res.L)
    ax.set_ylim(0, res.n_sites * 1.05)
    ax2 = ax.twinx()
    ax2.set_ylim(0, 105)
    ax2.set_ylabel("% of full probe set")
    ax.spines["top"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    return fig


def figure_bytes(fig, fmt: str = "png") -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt, bbox_inches="tight")
    return buf.getvalue()
