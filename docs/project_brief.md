# TransQuant probe weight factor (W) — project brief

## 1. Objective

Build a small, standalone tool that computes the TransQuant **probe weight factor W** for an smFISH probe library against its target gene, and shows the probe localisation plot. Nothing else from TransQuant is reimplemented.

Deliverables:

1. **Browser app** — Streamlit running fully client-side via `stlite` (Pyodide/WebAssembly), hosted on GitHub Pages and deployed by GitHub Actions (Pages "Actions" source, not the `docs/` folder method).
2. **CLI script** — same core code, for coding-friendly users.
3. **Faithful translation** of the algorithm MATLAB → R → Python, verified against analytic test cases and the existing R script.
4. **README** doubling as user manual (structure in §9).

Guiding rule: **lean, no overengineering**. This is a one-function app. Three Python source files, two runtime dependencies.

## 2. Background

- TransQuant is the smFISH quantification methodology and software from Bahar Halpern & Itzkovitz, *Methods* 2016 (https://doi.org/10.1016/j.ymeth.2015.11.015). Full text: `archive/halpern_itzkovitz_2016/Single molecule approaches ... .md`. Original MATLAB code: `archive/halpern_itzkovitz_2016/TransQuant_codes/`.
- The method infers transcription and degradation rates from smFISH images of transcription sites (TS) and single mature mRNAs. The TS intensity, relative to a single mRNA, gives the number of Pol II molecules on the gene. Because a nascent transcript attached to a Pol II halfway along the gene only carries the probes that bind upstream of that point, the TS intensity has to be corrected for **where the probes sit along the gene**. That correction is W (called *eta* in the paper's Eq. 7; we use **W** throughout, matching the code).
- W is close to 1 when probes cluster at the 5' end, ~0.5 when probes are spread uniformly, and small when probes cluster at the 3' end. The plot the app produces is the paper's Figure 3C.
- In our current workflow, spot detection and counting are done with `bigfish`. This app **only** computes W, L, N and the plot.

## 3. Algorithm (source of truth)

Two MATLAB functions define the algorithm:

- `compute_correction_factor.m` — wrapper: picks files, calls the probe-mapping step, then the weight-factor step.
- `output_UCSC_probe_track.m` + `find_probe_fit.m` — reverse-complements each probe and finds its position in the gene.
- `compute_weight_factor.m` — builds N and computes W.

Definitions:

- `L` — length of the target sequence in bases (full genomic locus, see §4).
- For each probe, `mid` = midpoint of its binding site on the target.
- `N(i)` for i = 1..L — number of probe binding sites whose midpoint lies before position i. This is the number of probes bound to a nascent RNA whose Pol II has reached position i.
- `W = (1/L) * sum(N(i) for i in 1..L) / N(L)`.

Previous R translation: `archive/jeff_r_implementation/calculate_probe_weight_factor.R`.

### 3.1 MATLAB → R comparison (done during brainstorming)

| Aspect | MATLAB | R | Decision for Python |
|---|---|---|---|
| W formula | `(1/L)*sum(N)/N(end)` | identical | identical (matches paper Eq. 7) |
| Probe orientation | reverse-complement probe, search target | same | same, pure Python (`str.maketrans` + `[::-1]`) |
| Matching | mismatch-tolerant best fit; every probe always placed | exact match; unmatched probe → NA, silently dropped | **exact match**; unmatched probes reported, not silently dropped |
| Probe length for midpoint | hardcoded 20 nt | actual probe length | **actual probe length** (R behaviour) |
| Strand handling | reads `strand=` from UCSC header and flips direction | none; assumes sequence is already in transcription direction | **assume transcription direction** and state it in the UI (§4) |
| N(i) threshold | `mid < i-1` (0-based flavour) | `mid < i` (1-based) | `mid < i`, 1-based; ±1 base difference is irrelevant for kb-scale genes |
| Multiple matches of one probe | first best hit | first hit | **all hits counted** (§5) |
| Complexity | loop over L, `find` over probes | loop over L, dplyr filter over probes → very slow for long genes | vectorised: sort midpoints once, `numpy.searchsorted` gives all N(i) at once |

Conclusion: the R translation is faithful in the formula, and improves on probe length. Its two weaknesses (silent NA drop, strand assumption not surfaced) are fixed by the decisions above.

## 4. Inputs

### 4.1 Target sequence

- Text box (paste) **or** FASTA upload. CLI: FASTA file path.
- **Stated assumption, displayed in an info box directly above the input:**
  - Provide the **full genomic locus**: all exons, introns and UTRs from transcription start to transcription end. Not the spliced mRNA/cDNA. Intron spacing is what makes W meaningful.
  - Provide it **in the direction of transcription** (5' → 3' of the RNA; i.e. the coding/sense strand). For a minus-strand gene in the UCSC Genome Browser, tick "reverse complement" when downloading DNA. The app does not read strand information from headers.
- Normalisation: strip whitespace, digits and line numbers; uppercase; `U → T`.
- Allowed characters after normalisation: `A C G T N`. `N` (undetermined regions) is fine; no probe can span an N so it never matches there. Any other character → error with position reported.
- Multi-entry FASTA: alert the user. Toggle "concatenate entries in order" (default **on**), because UCSC's exon/intron-separated download produces one entry per segment in genomic order. Concatenated sequence is what is used for L.
- Size cap on upload (e.g. 5 MB; longest human gene DMD ≈ 2.3 Mb). Note: stlite runs in the user's browser, so the cap protects their memory, not a server.

### 4.2 Probe sequences

- Text box (paste) **or** tab-delimited `.txt` upload. CLI: file path.
- Users typically paste straight from Excel/Google Sheets/Stellaris designer output, so rows may carry extra columns (probe number, name, GC%, position). **Column detection rule:** split each row on tabs/commas/whitespace; use the column in which every non-empty cell matches `^[ACGTU]+$` after uppercasing. If exactly one such column exists, use it; otherwise error with a clear message.
- Header rows (cells containing non-ACGTU characters) are removed and reported in the log.
- Normalisation: uppercase, `U → T`, strip whitespace.
- Duplicate probe sequences: collapse to one, report count in log.
- Probes are the **antisense** oligos as ordered (complementary to the RNA). The app reverse-complements them before searching. If **zero** probes match in that orientation but some match without reverse-complementing, tell the user they likely pasted target-sense sequences and let them flip with a toggle. Beyond this hint, no automatic orientation guessing.

## 5. Matching rules and edge cases

- Exact match only, after reverse complement.
- **A probe that matches several sites is counted at every site** (each site is a real binding event and adds fluorophores to the image). Report these probes and their positions in the log. Consequently `N(L)` = total binding sites, which can exceed the number of probes; the W formula is unchanged.
- Unmatched probes: excluded from N, listed by sequence in a warning. `N(L)` counts only matched sites.
- Zero matched probes: stop with an error (W undefined), show the orientation hint from §4.2.
- Probe longer than target, empty inputs, Windows line endings, lowercase, mixed-case UCSC repeat masking: all handled by normalisation or explicit errors.

## 6. Outputs

- **W** (4 decimal places), **L** (bp), number of probes supplied / matched / binding sites.
- **Plot**: N(i) against position i (step/line plot, not scatter). Y-axis: number of bound probes; optionally secondary axis as % of N(L). Title shows L and W. Downloadable as PNG (and SVG) via `st.download_button`. Downsample to ≤ ~50k points for genes > 50 kb.
- **N profile**: not displayed as a table and no CSV download (not useful; probe counts are ≤ ~200 and the plot carries the information). Available from the CLI as an optional TSV flag only.
- **Probe table**: sequence, start, end, midpoint, matched (yes/no/multiple). Small enough to display in full.
- **Log**: header rows removed, duplicates collapsed, unmatched probes, multi-site probes, FASTA concatenation notice.
- **"How to use W" box** — equations written for biologists, no Greek letters, one idea per line. Note that the paper applies a ceiling operator and a second factor *kappa* to Pol II occupancy; we present the simplified continuous form and point to the paper for the full version.

  ```
  Pol II occupancy per transcription site:
      Pol2_per_TS = TS_intensity / (W × median single_mRNA_intensity)

  Transcription rate of one active site (mRNA per hour):
      rate_per_TS = Pol2_per_TS × elongation_speed / L
      (elongation_speed in bp/h, e.g. ~ 1.5–3 kb/min from literature; L from this app)

  Transcription rate per cell:
      transcription_rate = gene_copies × fraction_active_sites × rate_per_TS

  Degradation rate (per hour), at steady state:
      degradation_rate = transcription_rate / mRNA_per_cell

  Half-life (hours):
      half_life = ln(2) / degradation_rate
  ```

- A one-sentence copy-pasteable methods statement, e.g. "Probe weight factor W = 0.53 was computed for 48 probes over a 12,345 bp locus using TransQuant-probe-weight (Bahar Halpern & Itzkovitz, 2016)."

## 7. Implementation

### 7.1 Code layout (single code path for app and CLI)

```
src/transquant_w/
    core.py   # pure functions: parse_target, parse_probes, revcomp, locate_probes,
              # probe_profile (N), weight_factor (W), make_plot. No Streamlit imports.
    cli.py    # argparse wrapper around core.py
app.py        # Streamlit UI, imports core.py
web/index.html   # stlite loader (see 7.2)
tests/test_core.py
```

- Dependencies: `numpy`, `matplotlib`. No Biopython (FASTA parsing is a few lines). Both are available as Pyodide packages.
- CLI: `transquant-w --target gene.fa --probes probes.txt [--out prefix] [--revcomp-target] [--no-concat] [--profile-tsv]`. Prints W, L, counts and warnings to stdout; writes plot PNG. Non-zero exit on error.

### 7.2 stlite hosting (designed from scratch, not copied from other apps)

- `web/index.html` loads a **pinned** version of `@stlite/browser` from jsDelivr and calls `mount()` with `entrypoint: "app.py"`, `files` (app.py + core.py), `requirements: ["numpy", "matplotlib"]`, and `streamlitConfig` for the theme.
- Verified during implementation (stlite 1.8.1, Pyodide 0.29.3, Python 3.13):
  - `files` entries may be `{ url: "./app.py" }` and are fetched relative to `index.html`, so the Python sources are served as plain files next to the page. No single-file exporter is needed.
  - `streamlitConfig` accepts only string/number/boolean/null values. List-valued options such as `theme.headingFontSizes` are rejected at startup.
  - A font URL inside `theme.font` / `theme.codeFont` crashes the frontend (react-helmet error). Fonts are loaded with `<link>` (Fira Code) and `@font-face` (Styrene B) in `index.html`; the theme only names the families.
- GitHub Actions workflow on push to `main`: assemble `dist/` (index.html, app.py, core.py, fonts), `actions/upload-pages-artifact`, `actions/deploy-pages`. Pages source set to "GitHub Actions".
- Theme: colours, radius, and typography from the other Streamlit app's `config.toml` are reused as values, no sidebar. **Styrene B trial OTF fonts are kept** (internal use), shipped in `static/fonts/` and declared with `@font-face` in `index.html` (for local `streamlit run`, `.streamlit/config.toml` uses `theme.fontFaces` with `enableStaticServing`).
- Upload size cap via `server.maxUploadSize`.
- Pin all versions (stlite, numpy, matplotlib) for reproducibility.

## 8. Validation (no MATLAB required)

- **Analytic toy cases** in `tests/test_core.py`:
  - 10 probes evenly spread across a synthetic target → W ≈ 0.5.
  - All probes in the first 5% of the target → W ≈ 1.
  - All probes in the last 5% → W ≈ 0.
  - One probe present twice in the target → N(L) = 2 and both sites counted.
  - Unmatched probe → reported and excluded; all unmatched → error.
  - Multi-FASTA concatenation, U→T, lowercase, Excel-style multi-column paste, header row removal.
- **Cross-check against R**: run the existing R script on one real gene + probe set and compare W to the Python result (tolerance 1e-3 to absorb the ±1 base threshold difference). Save the inputs and expected W under `tests/data/` as the golden case.

## 9. Documentation (README order)

1. What the tool does, in two paragraphs, with the paper citation and the statement W ≡ eta (Eq. 7).
2. **Browser app manual** for non-coders: link to the Pages URL, how to get the genomic locus from UCSC (reverse complement for minus-strand genes, exon/intron download), how to paste probes from a spreadsheet, how to read W and the plot, the "how to use W" equations.
3. **Interpretation and pitfalls**: W near 1 / 0.5 / 0 meaning; spliced vs genomic input; wrong strand; unmatched probes.
4. **Advanced**: run locally (`streamlit run app.py`), CLI usage and flags, how the stlite build/deploy works.
5. Citation, licence (MIT, this repo) and note that the algorithm is reimplemented from the paper's supplementary MATLAB code.

Also add `CITATION.cff`.

## 10. Housekeeping (before coding)

- `.gitignore` is the R/pkgdown template and **ignores `docs/`**, so this brief is currently untracked. Replace with a Python-oriented `.gitignore` (`__pycache__/`, `.venv/`, `dist/`, `.DS_Store`) and un-ignore `docs/`.
- Rename `archive/halpern_itzkovitz_2026` → `archive/halpern_itzkovitz_2016` (paper year).
- `archive/`: commit the MATLAB `.m` files and the R script for provenance. **Do not commit the paper markdown** (`archive/halpern_itzkovitz_2016/*.md`); add it to `.gitignore`.

## 11. Non-goals

- No spot detection, image handling, or full TransQuant workflow.
- No mismatch-tolerant probe matching.
- No strand inference from FASTA headers.
- No server component; everything runs in the browser or locally.
