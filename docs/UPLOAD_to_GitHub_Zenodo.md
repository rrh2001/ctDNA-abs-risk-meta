# Uploading to GitHub and minting a Zenodo DOI

This repository is a self-contained snapshot. No build step is required. The
steps below publish the code and data on GitHub and obtain a permanent,
citable Zenodo DOI for the manuscript.

## 1. What must NOT be uploaded

The archive already excludes third-party patient-level data and local/working
files via `.gitignore`. Before uploading, double-check that **none** of the
following are present anywhere in the repository (the supplied audit reports
`CLEAN`):

* cleaned appendix IPD (`ipd_clean_*.csv`), the raw `85_Li2022.tsv` /
  `24_Ogaard2022.tsv`, or any `real_IPD/` folder;
* the original Chinese-language extraction workbooks (e.g. `全表.xlsx`);
* personal phone numbers, email addresses beyond the author contact addresses,
  or internal file paths;
* caches and transient outputs (`__pycache__/`, `*.pyc`, `.Rhistory`,
  `.Rproj.user/`, `M*_out/`, `*_chunk*`, large `M4_all_iterations` files).

These data are from published appendices of other groups and are used only for
local validation; they must not be redistributed (see `LICENSE` and the Data
availability statement in the manuscript).

## 2. Create the GitHub repository

1. Sign in to GitHub and click **New repository**.
2. Name it `ctdna-abs-risk-meta` (the README uses the URL
   `https://github.com/renruhua/ctdna-abs-risk-meta`; keep this owner/name or
   update every occurrence of that URL).
3. Set it to **Public** (required for Zenodo).
4. Do **not** add a README, `.gitignore`, or license through the web form — the
   snapshot already contains them.

## 3. Upload the files (choose one)

### Option A — GitHub Desktop (recommended, no command line)

1. Install GitHub Desktop and sign in.
2. **File → Add local repository**, choose the unzipped `ctdna-abs-risk-meta`
   folder. If it is not yet a git repository, choose **create a repository here**
   and accept the existing `.gitignore`.
3. Write a summary such as `Initial release: reconstructed-IPD meta-analysis`,
   commit to `main`, then **Publish repository** / **Push** to the new GitHub
   repository.

### Option B — Git command line

```bash
cd ctdna-abs-risk-meta
git init -b main
git add .
git commit -m "Initial release: reconstructed-IPD meta-analysis with open risk calculator"
git remote add origin https://github.com/renruhua/ctdna-abs-risk-meta.git
git push -u origin main
```

### Option C — Web upload (no software)

On the new repository page choose **Add file → Upload files**, drag in the
contents of the unzipped folder **including the hidden `.gitignore` and the
subfolder structure** (`R/`, `python/`, `data/`, `outputs/`, `figures/`,
`app/`, `docs/`, `supplementary/`, `data-raw/`), then commit. GitHub Desktop
or git is preferable because web uploads can flatten empty folders.

## 4. Verify the public repository

* The README renders on the repository home page with the directory tree and
  badge placeholders.
* `figures/` shows the 12 high-resolution image files, `outputs/` shows 35 CSVs,
  and `app/app.R` is present.
* Reviewers can launch the calculator locally with
  `shiny::runApp("app")`; the hosted copy remains at
  `https://renruhua.shinyapps.io/ctDNA_risk_calculator/`.
* Reproducibility: `run_all.R` followed by the five Python steps reproduces the
  released tables and figures (see `docs/runbook.md`).

## 5. Mint the Zenodo DOI

1. Sign in to <https://zenodo.org> using the **GitHub** login.
2. Under **Account → Linked accounts (GitHub)**, connect GitHub and grant
   access to `ctdna-abs-risk-meta`.
3. In Zenodo's GitHub list, toggle the repository **on**.
4. On GitHub, create a release: **Releases → Draft a new release**, tag
   `v1.0.0`, title `v1.0.0`, publish the release.
5. Zenodo archives the release automatically and issues a DOI of the form
   `10.5281/zenodo.XXXXXXXXX`.
6. Replace the placeholder `10.5281/zenodo.ZENODO_DOI` in **`README.md`** and
   **`CITATION.cff`** with the real DOI (search the files for `ZENODO_DOI`),
   commit, and either re-tag `v1.0.1` or update the release. The DOI always
   resolves to the archived version; updating the files is a documentation
   nicety, not a requirement for the DOI to work.

## 6. After publication

* Keep the repository and the Zenodo deposit consistent with the manuscript
  version submitted (the same tables and figures).
* If the manuscript is revised, commit the changes and create a new release tag
  (`v1.1.0`); Zenodo issues a new version DOI while preserving the old one.
* Report the GitHub URL and Zenodo DOI in the manuscript data/code availability
  statement and the cover letter.
