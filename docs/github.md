# Updating the GitHub repository

Repository: [Homabazireh/fine-mapping](https://github.com/Homabazireh/fine-mapping).
The working branch is `main`. Open [panel c and its report](figures/README.md) to
view the shared figure, with PNG, PDF and SVG download links. On a GitHub PDF
page, use the download button to save the file.

## Push future changes from the analysis workspace

```bash
cd /mnt/storage/ahbazih1/Fine_mapping/new_analysis
git status
.venv/bin/python code/run_tests.py
git diff --check
git add README.md .gitignore code/ config/ docs/ literature_review/ tests/
git diff --cached --stat
git commit -m "Describe the changes"
git push origin main
```

`git add` selects files for a commit. `git commit` saves a local version with
a description. `git push` sends that commit to GitHub. Review the selected
files before committing. If you also change tracked documentation inside
`data/` or `results/`, stage those specific README files by name.

Raw genomic files, participant tables, full analysis outputs and publisher
PDF reading copies are excluded by `.gitignore`. The selected aggregate figure
copies and the focused report under `docs/figures/` are tracked. Regenerating the local report creates
a new working figure directory; it does not update the published copies
automatically. Copy the reviewed figure and focused report exports into `docs/figures/`, update
its captions and provenance manifest, then commit and push them.

If Git reports that the remote has newer commits, use `git fetch origin` and
review the differences before integrating them. Do not force-push over another
update.
