# Updating the public repository

Repository: [Homabazireh/fine-mapping](https://github.com/Homabazireh/fine-mapping).
The working branch is `main`. Publish Python and R code, documentation and
reviewed aggregate tables only. The user requested that all figures and reports
containing embedded figures be removed from this public repository.

All saved figures and full results belong in private `Homabazireh/fine-mapping-2`.
Its private release contains the genomic project data; its `all_figures/` folder
contains browsable copies of the saved figures. Never upload participant-level
data, genotype files, aliases, credentials or private archives to the public repo.

Before updating, check `git remote -v` and `git status`. Run the Python and R
checks documented in the project README and `code/R/README.md`, then stage
specific reviewed files. Inspect `git diff --cached --stat` and
`git diff --cached --check`, commit, and push with `git push origin main`.
Do not force-push. Copy new figure exports only to the verified private destination.

Removing figure files from the current branch does not rewrite older Git commits.
