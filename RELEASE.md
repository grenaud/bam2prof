# Release checklist: `-is` / `-comp` feature addition

Commands to run from the repo root (`/home/garen18/projects/bam2prof`) to commit this
session's changes and cut a release. Nothing here has been run for you — copy/paste as you go.

## 0. What's actually part of this change

```
git status --short
```

should show:

```
 M README.md
 M src/bam2prof.cpp
 M src/utilityFunctions.cpp
 M src/utilityFunctions.h
?? img/
?? src/plot_bam2prof_comp.py
```

Everything else `git status` lists is **not** part of this work and shouldn't be swept in:

- `src/bam2prof`, `src/diffprof`, `src/*.o` — build artifacts. There's no `.gitignore` in
  this repo, so they show as untracked; don't add them. (Optional: add one — see step 4.)
- ` M lib/htslib`, ` M lib/samtools`, `? lib/libgab` — pre-existing dirty submodule state
  from before this session, unrelated to `-is`/`-comp`. Leave these alone unless you mean
  to deal with them separately.

## 1. Stage only the real changes

```bash
git add README.md src/bam2prof.cpp src/utilityFunctions.cpp src/utilityFunctions.h \
        src/plot_bam2prof_comp.py img/base_composition_example.png
git status --short   # double check the list matches the one above
```

## 2. Commit

```bash
git commit -m "$(cat <<'EOF'
Add -is fragment size and -comp base composition profiling

-is integrates the standalone insertsize/insize tool as an option: reports
count<TAB>length over the same filtered read stream bam2prof already builds,
using a map rather than a fixed-size histogram.

-comp reproduces mapDamage's base composition plot around the fragment ends.
Without -fa it reports composition inside the fragment only (from each read's
MD tag); with -fa it also reports -around bp of true flanking reference
sequence on either side, strand-corrected. Includes src/plot_bam2prof_comp.py
to render it.
EOF
)"
```

## 3. Tag and push

No tags exist yet in this repo (`git tag -l` is empty), so pick a version. This is a
backwards-compatible feature addition (new opt-in flags, no change to existing output) —
that's a **minor** bump under semver. If you've been treating current `master` as an
unreleased `v1.0.0`, tag this `v1.1.0`; if this is genuinely the first tagged release,
use `v1.0.0` instead. Swap `v1.1.0` below for whichever you pick:

```bash
git push origin master

git tag -a v1.1.0 -m "Add -is fragment size and -comp base composition profiling"
git push origin v1.1.0
```

## 4. (Optional) stop build artifacts from showing as untracked

```bash
cat > .gitignore <<'EOF'
src/*.o
src/bam2prof
src/diffprof
src/__pycache__/
EOF
git add .gitignore
git commit -m "Add .gitignore for build artifacts"
git push origin master
```

## 5. Create the GitHub release

Requires `gh` (already authenticated as `grenaud` in this environment):

```bash
gh release create v1.1.0 \
  --repo grenaud/bam2prof \
  --title "v1.1.0" \
  --notes "$(cat <<'EOF'
## Added
- \`-is <file>\`: fragment/insert size distribution (count, length), integrating the
  standalone insertsize tool as an option.
- \`-comp\` / \`-around <N>\`: mapDamage-style base composition profile around the
  fragment ends, with or without a reference genome (\`-fa\`).
- \`src/plot_bam2prof_comp.py\`: plots the composition profile (PDF + PNG).
EOF
)"
```

Drop `--repo grenaud/bam2prof` if you're running this from inside the repo (`gh` infers
it from the current directory's remote).
