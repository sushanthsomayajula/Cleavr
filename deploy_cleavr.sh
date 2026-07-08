#!/usr/bin/env bash
# Deploys Cleavr to GitHub Pages. Run this from Terminal on your own Mac:
#   cd /Volumes/sushanthusb/tnbc-project
#   bash deploy_cleavr.sh
#
# What it does: installs the GitHub CLI if missing, logs you into GitHub in
# your browser (no password ever typed here), creates a public repo named
# "cleavr", pushes everything except the large raw TCGA data, and prints
# the exact steps to flip on GitHub Pages.

set -e
cd "$(dirname "$0")"

echo "== Checking for GitHub CLI =="
if ! command -v gh >/dev/null 2>&1; then
  if command -v brew >/dev/null 2>&1; then
    echo "Installing gh via Homebrew..."
    brew install gh
  else
    echo "GitHub CLI (gh) isn't installed and Homebrew isn't available."
    echo "Install Homebrew first: https://brew.sh"
    echo "Then re-run: brew install gh && bash deploy_cleavr.sh"
    exit 1
  fi
fi

echo "== Checking GitHub auth =="
if ! gh auth status >/dev/null 2>&1; then
  echo "Opening your browser to log into GitHub..."
  gh auth login --web --git-protocol https
fi

echo "== Preparing local repo =="
if [ ! -d .git ]; then
  git init
  git branch -m main
fi
git add -A
git commit -m "Cleavr: GNRHR expression across TNBC subtypes" --allow-empty

echo "== Creating/pushing GitHub repo =="
if gh repo view cleavr >/dev/null 2>&1; then
  echo "Repo 'cleavr' already exists on your account — pushing to it."
  git remote remove origin 2>/dev/null || true
  git remote add origin "https://github.com/$(gh api user -q .login)/cleavr.git"
  git push -u origin main
else
  gh repo create cleavr --public --source=. --remote=origin --push
fi

GH_USER=$(gh api user -q .login)
echo ""
echo "== Done pushing. Last step (30 seconds, one time only): =="
echo "1. Open: https://github.com/${GH_USER}/cleavr/settings/pages"
echo "2. Under 'Build and deployment' > Source, choose 'Deploy from a branch'"
echo "3. Branch: main, folder: / (root) -> Save"
echo ""
echo "Your site will be live in ~1 minute at:"
echo "  https://${GH_USER}.github.io/cleavr/"
