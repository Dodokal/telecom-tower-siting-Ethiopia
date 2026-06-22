#!/usr/bin/env bash
# ============================================================================
# setup.sh - one-shot Git Bash setup for the GitHub repository
# ----------------------------------------------------------------------------
# Usage from Git Bash:
#   cd "/k/telecom-tower-siting-ethiopia"
#   bash setup.sh
#
# What it does:
#   1. Verifies you are inside the project directory
#   2. Configures Git globally (user.name, user.email) if not set
#   3. Initialises the local Git repo
#   4. Stages all files (respecting .gitignore)
#   5. Makes the initial commit
#   6. Adds the GitHub remote
#   7. Pushes to GitHub
#
# Before running, edit the variables in the CONFIG block below.
# ============================================================================

set -e   # exit on first error

# ----------------------------- CONFIG -------------------------------------
GITHUB_USERNAME="Dodokal"
REPO_NAME="telecom-tower-siting-ethiopia"
YOUR_NAME="Kalid Hassen Yasin"
YOUR_EMAIL="kalid.yasin@plus.ac.at"
INITIAL_COMMIT_MSG="Initial commit: pipeline skeleton, configs, sensitivity results"
# ---------------------------------------------------------------------------

REMOTE_URL="https://github.com/${GITHUB_USERNAME}/${REPO_NAME}.git"

echo "=========================================="
echo "  GitHub setup for ${REPO_NAME}"
echo "=========================================="

# 1. Verify we're in the right directory
if [[ ! -f "README.md" ]] || [[ ! -d "R" ]]; then
  echo "ERROR: This does not look like the project directory."
  echo "  Expected to find README.md and an R/ folder here."
  echo "  cd into the project root and run this script again."
  exit 1
fi

# 2. Configure Git globally (only if not already set)
if [[ -z "$(git config --global user.name)" ]]; then
  echo "Setting global Git user.name -> ${YOUR_NAME}"
  git config --global user.name "${YOUR_NAME}"
fi
if [[ -z "$(git config --global user.email)" ]]; then
  echo "Setting global Git user.email -> ${YOUR_EMAIL}"
  git config --global user.email "${YOUR_EMAIL}"
fi
git config --global init.defaultBranch main 2>/dev/null || true
git config --global core.autocrlf true 2>/dev/null || true

# 3. Initialise local Git repo (skip if already done)
if [[ ! -d ".git" ]]; then
  echo "Initialising Git repository..."
  git init
  git branch -M main
else
  echo "Git repository already initialised — skipping init."
fi

# 4. Stage everything
echo ""
echo "Staging files..."
git add .

# 5. Show what will be committed
echo ""
echo "Files that will be committed:"
git diff --cached --name-only | head -40
echo ""
echo "...and any others below this (truncated to first 40)."
echo ""

# 6. Commit
if git diff --cached --quiet; then
  echo "Nothing to commit — working tree clean."
else
  echo "Creating initial commit..."
  git commit -m "${INITIAL_COMMIT_MSG}"
fi

# 7. Add remote (if not present)
if ! git remote get-url origin > /dev/null 2>&1; then
  echo ""
  echo "Adding remote origin -> ${REMOTE_URL}"
  git remote add origin "${REMOTE_URL}"
else
  echo ""
  echo "Remote origin already configured:"
  git remote get-url origin
fi

# 8. Push to GitHub
echo ""
echo "=========================================="
echo "  Pushing to GitHub..."
echo "=========================================="
echo ""
echo "If this is your first push, Git will prompt for credentials:"
echo "  Username: ${GITHUB_USERNAME}"
echo "  Password: paste a personal access token (NOT your GitHub password)"
echo "  Generate one at https://github.com/settings/tokens/new (scope: repo)"
echo ""

git push -u origin main

echo ""
echo "=========================================="
echo "  DONE!"
echo "=========================================="
echo ""
echo "Open in your browser:"
echo "  https://github.com/${GITHUB_USERNAME}/${REPO_NAME}"
echo ""
echo "Next steps:"
echo "  1. Verify the repo is visible on GitHub (link above)"
echo "  2. Add the manuscript PDF to docs/ when ready"
echo "  3. On paper acceptance: flip repo to Public, add DOI to README"
