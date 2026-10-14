#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if ! repo_root="$(git -C "$script_dir" rev-parse --show-toplevel 2>/dev/null)"; then
  echo "Error: sync-wiki.sh must be inside a Git checkout." >&2
  exit 1
fi
remote_url="$(git -C "$repo_root" remote get-url origin)"
source_branch="$(git -C "$repo_root" branch --show-current)"

if [[ -z "$source_branch" ]]; then
  echo "Error: cannot determine the current branch (detached HEAD)." >&2
  exit 1
fi

case "$remote_url" in
  https://github.com/*)
    repository_slug="${remote_url#https://github.com/}"
    ;;
  http://github.com/*)
    repository_slug="${remote_url#http://github.com/}"
    ;;
  git@github.com:*)
    repository_slug="${remote_url#git@github.com:}"
    ;;
  ssh://git@github.com/*)
    repository_slug="${remote_url#ssh://git@github.com/}"
    ;;
  *)
    echo "Error: origin must be a GitHub HTTPS or SSH URL." >&2
    exit 1
    ;;
esac

repository_slug="${repository_slug%.git}"
if [[ "$repository_slug" != */* || "$repository_slug" == */*/* || "$repository_slug" == *[!A-Za-z0-9._/-]* ]]; then
  echo "Error: could not determine a valid GitHub owner/repository from origin." >&2
  exit 1
fi

wiki_url="${remote_url%.git}.wiki.git"
source_branch_url="${source_branch//\//%2F}"
source_doc="$repo_root/docs/GenSDF.md"

if [[ ! -f "$source_doc" ]]; then
  echo "Error: documentation file not found: $source_doc" >&2
  exit 1
fi

work_dir="$(mktemp -d)"
trap 'rm -rf "$work_dir"' EXIT
wiki_dir="$work_dir/wiki"

if ! git clone "$wiki_url" "$wiki_dir"; then
  echo "Error: could not clone the wiki repository." >&2
  echo "Ensure the wiki is enabled and initialized, and that your Git credentials can push to it." >&2
  exit 1
fi

cp "$source_doc" "$wiki_dir/Home.md"
REPOSITORY_SLUG="$repository_slug" SOURCE_BRANCH_URL="$source_branch_url" \
  perl -0pi -e '
    s{(https://github\.com/\Q$ENV{REPOSITORY_SLUG}\E/(?:blob|tree)/)main(?=/)}
     {$1 . $ENV{SOURCE_BRANCH_URL}}gex
  ' "$wiki_dir/Home.md"

cd "$wiki_dir"
git add Home.md
if git diff --cached --quiet -- Home.md; then
  echo "Wiki Home.md is already up to date."
  exit 0
fi

if ! author_name="$(git -C "$repo_root" config user.name)" || [[ -z "$author_name" ]]; then
  echo "Error: configure git user.name before syncing the wiki." >&2
  exit 1
fi
if ! author_email="$(git -C "$repo_root" config user.email)" || [[ -z "$author_email" ]]; then
  echo "Error: configure git user.email before syncing the wiki." >&2
  exit 1
fi

git config user.name "$author_name"
git config user.email "$author_email"
git commit -m "Sync GenSDF documentation"
git push
echo "Wiki Home.md updated from branch ${source_branch}."
