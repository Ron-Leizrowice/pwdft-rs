#!/usr/bin/env bash
# One-shot bootstrap for the pwdft-rs developer environment.
#
# - Loads .env (copy from .env.example first).
# - Symlinks qe-7.5 → $QE_PATH so Rust/Python validation scripts can reach pw.x.
# - Runs `uv sync` to materialize the Python venv for the validation package.
#
# Idempotent: re-run any time .env changes.

set -euo pipefail

cd "$(dirname "$0")"

if [[ ! -f .env ]]; then
    echo "ERROR: .env is missing. Copy .env.example to .env and edit it." >&2
    exit 1
fi

set -a
# shellcheck disable=SC1091
source .env
set +a

: "${QE_PATH:?QE_PATH must be set in .env}"

if [[ ! -d "$QE_PATH" ]]; then
    echo "ERROR: QE_PATH=$QE_PATH is not a directory." >&2
    exit 1
fi
if [[ ! -x "$QE_PATH/build/bin/pw.x" ]]; then
    echo "WARN: $QE_PATH/build/bin/pw.x not found or not executable — QE-backed tests will fail until you build QE." >&2
fi

link=qe-7.5
if [[ -L $link ]]; then
    current=$(readlink "$link")
    if [[ $current != "$QE_PATH" ]]; then
        echo "Retargeting $link: $current → $QE_PATH"
        rm "$link"
        ln -s "$QE_PATH" "$link"
    else
        echo "$link → $QE_PATH (already correct)"
    fi
elif [[ -e $link ]]; then
    echo "ERROR: $link exists and is not a symlink; refusing to overwrite." >&2
    exit 1
else
    ln -s "$QE_PATH" "$link"
    echo "Created $link → $QE_PATH"
fi

echo "Running uv sync..."
uv sync

echo "Installing pre-commit hooks..."
uv run prek install

echo
echo "Setup complete."
echo "  QE:     $link → $QE_PATH"
echo "  Python: .venv/ (via uv sync)"
echo "  Hooks:  .git/hooks/pre-commit (via prek)"
echo
echo "Next steps:"
echo "  cargo test                             # Tier-1"
echo "  cargo test -- --ignored                # Tier-2 heavy"
echo "  cargo build --release                  # production build"
