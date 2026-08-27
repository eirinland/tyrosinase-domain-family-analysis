#!/usr/bin/env python3
"""Resolve the analysis root so scripts run from a fresh clone of this repository.

Several scripts were written against an absolute path to a working directory
named Super_reference_pipeline, which is not part of this deposit. Every input
those scripts read is tracked here, so the only thing that was unreachable was
the path itself.

Resolution order, first hit wins.
  1. $PPO_BASE, if set. Lets you point at a working tree elsewhere.
  2. The repository root, found by walking up from this file. This is the case
     that makes a fresh clone work with no configuration.
  3. The legacy Super_reference_pipeline directory, if it happens to exist
     alongside the clone. Kept so existing local setups behave as before.

Usage in a script that previously hardcoded BASE:

    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from repo_paths import BASE

or, from any depth, use `find_base()` directly.
"""
import os

__all__ = ["BASE", "find_base", "resolve"]

# Files that identify the analysis root. Chosen because they are tracked, sit at
# the top level, and are read by the scripts that need BASE.
_MARKERS = (
    "canonical_criteria_all_ca.csv",
    "taxonomy_lookup.csv",
    "1_filtering/final_pools/three_pool_assignment_final.csv",
)


def _looks_like_base(path):
    if not path or not os.path.isdir(path):
        return False
    return any(os.path.exists(os.path.join(path, m)) for m in _MARKERS)


def find_base():
    """Return the analysis root, or raise with a message naming what to set."""
    env = os.environ.get("PPO_BASE")
    if env:
        if not _looks_like_base(env):
            raise SystemExit(
                f"PPO_BASE is set to {env!r} but no expected data file was found "
                f"there. Expected one of {', '.join(_MARKERS)}."
            )
        return os.path.abspath(env)

    here = os.path.dirname(os.path.abspath(__file__))
    candidate = here
    for _ in range(6):
        if _looks_like_base(candidate):
            return candidate
        parent = os.path.dirname(candidate)
        if parent == candidate:
            break
        candidate = parent

    legacy = os.path.join(os.path.dirname(here), "Super_reference_pipeline")
    if _looks_like_base(legacy):
        return legacy

    raise SystemExit(
        "Could not locate the analysis root. Run from inside a clone of this "
        "repository, or set PPO_BASE to a directory containing "
        f"{_MARKERS[0]}."
    )


def resolve(*parts):
    """Join `parts` onto the analysis root."""
    return os.path.join(find_base(), *parts)


BASE = find_base()

if __name__ == "__main__":
    print(BASE)
    for m in _MARKERS:
        p = os.path.join(BASE, m)
        print(f"  {'found  ' if os.path.exists(p) else 'MISSING'}  {m}")
