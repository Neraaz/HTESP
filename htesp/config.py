#!/usr/bin/env python
"""Configuration loading for HTESP.

This module replaces the original ``check_json.config()``, which

* looked only in ``./config.json`` and then blindly at ``../../config.json``,
* returned ``None`` when nothing was found (every caller then raised
  ``TypeError`` several frames later),
* re-parsed the file on each of its ~40 call sites, and
* carried the Materials Project API key in the file itself.

The replacement

* searches the current directory and its parents (up to
  :data:`SEARCH_DEPTH` levels) plus ``$HTESP_CONFIG``,
* deep-merges whatever it finds over the packaged default
  (``htesp/data/config.json``), so a configuration written for an older
  schema still has every key the current code reads,
* caches the result per resolved path, and
* takes the API key from ``$MP_API_KEY`` or
  ``~/.config/htesp/credentials`` in preference to the file.

The public entry point keeps the old name and the old return type::

    from htesp.config import config
    cfg = config()          # plain dict, never None
"""
from __future__ import annotations

import copy
import functools
import json
import os
from pathlib import Path
from typing import Any, Mapping

#: file name looked for in the working directory and its parents
CONFIG_FILENAME = "config.json"

#: how many parent directories are searched above the working directory.
#: The bash layer runs helpers from ``R<id>-<comp>/<stage>/`` (two levels),
#: from ``R<id>-<comp>/pressure/R.../relax`` (four levels) and one level
#: deeper again for phonopy supercells, hence five.
SEARCH_DEPTH = 5

#: packaged fallback -- the canonical, fully-populated schema
DEFAULT_CONFIG_PATH = Path(__file__).resolve().parent / "data" / CONFIG_FILENAME

#: environment variable holding the Materials Project API key
API_KEY_ENV = "MP_API_KEY"

#: file consulted when :data:`API_KEY_ENV` is unset
CREDENTIALS_PATH = Path.home() / ".config" / "htesp" / "credentials"

#: placeholder shipped in every example configuration
API_KEY_PLACEHOLDER = "use_your_API_KEY"


class ConfigError(RuntimeError):
    """Raised by :func:`config` in strict mode and by :func:`require_api_key`."""


def _deep_merge(base: Mapping[str, Any], override: Mapping[str, Any]) -> dict:
    """Recursively merge ``override`` onto a copy of ``base``."""
    merged = dict(base)
    for key, value in override.items():
        if key in merged and isinstance(merged[key], Mapping) and isinstance(value, Mapping):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def find_config(start: os.PathLike | str | None = None) -> Path | None:
    """Return the first ``config.json`` at or above ``start``.

    ``$HTESP_CONFIG`` wins over the directory search, so a whole campaign can
    share one file without copying it into every working directory.
    """
    override = os.environ.get("HTESP_CONFIG")
    if override:
        candidate = Path(override).expanduser()
        if candidate.is_dir():
            candidate = candidate / CONFIG_FILENAME
        return candidate if candidate.is_file() else None

    here = Path(start or Path.cwd()).resolve()
    for directory in [here, *here.parents][: SEARCH_DEPTH + 1]:
        candidate = directory / CONFIG_FILENAME
        if candidate.is_file():
            return candidate
    return None


@functools.lru_cache(maxsize=None)
def _load_defaults() -> dict:
    with open(DEFAULT_CONFIG_PATH, "r") as handle:
        return json.load(handle)


@functools.lru_cache(maxsize=None)
def _load_merged(resolved: str | None) -> dict:
    defaults = _load_defaults()
    if resolved is None:
        return copy.deepcopy(defaults)
    with open(resolved, "r") as handle:
        user = json.load(handle)
    return _deep_merge(defaults, user)


def config(start: os.PathLike | str | None = None, strict: bool = False) -> dict:
    """Return the effective configuration as a plain ``dict``.

    Always returns a usable mapping: when no ``config.json`` is found the
    packaged default is returned, so importing a module no longer depends on
    the working directory.  ``strict=True`` raises :class:`ConfigError`
    instead, for commands that genuinely need the user's own settings.
    """
    found = find_config(start)
    if found is None and strict:
        raise ConfigError(
            f"no {CONFIG_FILENAME} found in {Path(start or Path.cwd()).resolve()} "
            f"or its {SEARCH_DEPTH} parent directories, and $HTESP_CONFIG is unset"
        )
    return copy.deepcopy(_load_merged(str(found) if found else None))


def config_path(start: os.PathLike | str | None = None) -> Path | None:
    """The file :func:`config` would read, or ``None`` for the packaged default."""
    return find_config(start)


def clear_cache() -> None:
    """Forget every parsed configuration (used by the tests)."""
    _load_merged.cache_clear()
    _load_defaults.cache_clear()


def api_key(cfg: Mapping[str, Any] | None = None) -> str | None:
    """Materials Project API key: environment, then credentials file, then config.

    Returns ``None`` when only the shipped placeholder is available, so
    callers can give a useful message rather than sending
    ``use_your_API_KEY`` to the MP servers.
    """
    from_env = os.environ.get(API_KEY_ENV, "").strip()
    if from_env:
        return from_env

    if CREDENTIALS_PATH.is_file():
        for line in CREDENTIALS_PATH.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            key, _, value = line.partition("=")
            if key.strip() in (API_KEY_ENV, "api_key", "key"):
                value = value.strip().strip("'\"")
                if value:
                    return value

    cfg = config() if cfg is None else cfg
    try:
        from_file = str(cfg["mpi_key"]["API_KEY"]["key"]).strip()
    except (KeyError, TypeError):
        return None
    if not from_file or from_file == API_KEY_PLACEHOLDER:
        return None
    return from_file


def require_api_key(cfg: Mapping[str, Any] | None = None) -> str:
    """:func:`api_key`, but raise a message that says how to set it."""
    key = api_key(cfg)
    if key:
        return key
    raise ConfigError(
        "No Materials Project API key.  Set it with\n"
        f"    export {API_KEY_ENV}=<your key>\n"
        f"or put '{API_KEY_ENV}=<your key>' in {CREDENTIALS_PATH}.\n"
        "Get a key at https://next-gen.materialsproject.org/api ."
    )


#: keys every configuration must resolve to something usable
REQUIRED_SECTIONS = (
    "job_script", "mpi_key", "download", "conv_test", "magmom", "pseudo",
    "substitute", "pwscf_in", "strain", "wanniertools_input", "kptden",
    "chull_cutoff", "kpt_opt", "elph_mode", "plot",
)


def validate(cfg: Mapping[str, Any] | None = None) -> list[str]:
    """Return a list of human-readable problems with ``cfg`` (empty == fine)."""
    cfg = config() if cfg is None else cfg
    problems: list[str] = []

    for section in REQUIRED_SECTIONS:
        if section not in cfg:
            problems.append(f"missing top-level section {section!r}")

    if api_key(cfg) is None:
        problems.append(
            f"no Materials Project API key (set ${API_KEY_ENV}); "
            "database search and download will not work"
        )

    calc = str(cfg.get("download", {}).get("inp", {}).get("calc", "")).lower()
    if calc not in ("qe", "vasp"):
        problems.append(
            f"download.inp.calc is {calc!r}; expected 'QE' or 'VASP'"
        )

    mode = cfg.get("download", {}).get("mode", "")
    if mode not in ("element", "chemsys", "fromcif", "fromvasp", ""):
        problems.append(
            f"download.mode is {mode!r}; expected element, chemsys, fromcif or fromvasp"
        )

    pseudo = cfg.get("pseudo", {}).get("PSEUDO", {})
    if not isinstance(pseudo, Mapping) or not pseudo:
        problems.append("pseudo.PSEUDO is empty; plane-wave cutoffs cannot be chosen")

    ordering = cfg.get("download", {}).get("element", {}).get("ordering", "NM")
    if ordering is not None and not isinstance(ordering, (str, list, tuple)):
        problems.append(
            f"download.element.ordering must be null, a string or a list, "
            f"got {ordering!r}"
        )
    elif isinstance(ordering, (list, tuple)) and not all(
            isinstance(item, str) for item in ordering):
        problems.append(
            f"download.element.ordering list must contain strings, got {ordering!r}"
        )

    plot = cfg.get("plot", {})
    for limit in ("xlim", "ylim"):
        value = plot.get(limit)
        if value is not None and (not isinstance(value, (list, tuple)) or len(value) != 2):
            problems.append(f"plot.{limit} must be null or a two-element list, got {value!r}")

    return problems
