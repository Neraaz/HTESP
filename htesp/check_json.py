#!/usr/bin/env python
"""Backwards-compatible wrapper around :mod:`htesp.config`.

Historically every module did ``from check_json import config``.  That name is
kept so the science modules read the same as before, but the implementation is
now the searching, caching, default-merging loader in :mod:`htesp.config`.
"""
from htesp.config import (  # noqa: F401  (re-exported for callers)
    CONFIG_FILENAME,
    ConfigError,
    api_key,
    clear_cache,
    config,
    config_path,
    find_config,
    require_api_key,
    validate,
)

__all__ = [
    "CONFIG_FILENAME", "ConfigError", "api_key", "clear_cache", "config",
    "config_path", "find_config", "require_api_key", "validate",
]
