"""Vendored moderated-linear-model ("limma"-style) statistics.

Forked verbatim from proteomics-toolkit v26.4.0 (commit 25a0133, 2026-07-02):
  statistical_analysis.py, normalization.py, preprocessing.py
so this tool runs standalone, decoupled from that repo's ongoing changes.

Public API used by intensity_prior_limma.py:
  - run_moderated_linear_model(feature_data, metadata_df, config)
  - StatisticalConfig
"""
from .statistical_analysis import StatisticalConfig, run_moderated_linear_model

__all__ = ["run_moderated_linear_model", "StatisticalConfig"]
__vendored_from__ = "proteomics-toolkit v26.4.0 @ 25a0133"
