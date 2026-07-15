"""Regression tests for the LOWESS delta speedup knob (config.lowess_delta_frac).

Covers (a) the _lowess_delta helper's edge cases and (b) an equivalence check
that delta-interpolated LOWESS (frac=0.01, the default) reproduces the exact
O(n^2) fit (frac=0.0) to a tight tolerance on the intensity_trend variance
prior, with logFC bit-identical and no significance-call flips.

The import falls back to the vendored `moderated_limma` package so this file
runs both inside proteomics-toolkit (the PR target) and against the standalone
fork in PostdocToolbox. For an upstream PR the fallback branch can be dropped.
Run locally:  python -m pytest tests/ -q   (from the tool directory)
"""
import numpy as np
import pandas as pd
import pytest

try:  # proteomics-toolkit (PR target)
    from proteomics_toolkit.statistical_analysis import (
        StatisticalConfig,
        _lowess_delta,
        run_moderated_linear_model,
    )
except ImportError:  # vendored fork (PostdocToolbox)
    from moderated_limma.statistical_analysis import (
        StatisticalConfig,
        _lowess_delta,
        run_moderated_linear_model,
    )


class _Cfg:
    """Minimal stand-in exposing only lowess_delta_frac."""

    def __init__(self, frac):
        self.lowess_delta_frac = frac


class TestLowessDeltaHelper:
    def test_default_is_one_percent_of_range(self):
        x = np.array([0.0, 10.0, 5.0, 7.0])
        assert _lowess_delta(x, StatisticalConfig()) == pytest.approx(0.1)

    def test_zero_frac_forces_exact(self):
        assert _lowess_delta(np.array([0.0, 10.0]), _Cfg(0.0)) == 0.0

    def test_negative_frac_forces_exact(self):
        assert _lowess_delta(np.array([0.0, 10.0]), _Cfg(-1.0)) == 0.0

    def test_nan_frac_forces_exact(self):
        assert _lowess_delta(np.array([0.0, 10.0]), _Cfg(float("nan"))) == 0.0

    def test_single_unique_value_gives_zero(self):
        assert _lowess_delta(np.array([3.0, 3.0, 3.0]), _Cfg(0.01)) == 0.0

    def test_empty_after_nonfinite_filter_gives_zero(self):
        assert _lowess_delta(np.array([np.nan, np.inf, -np.inf]), _Cfg(0.01)) == 0.0

    def test_ignores_nonfinite_in_range(self):
        x = np.array([0.0, np.nan, 10.0, np.inf])
        assert _lowess_delta(x, _Cfg(0.01)) == pytest.approx(0.1)

    def test_none_config_uses_default(self):
        # The deqms path calls _fit_count_dependent_prior(..., config=None),
        # which forwards config=None here; getattr default must be 0.01.
        assert _lowess_delta(np.array([0.0, 100.0]), None) == pytest.approx(1.0)


def _make_intensity_fixture(total_rows=300, with_effect_rows=30, seed=42):
    """Unpaired intensity_trend fixture with a wide dynamic range.

    The wide range (raw means spanning ~2^6..2^16) makes log-mean-intensity
    spread large enough that delta=0.01*range genuinely skips points, so the
    equivalence test is not trivially exact.
    """
    rng = np.random.default_rng(seed)
    samples_a = [f"A{i}" for i in range(6)]
    samples_b = [f"B{i}" for i in range(6)]
    all_samples = samples_a + samples_b

    base = rng.uniform(6.0, 16.0, size=(total_rows, 1))
    values = base + rng.normal(0.0, 0.5, size=(total_rows, 12))
    values[:with_effect_rows, 6:] += 1.5  # planted effect in the treatment group

    features = [f"P{i:04d}" for i in range(total_rows)]
    feature_data = pd.DataFrame(values, index=features, columns=all_samples)  # log2 scale
    metadata_df = pd.DataFrame({"Sample": all_samples, "Group": ["Control"] * 6 + ["Treatment"] * 6})

    config = StatisticalConfig()
    config.analysis_type = "unpaired"
    config.group_column = "Group"
    config.group_labels = ["Control", "Treatment"]
    config.log_transform_before_stats = False
    config.statistical_test_method = "moderated_linear_model"
    config.moderation = "intensity_trend"
    config._raw_feature_data = 2 ** feature_data  # trend prior expects raw intensities
    return feature_data, metadata_df, config


class TestLowessDeltaEquivalence:
    def test_delta_actually_engages(self):
        """Guard: the fixture's log-mean spread yields a positive delta, so the
        equivalence below really exercises the interpolation path."""
        _, _, config = _make_intensity_fixture()
        log_mean = np.log(config._raw_feature_data.mean(axis=1).to_numpy())
        assert _lowess_delta(log_mean, config) > 0.0

    def test_delta_matches_exact_intensity_trend(self):
        feat, meta, cfg_exact = _make_intensity_fixture()
        cfg_exact.lowess_delta_frac = 0.0
        exact = run_moderated_linear_model(feat, meta, cfg_exact).set_index("Protein")

        feat2, meta2, cfg_fast = _make_intensity_fixture()
        cfg_fast.lowess_delta_frac = 0.01
        fast = run_moderated_linear_model(feat2, meta2, cfg_fast).set_index("Protein")
        fast = fast.reindex(exact.index)

        # logFC is the OLS contrast; it never depends on the variance prior.
        assert np.array_equal(exact["logFC"].to_numpy(), fast["logFC"].to_numpy())
        # P-values agree to a tight tolerance on the smooth variance trend.
        dp = (exact["P.Value"] - fast["P.Value"]).abs()
        assert dp.max() < 1e-3
        # No significance-call flips at the usual 0.05 threshold.
        flips = ((exact["P.Value"] < 0.05) != (fast["P.Value"] < 0.05)).sum()
        assert flips == 0
