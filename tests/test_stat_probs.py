"""Tests for the compiled Cython module get_stat_probs."""
import math
import numpy as np
import pytest
from scipy.stats import kendalltau as scipy_kt

from get_stat_probs import (
    kt,
    generate_base_reference,
    get_waveform_list,
    make_references,
    farctanh,
    periodic,
)


class TestKendallTau:
    def test_perfect_agreement(self):
        tau, _ = kt([1, 2, 3, 4, 5], [1, 2, 3, 4, 5])
        assert tau == pytest.approx(1.0)

    def test_perfect_disagreement(self):
        tau, _ = kt([1, 2, 3, 4, 5], [5, 4, 3, 2, 1])
        assert tau == pytest.approx(-1.0)

    def test_matches_scipy(self):
        x = [1, 3, 2, 5, 4, 6, 8, 7]
        y = [2, 1, 4, 3, 5, 7, 6, 8]
        tau, _ = kt(x, y)
        scipy_tau, _ = scipy_kt(x, y)
        assert tau == pytest.approx(scipy_tau, abs=1e-5)

    def test_returns_two_values(self):
        result = kt([1, 2, 3], [1, 2, 3])
        assert len(result) == 2

    def test_p_value_in_range(self):
        _, p = kt([1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6])
        assert 0.0 <= p <= 1.0

    def test_uncorrelated_tau_near_zero(self):
        # A shuffled sequence should give tau close to 0 on average
        rng = np.random.default_rng(42)
        taus = [kt(list(range(20)), list(rng.permutation(20)))[0] for _ in range(50)]
        assert abs(np.mean(taus)) < 0.3


class TestGenerateBaseReference:
    HEADER = list(range(0, 24, 2))  # 12 evenly-spaced ZT timepoints

    def test_cosine_correct_length(self):
        ref = generate_base_reference(self.HEADER, 'cosine', 24.0, 0.0, 12.0)
        assert len(ref) == len(self.HEADER)

    def test_trough_correct_length(self):
        ref = generate_base_reference(self.HEADER, 'trough', 24.0, 0.0, 12.0)
        assert len(ref) == len(self.HEADER)

    def test_cosine_peak_at_phase_zero(self):
        # With phase=0 the cosine peak lands at ZT0 (first element)
        ref = generate_base_reference(self.HEADER, 'cosine', 24.0, 0.0, 12.0)
        assert ref[0] == pytest.approx(max(ref), abs=1e-6)

    def test_cosine_peak_shifts_with_phase(self):
        ref0 = generate_base_reference(self.HEADER, 'cosine', 24.0, 0.0, 12.0)
        ref6 = generate_base_reference(self.HEADER, 'cosine', 24.0, 6.0, 12.0)
        # Different phases should produce different waveforms
        assert not np.allclose(ref0, ref6)

    def test_cosine_values_bounded(self):
        ref = generate_base_reference(self.HEADER, 'cosine', 24.0, 0.0, 12.0)
        assert all(-1.001 <= v <= 1.001 for v in ref)

    def test_trough_values_bounded(self):
        ref = generate_base_reference(self.HEADER, 'trough', 24.0, 0.0, 12.0)
        assert all(-1.001 <= v <= 1.001 for v in ref)


class TestGetWaveformList:
    # get_waveform_list pre-allocates int(n_phases * n_widths / 2) slots and
    # deduplicates (phase, nadir) pairs. Use the canonical reference-file inputs
    # (12 phases × 11 widths) so the formula holds.
    PHASES = np.array(list(range(0, 24, 2)))   # 12 values
    WIDTHS = np.array(list(range(2, 24, 2)))   # 11 values

    def test_returns_2d_array(self):
        triples = get_waveform_list(np.array([24]), self.PHASES, self.WIDTHS)
        assert triples.ndim == 2
        assert triples.shape[1] == 3

    def test_positive_row_count(self):
        triples = get_waveform_list(np.array([24]), self.PHASES, self.WIDTHS)
        assert triples.shape[0] > 0

    def test_period_value_present(self):
        triples = get_waveform_list(np.array([24]), self.PHASES, self.WIDTHS)
        assert all(t[0] == 24 for t in triples)

    def test_phase_values_in_input_range(self):
        triples = get_waveform_list(np.array([24]), self.PHASES, self.WIDTHS)
        for t in triples:
            assert t[1] in self.PHASES

    def test_row_count_is_nonzero(self):
        # With the standard 12-phase × 11-width inputs the result is non-empty
        triples = get_waveform_list(np.array([24]), self.PHASES, self.WIDTHS)
        assert triples.shape[0] > 0


class TestMakeReferences:
    HEADER = list(range(0, 24, 2))  # 12 ZT timepoints
    PHASES = np.array(list(range(0, 24, 2)))
    WIDTHS = np.array(list(range(2, 24, 2)))

    def test_returns_dict(self):
        triples = get_waveform_list(np.array([24]), self.PHASES, self.WIDTHS)
        dref = make_references(self.HEADER, triples)
        assert isinstance(dref, dict)

    def test_key_count_matches_triples(self):
        triples = get_waveform_list(np.array([24]), self.PHASES, self.WIDTHS)
        dref = make_references(self.HEADER, triples)
        assert len(dref) == triples.shape[0]

    def test_each_reference_correct_length(self):
        triples = get_waveform_list(np.array([24]), self.PHASES, self.WIDTHS)
        dref = make_references(self.HEADER, triples)
        for ref in dref.values():
            assert len(ref) == len(self.HEADER)

    def test_keys_are_tuples_of_three(self):
        triples = get_waveform_list(np.array([24]), self.PHASES, self.WIDTHS)
        dref = make_references(self.HEADER, triples)
        for key in dref.keys():
            assert len(key) == 3

    def test_trough_waveform_differs_from_cosine(self):
        triples = get_waveform_list(np.array([24]), self.PHASES, self.WIDTHS)
        dref_cos = make_references(self.HEADER, triples, waveform='cosine')
        dref_tr = make_references(self.HEADER, triples, waveform='trough')
        key = list(dref_cos.keys())[0]
        assert not np.allclose(dref_cos[key], dref_tr[key])


class TestHelpers:
    def test_farctanh_at_zero(self):
        assert farctanh(0.0) == pytest.approx(0.0)

    def test_farctanh_positive(self):
        result = farctanh(0.5)
        assert result == pytest.approx(math.atanh(0.5), rel=1e-4)

    def test_farctanh_clips_near_one(self):
        # Should not blow up at 1.0 (clips to 0.99)
        result = farctanh(1.0)
        assert math.isfinite(result)
        assert result == pytest.approx(math.atanh(0.99), rel=1e-4)

    def test_periodic_in_range(self):
        for x in [-24, -12, 0, 12, 24, 36]:
            result = periodic(float(x))
            assert -12 < result <= 12
