#!/usr/bin/env python
"""
Unit test script for the functions in utils/histogram.py
"""

import unittest

import numpy as np

from puma.utils import logger, set_log_level
from puma.utils.histogram import hist_ratio, hist_w_unc, save_divide

set_log_level(logger, "DEBUG")


class hist_w_unc_TestCase(unittest.TestCase):
    def setUp(self):
        self.input = np.array([0, 1, 1, 3])
        self.weights = np.array([1, 2, 1, 1])
        self.n_bins = 3
        self.bin_edges = np.array([0, 1, 2, 3])

        self.hist = np.array([1, 2, 1])
        self.unc = np.sqrt(np.array([1, 2, 1]))
        self.band = self.hist - self.unc

        self.hist_normed = np.array([1, 2, 1]) / len(self.input)
        self.unc_normed = np.sqrt(np.array([1, 2, 1])) / len(self.input)
        self.band_normed = self.hist_normed - self.unc_normed

        # --- weighted cases ---
        # 3 counts in second bin due to weights
        self.hist_weighted_normed = np.array([1, 3, 1]) / np.sum(self.weights)
        # use sqrt(sum of squared weights) for error calculation
        self.unc_weighted_normed = np.sqrt(np.array([1, 2**2 + 1, 1])) / np.sum(
            self.weights
        )
        self.band_weighted_normed = self.hist_weighted_normed - self.unc_weighted_normed

        self.hist_weighted = np.array([1, 3, 1])
        # use sqrt(sum of squared weights) for error calculation
        self.unc_weighted = np.sqrt(np.array([1, 2**2 + 1, 1]))
        self.band_weighted = self.hist_weighted - self.unc_weighted

    def test_hist_w_unc_zero_case(self):
        """Test what happens if empty array is provided as input."""
        bins, hist, unc, band = hist_w_unc(
            arr=[],
            bins=[],
        )

        np.testing.assert_almost_equal(bins, [])
        np.testing.assert_almost_equal(hist, [])
        np.testing.assert_almost_equal(unc, [])
        np.testing.assert_almost_equal(band, [])

    def test_hist_w_unc_normed(self):
        """Test normalised case."""
        bins, hist, unc, band = hist_w_unc(
            arr=self.input,
            bins=self.bin_edges,
        )

        np.testing.assert_almost_equal(bins, self.bin_edges)
        np.testing.assert_almost_equal(hist, self.hist_normed)
        np.testing.assert_almost_equal(unc, self.unc_normed)
        np.testing.assert_almost_equal(band, self.band_normed)

    def test_hist_w_unc_not_normed(self):
        """Test not normalised case."""
        bins, hist, unc, band = hist_w_unc(
            arr=self.input,
            bins=self.bin_edges,
            normed=False,
        )

        np.testing.assert_almost_equal(bins, self.bin_edges)
        np.testing.assert_almost_equal(hist, self.hist)
        np.testing.assert_almost_equal(unc, self.unc)
        np.testing.assert_almost_equal(band, self.band)

    def test_histogram_weighted_normalised(self):
        """Test weighted histogram (normalised)."""

        bin_edges, hist, unc, band = hist_w_unc(
            self.input, weights=self.weights, bins=self.n_bins, normed=True
        )

        np.testing.assert_array_almost_equal(self.bin_edges, bin_edges)
        np.testing.assert_array_almost_equal(self.hist_weighted_normed, hist)
        np.testing.assert_array_almost_equal(self.unc_weighted_normed, unc)
        np.testing.assert_array_almost_equal(self.band_weighted_normed, band)

    def test_histogram_weighted_not_normalised(self):
        """Test weighted histogram (not normalised)."""

        bin_edges, hist, unc, band = hist_w_unc(
            self.input, weights=self.weights, bins=self.n_bins, normed=False
        )

        np.testing.assert_array_almost_equal(self.bin_edges, bin_edges)
        np.testing.assert_array_almost_equal(self.hist_weighted, hist)
        np.testing.assert_array_almost_equal(self.unc_weighted, unc)
        np.testing.assert_array_almost_equal(self.band_weighted, band)

    def test_range_argument_ignored(self):
        """Test if the hist_range argument is ignored when bin_edges are provided."""

        bins_range = (1, 2)

        bin_edges, hist, _, _ = hist_w_unc(
            self.input,
            bins=self.bin_edges,
            bins_range=bins_range,
            normed=False,
        )

        # check if we end up with the same bin edges anyway
        np.testing.assert_array_almost_equal(self.bin_edges, bin_edges)
        np.testing.assert_array_almost_equal(self.hist, hist)

    def test_range_argument(self):
        """Test if the hist_range argument is used when bins is an integer."""

        # we test with range from 0 to 2, with 3 bins -> [0, 0.66, 1.33, 2] exp. bins
        bins_range = (0, 2)
        bins_exp = np.array([0, 2 / 3, 1 + 1 / 3, 2])
        hist_exp = np.array([1, 2, 0])

        bin_edges, hist, _, _ = hist_w_unc(
            self.input,
            bins=self.n_bins,
            bins_range=bins_range,
            normed=False,
        )

        # check if we end up with the same bin edges anyway
        np.testing.assert_array_almost_equal(bins_exp, bin_edges)
        np.testing.assert_array_almost_equal(hist_exp, hist)

    def test_negative_weights(self):
        """Test if negative weights are properly handled."""

        values = np.array([0, 1, 2, 2, 3])
        weights = np.array([1, -1, 3, -2, 1])

        hist_exp = np.array([1, -1, 2])
        # uncertainties are the sqrt(sum of squared weights)
        unc_exp = np.sqrt(np.array([1, (-1) ** 2, 3**2 + (-2) ** 2 + 1]))

        _, hist, unc, _ = hist_w_unc(values, weights=weights, bins=3, normed=False)
        np.testing.assert_array_almost_equal(hist_exp, hist)
        np.testing.assert_array_almost_equal(unc_exp, unc)

    # TODO: Add unit tests for hist_ratio


class save_divide_TestCase(unittest.TestCase):
    def test_zero_case(self):
        steps = save_divide(np.zeros(2), np.zeros(2))
        np.testing.assert_equal(steps, np.ones(2))

    def test_ones_case(self):
        steps = save_divide(np.ones(2), np.ones(2))
        np.testing.assert_equal(steps, np.ones(2))

    def test_half_case(self):
        steps = save_divide(np.ones(2), 2 * np.ones(2))
        np.testing.assert_equal(steps, 0.5 * np.ones(2))

    def test_denominator_float(self):
        steps = save_divide(np.ones(2), 2)
        np.testing.assert_equal(steps, 0.5 * np.ones(2))

    def test_numerator_float(self):
        steps = save_divide(1, np.ones(2) * 2)
        np.testing.assert_equal(steps, 0.5 * np.ones(2))


class hist_ratio_TestCase(unittest.TestCase):
    def setUp(self):
        self.numerator = np.array([5, 3, 2, 5, 6, 2])
        self.denominator = np.array([3, 6, 2, 7, 10, 12])
        self.numerator_unc = np.array([0.5, 1, 0.3, 0.2, 0.5, 0.3])
        self.denominator_unc = np.array([1, 0.3, 2, 1, 5, 3])
        self.step = np.array([1.6666667, 1.6666667, 0.5, 1, 0.7142857, 0.6, 0.1666667])
        self.step_unc = np.array(
            [
                0.580017,
                0.580017,
                0.1685312,
                1.0111874,
                0.1059653,
                0.3041381,
                0.0485913,
            ]
        )

    def test_hist_ratio(self):
        step, step_unc = hist_ratio(
            numerator=self.numerator,
            denominator=self.denominator,
            numerator_unc=self.numerator_unc,
            denominator_unc=self.denominator_unc,
        )

        np.testing.assert_almost_equal(step, self.step)
        np.testing.assert_almost_equal(step_unc, self.step_unc)

    def test_hist_not_same_length_nominator_denominator(self):
        with self.assertRaises(AssertionError):
            _, _ = hist_ratio(
                numerator=np.ones(2),
                denominator=np.ones(3),
                numerator_unc=np.ones(3),
                denominator_unc=np.ones(3),
            )

    def test_hist_not_same_length_nomiantor_and_unc(self):
        with self.assertRaises(AssertionError):
            _, _ = hist_ratio(
                numerator=np.ones(3),
                denominator=np.ones(3),
                numerator_unc=np.ones(2),
                denominator_unc=np.ones(3),
            )

    def test_hist_not_same_length_denomiantor_and_unc(self):
        with self.assertRaises(AssertionError):
            _, _ = hist_ratio(
                numerator=np.ones(3),
                denominator=np.ones(3),
                numerator_unc=np.ones(3),
                denominator_unc=np.ones(2),
            )
