#!/usr/bin/env python

"""
Unit test script for the functions in plot_base.py
"""

import unittest

from puma.utils import logger, set_log_level
from puma import PlotObject

set_log_level(logger, "DEBUG")


class plot_object_TestCase(unittest.TestCase):
    """Test class for the puma.PlotObject dataclass."""

    def test_only_one_input_figsize(self):
        with self.assertRaises(ValueError):
            PlotObject(figsize=1)

    def test_only_tuple_three_inputs_figsize(self):
        with self.assertRaises(ValueError):
            PlotObject(figsize=(1, 2, 3))

    def test_tuple_input_figsize(self):
        figsize = PlotObject(figsize=(1, 2)).figsize
        self.assertEqual(figsize, (1, 2))

    def test_list_input_figsize(self):
        figsize = PlotObject(figsize=[1, 2]).figsize
        self.assertEqual(figsize, (1, 2))

    def test_list_input_wrong_len_figsize(self):
        with self.assertRaises(ValueError):
            PlotObject(figsize=[1, 2, 3])

    def test_wrong_n_ratio_panels(self):
        with self.assertRaises(ValueError):
            PlotObject(n_ratio_panels=5)
