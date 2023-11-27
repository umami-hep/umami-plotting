#!/usr/bin/env python

"""Unit test script for the functions in hlplots/tagger.py."""
from __future__ import annotations

import tempfile
import unittest

import h5py
import numpy as np
import pandas as pd
from numpy.lib.recfunctions import structured_to_unstructured as s2u
from numpy.lib.recfunctions import unstructured_to_structured as u2s

from puma.hlplots import Tagger
from puma.utils import logger, set_log_level

set_log_level(logger, "DEBUG")


class TaggerBasisTestCase(unittest.TestCase):
    """Test class for the Tagger class."""

    def test_empty_string_tagger_name(self):
        """Test empty string as model name."""
        tagger = Tagger("")
        self.assertEqual(tagger.name, "")

    def test_wrong_template(self):
        """Test wrong template."""
        template_wrong = {"test": 1}
        with self.assertRaises(TypeError):
            Tagger("dummy", **template_wrong)  # pylint: disable=W0212,E1123

    def test_label_template(self):
        """Test template with label."""
        template_label = {"label": 1.5}
        tagger = Tagger("dummy", **template_label)  # pylint: disable=W0212
        self.assertEqual(tagger.label, 1.5)

    def test_n_jets(self):
        """Test if number of n_jets correctly calculated."""
        tagger = Tagger("dummy", output_flavours=["ujets", "cjets", "bjets"])
        labels = np.concatenate([np.zeros(80), np.ones(5) * 4, np.ones(15) * 5])
        tagger.labels = np.array(labels, dtype=[("HadronConeExclTruthLabelID", "i4")])
        with self.subTest():
            self.assertEqual(tagger.n_jets("ujets"), 80)
        with self.subTest():
            self.assertEqual(tagger.n_jets("cjets"), 5)
        with self.subTest():
            self.assertEqual(tagger.n_jets("bjets"), 15)


class TaggerScoreExtractionTestCase(unittest.TestCase):
    """Test extract_tagger_scores function in Tagger class."""

    def setUp(self) -> None:
        """Set up for tests."""
        self.df_dummy = pd.DataFrame(
            {
                "dummy_pc": np.zeros(10),
                "dummy_pu": np.ones(10),
                "dummy_pb": np.zeros(10),
            }
        )
        self.scores_expected = np.column_stack(
            (np.ones(10), np.zeros(10), np.zeros(10))
        )

    def test_wrong_source_type(self):
        """Test using wrong source type."""
        tagger = Tagger("dummy")
        with self.assertRaises(ValueError):
            tagger.extract_tagger_scores(self.df_dummy, source_type="dummy")

    def test_data_frame_path_no_key(self):
        """Test passing data frame path but no key."""
        tagger = Tagger("dummy")
        with self.assertRaises(ValueError):
            tagger.extract_tagger_scores(self.df_dummy, source_type="data_frame_path")

    def test_data_frame(self):
        """Test passing data frame."""
        tagger = Tagger("dummy")
        tagger.extract_tagger_scores(self.df_dummy)
        np.testing.assert_array_equal(tagger.scores, self.scores_expected)

    def test_data_frame_path(self):
        """Test passing data frame path."""
        tagger = Tagger("dummy")
        with tempfile.TemporaryDirectory() as tmp_dir:
            file_name = f"{tmp_dir}/dummy_df.h5"
            self.df_dummy.to_hdf(file_name, key="dummy_tagger")

            tagger.extract_tagger_scores(
                file_name, key="dummy_tagger", source_type="data_frame_path"
            )
        np.testing.assert_array_equal(tagger.scores, self.scores_expected)

    def test_h5_structured_numpy_path(self):
        """Test passing structured h5 path."""
        tagger = Tagger("dummy")
        with tempfile.TemporaryDirectory() as tmp_dir:
            file_name = f"{tmp_dir}/dummy_df.h5"
            with h5py.File(file_name, "w") as f_h5:
                f_h5.create_dataset(
                    data=self.df_dummy.to_records(), name="dummy_tagger"
                )

            tagger.extract_tagger_scores(
                file_name, key="dummy_tagger", source_type="h5_file"
            )
        np.testing.assert_array_equal(s2u(tagger.scores), self.scores_expected)

    def test_structured_array(self):
        """Test passing structured numpy array."""
        tagger = Tagger("dummy")
        tagger.extract_tagger_scores(
            self.df_dummy.to_records(),
            key="dummy_tagger",
            source_type="structured_array",
        )
        np.testing.assert_array_equal(s2u(tagger.scores), self.scores_expected)


class TaggerTestCase(unittest.TestCase):
    """Test class for the Tagger class."""

    def setUp(self) -> None:
        """Set up for tests."""
        scores = np.column_stack((np.ones(10), np.ones(10), np.ones(10)))
        self.scores = u2s(
            scores,
            dtype=[("dummy_pu", "f4"), ("dummy_pc", "f4"), ("dummy_pb", "f4")],
        )

    def test_errors(self):
        tagger = Tagger("dummy")
        tagger.scores = self.scores
        with self.assertRaises(ValueError):
            tagger.discriminant("hbb", fx=0.1)
        with self.assertRaises(ValueError):
            tagger.discriminant("ujets")

    def test_disc_b_calc_no_fc(self):
        """Test b-disc calculation w/o f_c provided."""
        tagger = Tagger("dummy")
        tagger.scores = self.scores
        with self.assertRaises(TypeError):
            tagger.discriminant("bjets")

    def test_disc_b_calc(self):
        """Test b-disc calculation."""
        tagger = Tagger("dummy", f_c=0.5)
        tagger.scores = self.scores
        discs = tagger.discriminant("bjets")
        np.testing.assert_array_equal(discs, np.zeros(10))

    def test_disc_c_calc_no_fb(self):
        """Test c-disc calculation w/o f_c provided."""
        tagger = Tagger("dummy")
        tagger.scores = self.scores
        with self.assertRaises(TypeError):
            tagger.discriminant("cjets")

    def test_disc_c_calc(self):
        """Test c-disc calculation."""
        tagger = Tagger("dummy", f_b=0.5)
        tagger.scores = self.scores
        discs = tagger.discriminant("cjets")
        np.testing.assert_array_equal(discs, np.zeros(10))

    def test_disc_hbb_calc(self):
        """Test hbb-disc calculation."""
        from ftag import Flavours as F

        tagger = Tagger(
            "dummy", output_flavours=[F["hbb"], F["hcc"], F["tqqb"], F["qcd"]]
        )
        tagger.scores = u2s(
            np.column_stack((np.ones(10), np.ones(10), np.ones(10), np.ones(10))),
            dtype=[
                ("dummy_phbb", "f4"),
                ("dummy_phcc", "f4"),
                ("dummy_ptqqb", "f4"),
                ("dummy_pqcd", "f4"),
            ],
        )
        discs = tagger.discriminant("hbb")
        np.testing.assert_array_equal(discs, np.ones(10))
