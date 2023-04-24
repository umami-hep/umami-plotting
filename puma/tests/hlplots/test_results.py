#!/usr/bin/env python
"""Unit test script for the functions in hlplots/tagger.py."""

import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
from ftag import get_mock_file

from puma.hlplots import Results
from puma.hlplots.tagger import Tagger
from puma.utils import (
    get_dummy_2_taggers,
    logger,
    set_log_level,
)

set_log_level(logger, "DEBUG")


class ResultsTestCase(unittest.TestCase):
    """Test class for the Results class."""

    def test_add_duplicated(self):
        """Test empty string as model name."""
        dummy_tagger_1 = Tagger("dummy")
        dummy_tagger_2 = Tagger("dummy")
        results = Results(signal="bjets", sample="test")
        results.add(dummy_tagger_1)
        with self.assertRaises(KeyError):
            results.add(dummy_tagger_2)

    def test_add_2_taggers(self):
        """Test empty string as model name."""
        dummy_tagger_1 = Tagger("dummy")
        dummy_tagger_2 = Tagger("dummy_2")
        results = Results(signal="bjets", sample="test")
        results.add(dummy_tagger_1)
        results.add(dummy_tagger_2)
        self.assertEqual(
            list(results.taggers.keys()),
            ["dummy (dummy)", "dummy_2 (dummy_2)"],  # pylint: disable=W0212
        )

    def test_get_taggers(self):
        """Test empty string as model name."""
        dummy_tagger_1 = Tagger("dummy")
        dummy_tagger_2 = Tagger("dummy_2")
        results = Results(signal="bjets", sample="test")
        results.add(dummy_tagger_1)
        results.add(dummy_tagger_2)
        retrieved_dummy_tagger_2 = results["dummy_2 (dummy_2)"]
        self.assertEqual(retrieved_dummy_tagger_2.name, dummy_tagger_2.name)

    def test_add_taggers_from_file(self):
        """Test for Results.add_taggers_from_file function."""
        tmp_dir = tempfile.TemporaryDirectory()  # pylint: disable=R1732
        rng = np.random.default_rng(seed=16)
        with h5py.File(f"{tmp_dir.name}/test.h5", "w") as file:
            data = get_dummy_2_taggers()
            data["pt"] = rng.random(len(data))
            file.create_dataset("jets", data=data.to_records())
        results = Results(signal="bjets", sample="test")
        taggers = [Tagger("rnnip")]
        results.add_taggers_from_file(
            taggers, f"{tmp_dir.name}/test.h5", perf_var=data["pt"]
        )
        self.assertEqual(list(results.taggers.values()), taggers)


class ResultsPlotsTestCase(unittest.TestCase):
    """Test class for the Results class running plots."""

    def setUp(self) -> None:
        """Set up for unit tests."""
        f = get_mock_file()[1]
        dummy_tagger_1 = Tagger("MockTagger")
        dummy_tagger_1.labels = np.array(
            f["jets"]["HadronConeExclTruthLabelID"],
            dtype=[("HadronConeExclTruthLabelID", "i4")],
        )
        dummy_tagger_1.scores = f["jets"]
        dummy_tagger_1.label = "dummy tagger"
        self.dummy_tagger_1 = dummy_tagger_1

    def assertIsFile(self, path: str):
        """Check for file to exist.
        Taken from https://stackoverflow.com/a/59198749/10896585
        Parameters
        ----------
        path : str
            Path to file.

        Raises
        ------
        AssertionError
            if file does not exist
        """
        if not Path(path).resolve().is_file():
            raise AssertionError(f"File does not exist: {path}")

    def test_plot_roc_bjets(self):
        """Test that png file is being created."""
        self.dummy_tagger_1.reference = True
        self.dummy_tagger_1.f_c = 0.05
        with tempfile.TemporaryDirectory() as tmp_file:
            results = Results(signal="bjets", sample="test", output_dir=tmp_file)
            results.add(self.dummy_tagger_1)
            results.sig_eff = np.linspace(0.6, 0.95, 20)
            results.plot_rocs()
            self.assertIsFile(results.get_filename("roc"))

    def test_plot_roc_cjets(self):
        """Test that png file is being created."""
        self.dummy_tagger_1.reference = True
        self.dummy_tagger_1.f_b = 0.05
        with tempfile.TemporaryDirectory() as tmp_file:
            results = Results(signal="cjets", sample="test", output_dir=tmp_file)
            results.add(self.dummy_tagger_1)
            results.sig_eff = np.linspace(0.2, 0.95, 20)
            results.plot_rocs()
            self.assertIsFile(results.get_filename("roc"))

    def test_plot_var_perf_bjets(self):
        """Test that png file is being created."""
        self.dummy_tagger_1.reference = True
        self.dummy_tagger_1.f_c = 0.05
        self.dummy_tagger_1.disc_cut = 2
        rng = np.random.default_rng(seed=16)
        self.dummy_tagger_1.perf_var = rng.exponential(
            100, size=len(self.dummy_tagger_1.scores)
        )
        with tempfile.TemporaryDirectory() as tmp_file:
            results = Results(signal="bjets", sample="test", output_dir=tmp_file)
            results.add(self.dummy_tagger_1)
            results.plot_var_perf(
                bins=[20, 30, 40, 60, 85, 110, 140, 175, 250],
            )

            self.assertIsFile(Path(tmp_file) / "test_bjets_profile_fixed_bjets_eff.png")
            self.assertIsFile(Path(tmp_file) / "test_bjets_profile_fixed_cjets_rej.png")
            self.assertIsFile(Path(tmp_file) / "test_bjets_profile_fixed_ujets_rej.png")

    def test_plot_var_perf_cjets(self):
        """Test that png file is being created."""
        self.dummy_tagger_1.reference = True
        self.dummy_tagger_1.f_b = 0.05
        self.dummy_tagger_1.working_point = 0.5
        rng = np.random.default_rng(seed=16)
        self.dummy_tagger_1.perf_var = rng.exponential(
            100, size=len(self.dummy_tagger_1.scores)
        )
        with tempfile.TemporaryDirectory() as tmp_file:
            results = Results(signal="cjets", sample="test", output_dir=tmp_file)
            results.add(self.dummy_tagger_1)
            results.plot_var_perf(
                h_line=self.dummy_tagger_1.working_point,
                bins=[20, 30, 40, 60, 85, 110, 140, 175, 250],
            )
            self.assertIsFile(Path(tmp_file) / "test_cjets_profile_fixed_cjets_eff.png")
            self.assertIsFile(Path(tmp_file) / "test_cjets_profile_fixed_bjets_rej.png")
            self.assertIsFile(Path(tmp_file) / "test_cjets_profile_fixed_ujets_rej.png")

    def test_plot_discs_bjets(self):
        """Test that png file is being created."""
        self.dummy_tagger_1.reference = True
        self.dummy_tagger_1.f_c = 0.05
        with tempfile.TemporaryDirectory() as tmp_file:
            results = Results(signal="bjets", sample="test", output_dir=tmp_file)
            results.add(self.dummy_tagger_1)
            results.plot_discs()
            self.assertIsFile(results.get_filename("disc"))

    def test_plot_discs_cjets(self):
        """Test that png file is being created."""
        self.dummy_tagger_1.reference = True
        self.dummy_tagger_1.f_b = 0.05
        with tempfile.TemporaryDirectory() as tmp_file:
            results = Results(signal="cjets", sample="test", output_dir=tmp_file)
            results.add(self.dummy_tagger_1)
            results.plot_discs()
            self.assertIsFile(results.get_filename("disc"))
