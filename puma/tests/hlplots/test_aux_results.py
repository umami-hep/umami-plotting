#!/usr/bin/env python
"""Unit test script for the functions in hlplots/aux_results.py."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from puma.hlplots import VtxResults
from puma.hlplots.tagger import Tagger
from puma.utils import logger, set_log_level, get_dummy_tagger_aux
from puma.utils.vertexing import calculate_vertex_metrics

set_log_level(logger, "DEBUG")


class VtxResultsTestCase(unittest.TestCase):
    """Test class for the VtxResults class."""

    def test_add_taggers_from_file(self):
        """Test for Results.add_taggers_from_file function."""
        np.random.default_rng(seed=16)
        fname = get_dummy_tagger_aux()[0]
        results = VtxResults(signal="bjets", sample="test")
        taggers = [Tagger("GN2")]
        results.add_taggers_from_file(
            taggers,
            fname,
        )
        self.assertEqual(list(results.taggers.values()), taggers)

    def test_add_taggers_with_cuts(self):
        np.random.default_rng(seed=16)
        fname = get_dummy_tagger_aux()[0]
        cuts = [("eta", ">", 0)]
        tagger_cuts = [("pt", ">", 20)]
        results = VtxResults(signal="bjets", sample="test")
        taggers = [Tagger("GN2", cuts=tagger_cuts)]
        results.add_taggers_from_file(
            taggers,
            fname,
            cuts=cuts,
        )
        self.assertEqual(list(results.taggers.values()), taggers)


class VtxResultsPlotsTestCase(unittest.TestCase):
    """Test class for the VtxResults class running plots."""

    def setUp(self) -> None:
        """Set up for unit tests."""
        f = get_dummy_tagger_aux()[1]
        dummy_tagger_1 = Tagger("GN2")
        dummy_tagger_1.labels = np.array(
            f["jets"]["HadronConeExclTruthLabelID"],
            dtype=[("HadronConeExclTruthLabelID", "i4")],
        )
        dummy_tagger_1.aux_metrics = calculate_vertex_metrics(
            f["tracks"]["ftagTruthVertexIndex"],
            f["tracks"]["GN2_VertexIndex"],
        )
        dummy_tagger_1.perf_var = f["jets"]["pt"]
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

    def test_plot_var_vtx_eff(self):
        """Test that png file is being created."""
        self.dummy_tagger_1.reference = True
        with tempfile.TemporaryDirectory() as tmp_file:
            vtxresults = VtxResults(signal="bjets", sample="test", output_dir=tmp_file)
            vtxresults.add(self.dummy_tagger_1)
            vtxresults.plot_var_vtx_eff()
            self.assertIsFile(vtxresults.get_filename("vtx_eff_vs_pt"))

    def test_plot_var_vtx_fr(self):
        """Test that png file is being created."""
        self.dummy_tagger_1.reference = True
        with tempfile.TemporaryDirectory() as tmp_file:
            vtxresults = VtxResults(signal="bjets", sample="test", output_dir=tmp_file)
            vtxresults.add(self.dummy_tagger_1)
            vtxresults.plot_var_vtx_fr()
            self.assertIsFile(vtxresults.get_filename("vtx_fr_vs_pt"))

    def test_plot_var_vtx_trk_eff(self):
        """Test that png file is being created."""
        self.dummy_tagger_1.reference = True
        with tempfile.TemporaryDirectory() as tmp_file:
            vtxresults = VtxResults(signal="bjets", sample="test", output_dir=tmp_file)
            vtxresults.add(self.dummy_tagger_1)
            vtxresults.plot_var_vtx_trk_eff()
            self.assertIsFile(vtxresults.get_filename("vtx_trk_eff_vs_pt"))

    def test_plot_var_vtx_trk_fr(self):
        """Test that png file is being created."""
        self.dummy_tagger_1.reference = True
        with tempfile.TemporaryDirectory() as tmp_file:
            vtxresults = VtxResults(signal="bjets", sample="test", output_dir=tmp_file)
            vtxresults.add(self.dummy_tagger_1)
            vtxresults.plot_var_vtx_trk_fr()
            self.assertIsFile(vtxresults.get_filename("vtx_trk_fr_vs_pt"))
