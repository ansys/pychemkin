"""Tests for Pychemkin run environment."""

from pathlib import Path

import pytest

import ansys.chemkin.core as ck


class TestClassInstallation:
    """verify Ansys Chemkin and PyChemkin installations."""

    def check_ck_install(self):
        """Check Ansys Chemkin installation location."""
        return Path(ck.ansys_dir()).is_dir()

    def test_installations(self):
        """Check proper installations."""
        assert self.check_ck_install

    def test_minimum_version(self):
        """Check minimum version required to run PyChemkin."""
        minimum_version = 252
        assert ck.chemkin_version() >= minimum_version

    def test_chemistry_metadata_freeze(self):
        """Check chemistry metadata protection."""
        """
        Allow only Chemistry label and real-gas status changes after preprocessing.
        """
        chemistry = ck.Chemistry()
        chemistry._freeze_metadata()

        chemistry.label = "updated label"
        chemistry.userealgas = True

        assert chemistry.label == "updated label"
        assert chemistry.userealgas
        with pytest.raises(AttributeError, match="metadata is immutable"):
            chemistry.chemfile = "another-chemistry.inp"
        with pytest.raises(ValueError, match="read-only"):
            chemistry._wt[0] = 1.0
        with pytest.raises(AttributeError):
            chemistry.ksymbol.append("new-species")
        with pytest.raises(TypeError):
            chemistry._gas_species["new-species"] = 0
