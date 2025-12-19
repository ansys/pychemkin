"""Tests for Pychemkin run environment."""

from pathlib import Path

import ansys.chemkin.core as ck


class TestClassInstallation:
    """verify Ansys Chemkin and PyChemkin installations."""

    def check_ck_install(self):
        """Check Ansys Chemkin installation location."""
        return Path(ck.chemkin_dir).is_dir()

    def check_pyck_install(self):
        """Check PyChemkin module installation."""
        return Path(ck.pychemkin_dir).is_dir()

    def test_installations(self):
        """Check proper installations."""
        assert self.check_ck_install() and self.check_pyck_install()

    def test_minimum_version(self):
        """Check minimum version required to run PyChemkin."""
        minimum_version = 252
        assert ck.ansys_version >= minimum_version
