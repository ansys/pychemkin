import os

import ansys.chemkin as ck

class TestClassInstallation:
    """
    verify Ansys Chemkin and PyChemkin installations.
    """
    def check_CK_install(self):
        """
        Check Ansys Chemkin installation location.
        """
        return os.path.isdir(ck.chemkin_dir)


    def check_PyCK_install(self):
        """
        Check PyChemkin module installtion.
        """
        return os.path.isdir(ck.pychemkin_dir)


    def test_installations(self):
        """
        Check proper installations.
        """
        assert self.check_CK_install() and self.check_PyCK_install()


    def test_minimum_version(self):
        """
        Check minimum version required to run PyChemkin.
        """
        minimum_version = 252
        assert ck.ansys_version >= minimum_version
