"""Tests for PyChemkin cleanup behavior."""

import ansys.chemkin.core.chemistry as chemistry


class TestClassCleanup:
    """Test cleanup behavior."""

    def test_done_swallows_kinexit_oserror(self, monkeypatch):
        """An OS error from native shutdown does not escape ``done``."""

        class FailingChemkin:
            def KINExit(self):
                raise OSError("native shutdown failure")

        monkeypatch.setattr(chemistry, "chemkin_version", lambda: 271)
        monkeypatch.setattr(chemistry.ck_wrapper, "chemkin", FailingChemkin())
        monkeypatch.setattr(chemistry, "clear_hints", lambda: None)
        warnings = []
        monkeypatch.setattr(chemistry, "_log_warning_message", warnings.append)

        chemistry._chemset_identifiers[:] = ["test"]
        chemistry._CKInitialized["test"] = True

        chemistry.done()

        assert chemistry._chemset_identifiers == []
        assert chemistry._CKInitialized == {}
        assert warnings
        assert "KINExit() reported an OS error" in " ".join(warnings[0])