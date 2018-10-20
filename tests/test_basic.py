from context import sample
import pytest as pt
import numpy as np

# %% orbital_analyses/Transform_State.py


class TestTransformState(object):
    def test_jd2gd(self):
        from sample.orbital_analyses.Transform_State import JD2Gregorian
        jd = 2440422.622396  # First man on the moon
        gd = np.matrix([[1969], [7], [20], [2], [56], [15.0]])
        assert np.all(JD2Gregorian(jd) == pt.approx(gd))

    def test_gd2jd(self):
        from sample.orbital_analyses.Transform_State import Gregorian2JD
        jd = 2440422.622396  # First man on the moon
        gd = np.matrix([[1969], [7], [20], [2], [56], [15.0]])
        assert Gregorian2JD(gd) == pt.approx(jd)

    def test_date_equivalence(self):
        from sample.orbital_analyses.Transform_State import JD2Gregorian
        from sample.orbital_analyses.Transform_State import Gregorian2JD
        jd = 2440422.622396  # First man on the moon
        gd = np.matrix([[1969], [7], [20], [2], [56], [15.0]])
        gd_new = JD2Gregorian(jd)
        jd_new = Gregorian2JD(gd)
        assert Gregorian2JD(gd_new) == pt.approx(jd)
        assert np.all(JD2Gregorian(jd_new) == pt.approx(gd))

