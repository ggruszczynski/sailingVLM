import numpy as np
from numpy.testing import assert_almost_equal
from unittest import TestCase

from Solver.coeff_formulas import get_CL_CD_free_wing


class TestMesher(TestCase):
    def test_get_CL_CD_from_coeff(self):
        AR = 20
        AoA_deg = 10
        sweep_half_chord_deg = 5
        CL_expected, CD_ind_expected, a = get_CL_CD_free_wing(AR, AoA_deg, sweep_half_chord_deg)

        assert_almost_equal(CL_expected, 0.9890278418785605, decimal=6)
        assert_almost_equal(CD_ind_expected, 0.019460194634344816, decimal=6)
