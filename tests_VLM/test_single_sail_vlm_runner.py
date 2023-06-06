import shutil

from YachtGeometry.SailFactory import SailFactory
from YachtGeometry.SailSet import SailSet
from Rotations.CSYS_transformations import CSYS_transformations
from Solver.Interpolator import Interpolator
from Inlet.InletConditions import InletConditions
from Inlet.Winds import ExpWindProfile, FlatWindProfile

from Solver.vlm_solver import is_no_flux_BC_satisfied

from Solver.vlm_solver import calc_circulation
from ResultsContainers.InviscidFlowResults import prepare_inviscid_flow_results_llt, prepare_inviscid_flow_results_vlm
from Solver.vlm_solver import calculate_app_fs
from Solver.coeff_formulas import get_CL_CD_free_wing
import pandas as pd
from unittest import TestCase


from tests_VLM.InputFiles.swept_sail import *
# np.set_printoptions(precision=3, suppress=True)
from numpy.testing import assert_almost_equal

class TestVLM_Solver(TestCase):
    def setUp(self):
        self.interpolator = Interpolator(interpolation_type)

        self.csys_transformations = CSYS_transformations(
            heel_deg, leeway_deg,
            v_from_original_xyz_2_reference_csys_xyz=reference_level_for_moments)

        self.wind = FlatWindProfile(alpha_true_wind_deg, tws_ref, SOG_yacht)

    def _prepare_sail_set(self, n_spanwise, n_chordwise):

        sail_factory = SailFactory(n_spanwise=n_spanwise, n_chordwise=n_chordwise,
                                   csys_transformations=self.csys_transformations, rake_deg=rake_deg,
                                   sheer_above_waterline=sheer_above_waterline)

        main_sail_geometry = sail_factory.make_main_sail(
            main_sail_luff=main_sail_luff,
            boom_above_sheer=boom_above_sheer,
            main_sail_chords=self.interpolator.interpolate_girths(main_sail_girths, main_sail_chords, n_spanwise + 1),
            sail_twist_deg=self.interpolator.interpolate_girths(main_sail_girths, main_sail_centerline_twist_deg, n_spanwise + 1),
            LLT_twist=LLT_twist,
            interpolated_camber=self.interpolator.interpolate_girths(main_sail_girths, main_sail_camber, n_spanwise + 1),
            interpolated_distance_from_LE=self.interpolator.interpolate_girths(main_sail_girths,
                                                                               main_sail_max_camber_distance_from_luff,
                                                                               n_spanwise + 1))
        return main_sail_geometry

    def test_calc_forces_and_moments_single_sail(self):
        """
        Benchmark from "Aerodynamics for Engineers" 5th edition, John. J. Bertin, Pearson 2009
        p 364
        Example 7.2
        """

        main_sail_geometry = self._prepare_sail_set(n_spanwise=4, n_chordwise=1)
        sail_set = SailSet([main_sail_geometry])


        inlet_condition = InletConditions(self.wind, rho=rho, panels1D=sail_set.panels1d)
        gamma_magnitude, v_ind_coeff = calc_circulation(inlet_condition.V_app_infs, sail_set.panels1d)
        V_induced_at_ctrl_p, V_app_fs_at_ctrl_p = calculate_app_fs(inlet_condition, v_ind_coeff, gamma_magnitude)
        assert is_no_flux_BC_satisfied(V_app_fs_at_ctrl_p, sail_set.panels1d)
        inviscid_flow_results = prepare_inviscid_flow_results_vlm(gamma_magnitude,
                                                                      sail_set, inlet_condition,
                                                                      self.csys_transformations)

        q = 0.5 * rho * (np.linalg.norm(V_inf) ** 2)
        wing_span = 2 * half_wing_span
        S = chord * wing_span
        AR = wing_span / chord

        sail_mirror_multiplier = 2

        from Rotations.geometry_calc import rotation_matrix
        F_xyz_total = sail_mirror_multiplier*inviscid_flow_results.F_xyz_total
        A = rotation_matrix([0, 0, 1], np.deg2rad(-AoA_deg)) # rotate the wind instead of the lifting surface. Be aware the sail_twist rotates around the luff (leading edge)
        F = np.dot(A, F_xyz_total)  # By definition, the lift force is perpendicular to V_inf

        CX_vlm = F[0] / (q * S)
        CY_vlm = F[1] / (q * S)
        CZ_vlm = F[2] / (q * S)

        a_VLM = CY_vlm / np.deg2rad(AoA_deg)
        # CL_expected, CD_ind_expected, a_expected = get_CL_CD_free_wing(AR, AoA_deg, sweep_half_chord_deg=sweep_angle_deg)

        assert_almost_equal(sail_set.panels1d[7].get_vortex_ring_position()[1], np.array([0.425, 0., 0.375]), decimal=3)
        assert_almost_equal(sail_set.panels1d[7].get_vortex_ring_position()[2], np.array([0.550, 0., 0.500]), decimal=3)
        assert_almost_equal(sail_set.get_ctr_points1d()[7], np.array([0.5875, 0., 0.4375]), decimal=4)

        assert_almost_equal(inviscid_flow_results.gamma_magnitude[4:], np.array([0.0239377, 0.02521229, 0.02511784, 0.02187719]), decimal=6)

        assert_almost_equal(F, np.array([0.00031786, 0.02402581, 0.00135887]), decimal=6)
        assert_almost_equal(CY_vlm, 0.240258, decimal=6)
        assert_almost_equal(a_VLM, 3.441443, decimal=6)

        def get_gamma_reference(b, V, AoA_deg):
            # Eq 7.48, p368 from "Aerodynamics for Engineers" 5th edition, John. J. Bertin, Pearson 2009
            return 4*np.pi*b*V*np.deg2rad(AoA_deg)*np.array([0.0273, 0.0287, 0.0286, 0.0250])

        gamma_ref = get_gamma_reference(wing_span, V_inf, AoA_deg)
        dy = wing_span/(2*len(gamma_ref))
        L_ref = 2*rho*V_inf*V_inf*np.sum(gamma_ref*dy)    #eq 7.50b
        CL_ref = L_ref/(q*S)
        a_ref = CL_ref/np.deg2rad(AoA_deg)

        L_ref2 = rho * V_inf * V_inf * wing_span * wing_span * np.pi * np.deg2rad(AoA_deg) * 0.1096
        CL_ref2 = 1.096*np.pi*np.deg2rad(AoA_deg)

        assert_almost_equal(L_ref, L_ref2, decimal=6)
        assert_almost_equal(L_ref, F[1], decimal=4)
        assert_almost_equal(CL_ref, CL_ref2, decimal=6)
        assert_almost_equal(CL_ref, 0.240379, decimal=6)
        assert_almost_equal(a_ref, 3.443185, decimal=6)
