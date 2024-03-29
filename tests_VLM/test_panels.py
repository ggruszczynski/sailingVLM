"""
    Unit tests_LLT_optimizer of the Panel class and its methods

"""

import numpy as np
from numpy.testing import assert_almost_equal

from Solver.Panel import Panel
from unittest import TestCase


class TestPanels(TestCase):
    def setUp(self):
        p1 = np.array([10., 0., 0.])
        p2 = np.array([0., 0., 0.])
        p3 = np.array([0., 10., 0.])
        p4 = np.array([10., 10., 0.])

        p1next = p1 + (p1 - p2)
        p4next = p4 + (p4 - p3)

        self.points = [p1, p2, p3, p4]
        self.panel = Panel(*self.points, gamma_orientation=1, p1next=p1next, p4next=p4next)
        self.assertTrue(self.panel._are_points_coplanar())

    def test_area(self):
        expected_area = 100.0
        assert_almost_equal(self.panel.area, expected_area)

    def test_pressure(self):
        self.panel.force_xyz = np.array([3., 2., 1.])
        self.panel.calc_pressure()
        assert_almost_equal(self.panel.pressure, 0.01)

    def test_get_ctr_point_postion(self):
        ctr_point = self.panel.ctr_point_position
        expected_ctr_point = [7.5, 5, 0]

        assert_almost_equal(expected_ctr_point, ctr_point)

        points2 = [np.array([ 8., 2., 0]), np.array([0., 0., 0]),
                   np.array([-2., 6., 0]), np.array([6., 7., 0])]

        panel2 = Panel(*points2)
        ctr_point2 = panel2.ctr_point_position
        expected_ctr_point2 = [5., 4.125, 0]
        assert_almost_equal(expected_ctr_point2, ctr_point2)

    def test_get_cp_postion(self):
        cp = self.panel.cp_position
        expected_ctr_point = [2.5, 5, 0]

        assert_almost_equal(expected_ctr_point, cp)

        points2 = [np.array([ 8., 2., 0]), np.array([0., 0., 0]),
                   np.array([-2., 6., 0]), np.array([6., 7., 0])]

        panel2 = Panel(*points2)
        cp2 = panel2.cp_position
        expected_cp2 = [1., 3.375, 0]
        assert_almost_equal(expected_cp2, cp2)

    def test_get_vortex_ring_position(self):
        vortex_ring_position = self.panel.get_vortex_ring_position()
        expected_vortex_riing_position = [[12.5, 0., 0.],
                                          [2.5, 0., 0.],
                                          [2.5, 10., 0.],
                                          [12.5, 10., 0.]]

        assert_almost_equal(expected_vortex_riing_position, vortex_ring_position)

    def test_get_vortex_ring_induced_velocity(self):
        ctr_p = self.panel.ctr_point_position
        dummy_velocity = None
        v_ind = self.panel.get_induced_velocity(ctr_p, dummy_velocity)
        v_ind_expected = [0, 0, -0.09003163161571061]

        assert_almost_equal(v_ind, v_ind_expected)

    def test_panel_is_not_plane(self):
        points = [np.array([10, 0, 0]), np.array([0, 0, 666]),
                  np.array([0, 10, 0]), np.array([10, 10, 0])]

        panel = Panel(*points)
        self.assertFalse(panel._are_points_coplanar())
