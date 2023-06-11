import numpy as np
import pandas as pd

from abc import abstractmethod, ABC
from Rotations.geometry_calc import rotation_matrix
from Solver import Panel
from Rotations.CSYS_transformations import CSYS_transformations
from Solver.mesher import make_panels_from_le_te_points, make_panels_from_le_points_and_chords, make_panels_from_mesh_spanwise, make_airfoil_mesh
from typing import List
from YachtGeometry.SailBaseGeometry import SailBaseGeometry

from Rotations.geometry_calc import rotate_points_around_arbitrary_axis


class CamberedSailGeometry(SailBaseGeometry, ABC):
    def __init__(self, head_mounting: np.array, tack_mounting: np.array,
                 csys_transformations: CSYS_transformations,
                 n_spanwise=10, n_chordwise=1, chords=None,
                 initial_sail_twist_deg=None, name=None, LLT_twist=None,
                 interpolated_camber=None, interpolated_max_camber_distance_from_luff=None
                 ):

        super(CamberedSailGeometry, self).__init__(head_mounting, tack_mounting,
                                           csys_transformations,
                                           n_spanwise, n_chordwise, name)

        chords_vec = np.array([chords, np.zeros(len(chords)), np.zeros(len(chords))])
        chords_vec = chords_vec.transpose()
        fchords_vec = np.flip(chords_vec, axis=0)

        # state "zero" (i.e. yacht in an upright position, without any rotations)
        mesh = make_airfoil_mesh([self.le_SW, self.le_NW],
                                 [self._n_chordwise, self._n_spanwise], chords_vec,
                                 interpolated_max_camber_distance_from_luff, interpolated_camber)
        sh0, sh1, sh2 = mesh.shape
        mesh = mesh.reshape(sh0*sh1, sh2)
        mesh_underwater = make_airfoil_mesh([self.le_SW_underwater, self.le_NW_underwater],
                                            [self._n_chordwise, self._n_spanwise], fchords_vec,
                                            interpolated_max_camber_distance_from_luff, interpolated_camber).reshape(sh0 * sh1, sh2)
        mesh_underwater = mesh_underwater.reshape(sh0*sh1, sh2)

        # rotation
        # mesh:  le NW p1 p2 ... te
        #       le p1 p2 ... te
        #       le p1 p2 ... te
        #       le SW p2 ... te
        rmesh = np.array([self.csys_transformations.rotate_point_around_origin_with_mirror(point) for point in mesh])
        rmesh_underwater = np.array([self.csys_transformations.rotate_point_around_origin_with_mirror(point) for point in mesh_underwater])

        mesh = rmesh
        mesh_underwater = rmesh_underwater

        ## twist
        print(f"Applying initial_sail_twist_deg to {self.name}, mode: {LLT_twist}")
        twist_dict = {
            'sheeting_angle_const': np.full(len(initial_sail_twist_deg), np.min(initial_sail_twist_deg)),
            'average_const': np.full(len(initial_sail_twist_deg), np.average(initial_sail_twist_deg)),
            'real_twist': initial_sail_twist_deg
        }
        sail_twist_deg = twist_dict[LLT_twist]
        sail_twist_deg = np.hstack([sail_twist_deg] * (sh1))

        sail_twist_deg = sail_twist_deg.reshape(sh1, sh0).transpose().flatten()

        p2 = mesh[::sh1][-1]
        p1 = mesh[::sh1][0]
        trmesh = self.rotate_points_around_le(mesh, p1, p2, sail_twist_deg)
        # check if points on forestay are not rotated (they are on axis of rotation)
        np.testing.assert_almost_equal(trmesh[::sh1], mesh[::sh1])

        p2_u = mesh_underwater[::sh1][-1]
        p1_u = mesh_underwater[::sh1][0]
        trmesh_underwater = self.rotate_points_around_le(mesh_underwater, p1_u, p2_u, np.flip(sail_twist_deg))
        # check if points on forestay are not rotated (they are on axis of rotation)
        np.testing.assert_almost_equal(trmesh_underwater[::sh1], mesh_underwater[::sh1])

        mesh = trmesh
        mesh_underwater = trmesh_underwater

        # come back to original shape
        mesh = mesh.reshape(sh0, sh1, sh2)
        mesh = np.swapaxes(mesh, 0, 1)
        mesh_underwater = mesh_underwater.reshape(sh0, sh1, sh2)
        mesh_underwater = np.swapaxes(mesh_underwater, 0, 1)

        panels = make_panels_from_mesh_spanwise(mesh, gamma_orientation=-1)
        panels_mirror = make_panels_from_mesh_spanwise(mesh_underwater, gamma_orientation=-1)

        # https://stackoverflow.com/questions/33356442/when-should-i-use-hstack-vstack-vs-append-vs-concatenate-vs-column-stack
        # self.__panels = np.vstack((panels, panels_mirror))
        self.__panels = np.hstack((panels_mirror, panels))  # original version
        self.__panels1D = self.__panels.flatten()
        self.__spans = np.array([panel.get_panel_span_at_cp() for panel in self.panels1d])

    def rotate_points_around_le(self, points, p1, p2, sail_twist_deg_vec):
        rotated_points = points

        # if all elements are the same -> only once do the calculation -> less matrix multiplications
        if np.all(sail_twist_deg_vec == sail_twist_deg_vec[0]):
            rotated_points = rotate_points_around_arbitrary_axis(points, p1, p2, np.deg2rad(sail_twist_deg_vec[0]))
        else:
            # squezee removes "1" dimention
            # rotate_points_around_arbitrary_axis needs [[a, b, c]] or [[a, b, c], [e, f, g], ...]
           rotated_points = np.array(
                [np.squeeze(rotate_points_around_arbitrary_axis(np.array([point]), p1, p2, np.deg2rad(deg))) for
                 point, deg in zip(points, sail_twist_deg_vec)])

        return rotated_points
    def get_cp_points_upright(self):
        cp_points = self.get_cp_points1d()
        cp_straight_yacht = np.array([self.csys_transformations.reverse_point_rotations_around_origin_with_mirror(p) for p in cp_points])
        return cp_straight_yacht

    def sail_cp_to_girths(self):
        sail_cp_straight_yacht = self.get_cp_points_upright()
        tack_mounting = self.tack_mounting
        y = sail_cp_straight_yacht[:, 2]
        y_as_girths = (y - tack_mounting[2]) / (max(y) - tack_mounting[2])
        return y_as_girths

    @property
    def spans(self):
        return self.__spans

    @property
    def panels1d(self) -> np.array([Panel]):
        return self.__panels1D

    @property
    def panels(self) -> np.array([Panel]):
        return self.__panels


