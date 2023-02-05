import numpy as np
import pandas as pd

from abc import abstractmethod, ABC
from Rotations.geometry_calc import rotation_matrix
from Solver import Panel
from Rotations.CSYS_transformations import CSYS_transformations
from Solver.mesher import make_panels_from_le_te_points, make_panels_from_le_points_and_chords
from typing import List
from YachtGeometry.SailBaseGeometry import SailBaseGeometry


class SailGeometry(SailBaseGeometry, ABC):
    def __init__(self, head_mounting: np.array, tack_mounting: np.array,
                 csys_transformations: CSYS_transformations,
                 n_spanwise=10, n_chordwise=1, chords=None,
                 initial_sail_twist_deg=None, name=None, LLT_twist=None
                 ):

        self.__n_spanwise = n_spanwise  # number of panels (span-wise) - above the water
        self.__n_chordwise = n_chordwise  # number of panels (chord-wise) - in LLT there is line instead of panels
        self.name = name

        self.csys_transformations = csys_transformations

        self.head_mounting = head_mounting
        self.tack_mounting = tack_mounting
        """
            The geometry is described using the following CSYS.
            le - leading edge (luff) of the sail
            te - trailing edge (leech) of the sail
            below is example for the main sail.
            same for jib.

                        Z ^ (mast)     
                          |
                         /|
                        / |
                       /  |              ^ Y     
                   lem_NW +--+tem_NE    / 
                     /    |   \        /
                    /     |    \      /
                   /      |     \    /
                  /       |      \  /
                 /        |       \/
                /         |       /\
               /          |      /  \
              /           |     /    \
             /            |    /      \
            /      lem_SW |---/--------+tem_SE
           /              |  /          
  (bow) ------------------|-/-------------------------| (stern)
         \                |/                          |
    ------\---------------*---------------------------|-------------------------> X (water level)

        """

        le_NW = head_mounting
        le_SW = tack_mounting

        # mirror z coord in water surface
        # remember that direction of the lifting-line matters
        le_NW_underwater = np.array(
            [tack_mounting[0], tack_mounting[1], -tack_mounting[2]])  # leading edge South - West coordinate - mirror
        le_SW_underwater = np.array(
            [head_mounting[0], head_mounting[1], -head_mounting[2]])  # leading edge North - West coordinate - mirror

        le_NW = self.csys_transformations.rotate_point_with_mirror(le_NW)
        le_SW = self.csys_transformations.rotate_point_with_mirror(le_SW)
        le_SW_underwater = self.csys_transformations.rotate_point_with_mirror(le_SW_underwater)
        le_NW_underwater = self.csys_transformations.rotate_point_with_mirror(le_NW_underwater)

        if chords is not None:
            chords_vec = np.array([chords, np.zeros(len(chords)), np.zeros(len(chords))])
            chords_vec = chords_vec.transpose()

            rchords_vec = np.array([self.csys_transformations.rotate_point_with_mirror(c) for c in chords_vec])
            frchords_vec = np.flip(rchords_vec, axis=0)

            if initial_sail_twist_deg is not None and LLT_twist is not None:
                print(f"Applying initial_sail_twist_deg to {self.name} -  Lifting Line, mode: {LLT_twist}")
                twist_dict = {
                    'sheeting_angle_const': np.full(len(initial_sail_twist_deg), np.average(initial_sail_twist_deg)),
                    'average_const': np.full(len(initial_sail_twist_deg), np.average(initial_sail_twist_deg)),
                    'real_twist': initial_sail_twist_deg
                }
                sail_twist_deg = twist_dict[LLT_twist]
                axis = le_NW - le_SW  # head - tack
                rchords_vec = self.rotate_chord_around_le(axis, rchords_vec, sail_twist_deg)
                underwater_axis = le_NW_underwater - le_SW_underwater  # head - tack
                frchords_vec = self.rotate_chord_around_le(underwater_axis, frchords_vec,
                                                           np.flip(sail_twist_deg, axis=0))
                pass

            panels, mesh = make_panels_from_le_points_and_chords(
                [le_SW, le_NW],
                [self.__n_chordwise, self.__n_spanwise],
                rchords_vec, gamma_orientation=-1)

            panels_mirror, mesh_mirror = make_panels_from_le_points_and_chords(
                [le_SW_underwater, le_NW_underwater],
                [self.__n_chordwise, self.__n_spanwise],
                frchords_vec, gamma_orientation=-1)
        else:
            # make a lifting line instead of panels
            te_NE = le_NW  # trailing edge North - East coordinate
            te_SE = le_SW  # trailing edge South - East coordinate

            panels, mesh = make_panels_from_le_te_points(
                [le_SW, te_SE, le_NW, te_NE],
                [self.__n_chordwise, self.__n_spanwise], gamma_orientation=-1)

            te_NE_underwater = le_NW_underwater  # trailing edge North - East coordinate
            te_SE_underwater = le_SW_underwater  # trailing edge South - East coordinate

            panels_mirror, mesh_mirror = make_panels_from_le_te_points(
                [le_SW_underwater, te_SE_underwater, le_NW_underwater, te_NE_underwater],
                [self.__n_chordwise, self.__n_spanwise], gamma_orientation=-1)

        # https://stackoverflow.com/questions/33356442/when-should-i-use-hstack-vstack-vs-append-vs-concatenate-vs-column-stack
        # self.__panels = np.vstack((panels, panels_mirror))
        self.__panels = np.hstack((panels_mirror, panels))  # original version
        self.__panels1D = self.__panels.flatten()
        self.__spans = np.array([panel.get_panel_span_at_cp() for panel in self.panels1d])

    def rotate_chord_around_le(self, axis, chords_vec, sail_twist_deg_vec):
        # sail_twist = np.deg2rad(45.)
        # todo: dont forget to reverse rotations in postprocessing (plots)

        # m = rotation_matrix(axis, np.deg2rad(sail_twist_deg))
        rchords_vec = np.array([
            np.dot(rotation_matrix(axis, np.deg2rad(t)), c) for t, c in zip(sail_twist_deg_vec, chords_vec)])

        return rchords_vec

    def get_cp_points_upright(self):
        cp_points = self.get_cp_points1d()
        cp_straight_yacht = np.array([self.csys_transformations.reverse_rotations_with_mirror(p) for p in cp_points])
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


