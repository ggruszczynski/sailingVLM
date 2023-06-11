import numpy as np
import pandas as pd

from abc import abstractmethod, ABC
from Rotations.geometry_calc import rotation_matrix
from Solver import Panel
from Rotations.CSYS_transformations import CSYS_transformations
from Solver.mesher import make_panels_from_le_te_points, make_panels_from_le_points_and_chords
from typing import List

from Solver.forces import get_stuff_from_panels

# np.set_printoptions(precision=3, suppress=True)

class SailBaseGeometry:
    def __init__(self,
                 head_mounting: np.array, tack_mounting: np.array,
                 csys_transformations: CSYS_transformations,
                 n_spanwise=10, n_chordwise=1, name=None):
        self._n_spanwise = n_spanwise  # number of panels (span-wise) - above the water
        self._n_chordwise = n_chordwise  # number of panels (chord-wise) - in LLT there is line instead of panels
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

        self.le_NW = head_mounting
        self.le_SW = tack_mounting

        # mirror z coord in water surface
        # remember that direction of the lifting-line matters
        self.le_NW_underwater = np.array(
            [tack_mounting[0], tack_mounting[1], -tack_mounting[2]])  # leading edge South - West coordinate - mirror
        self.le_SW_underwater = np.array(
            [head_mounting[0], head_mounting[1], -head_mounting[2]])  # leading edge North - West coordinate - mirror

    @property
    @abstractmethod
    def spans(self):
        pass

    @property
    @abstractmethod
    def panels1d(self):
        pass

    @property
    @abstractmethod
    def panels(self):
        pass

    @abstractmethod
    def extract_data_above_water_to_df(self, data):
        pass

    @abstractmethod
    def sail_cp_to_girths(self):
        pass

    @abstractmethod
    def get_cp_points_upright(self):
        pass

    def get_ctr_points(self):
        return get_stuff_from_panels(self.panels, 'ctr_point_position', (self.panels.shape[0], self.panels.shape[1], 3))

    def get_ctr_points1d(self):
        return get_stuff_from_panels(self.panels1d, 'ctr_point_position', (self.panels1d.shape[0], 3))
    def get_cp_points(self):
        return get_stuff_from_panels(self.panels, 'cp_position', (self.panels.shape[0], self.panels.shape[1], 3))

    def get_cp_points1d(self):
        return get_stuff_from_panels(self.panels1d, 'cp_position', (self.panels1d.shape[0], 3))


    @property
    def pressures(self):
        return get_stuff_from_panels(self.panels, 'pressure', (self.panels.shape[0], self.panels.shape[1], 1))

    @property
    def coeffs_of_pressure(self):
        return get_stuff_from_panels(self.panels, 'coeff_of_pressure', (self.panels.shape[0], self.panels.shape[1], 1))

    @property
    def forces_xyz(self):
        return get_stuff_from_panels(self.panels, 'force_xyz', (self.panels.shape[0], self.panels.shape[1], 3))

    @property
    def V_app_fs_at_cp(self):
        return get_stuff_from_panels(self.panels, 'V_app_fs_at_cp', (self.panels.shape[0], self.panels.shape[1], 3))
        # return get_V_app_fs_at_cp_from_panels(self.panels)

    @property
    def V_induced_at_cp(self):
        return get_stuff_from_panels(self.panels, 'V_induced_at_cp', (self.panels.shape[0], self.panels.shape[1], 3))
        # return get_V_induced_at_cp_from_panels(self.panels)

