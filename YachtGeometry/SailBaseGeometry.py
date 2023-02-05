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

    def get_cp_points(self):
        return get_stuff_from_panels(self.panels, 'cp_position', (self.panels.shape[0], self.panels.shape[1], 1))

    def get_cp_points1d(self):
        return get_stuff_from_panels(self.panels1d, 'cp_position', (self.panels1d.shape[0], 3))


    @property
    def pressures(self):
        return get_stuff_from_panels(self.panels, 'pressure', (self.panels.shape[0], self.panels.shape[1], 1))

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

