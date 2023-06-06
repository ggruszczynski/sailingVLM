from YachtGeometry.SailGeometry import SailBaseGeometry
import numpy as np
import pandas as pd

from abc import abstractmethod, ABC
from Rotations.geometry_calc import rotation_matrix
from typing import List
from Solver import Panel
from YachtGeometry.SailGeometry import SailGeometry
from YachtGeometry.SailBaseGeometry import SailBaseGeometry

class SailSet(SailBaseGeometry):
    def __init__(self, sails: List[SailGeometry]):
        self.sails = sails
        # https://stackoverflow.com/questions/33356442/when-should-i-use-hstack-vstack-vs-append-vs-concatenate-vs-column-stack
        # self.__panels = np.vstack([sail.panels for sail in self.sails])
        self.__panels = np.hstack([sail.panels for sail in self.sails]) # original version
        self.__panels1D = self.__panels.flatten()
        self.__spans = np.array([panel.get_panel_span_at_cp() for panel in self.panels1d])

    @property
    def panels1d(self):
        return self.__panels1D

    @property
    def panels(self) -> np.array([Panel]):
        return self.__panels

    @property
    def spans(self):
        return self.__spans

    def get_cp_points_upright(self):
        cp_points_straight = np.array([[], [], []]).transpose()
        for sail in self.sails:
            cp_points_straight = np.append(cp_points_straight, sail.get_cp_points_upright(), axis=0)
        return cp_points_straight

    def sail_cp_to_girths(self):
        y_as_girths = np.array([])
        for sail in self.sails:
            y_as_girths = np.append(y_as_girths, sail.sail_cp_to_girths())
        return y_as_girths

    def extract_data_above_water_by_id(self, data, sail_no): # TODO: this is buggy - doesnot work for mainsail only
        all_cp_points = self.get_cp_points1d()[:, 2]
        reference_cp_points = self.sails[sail_no].get_cp_points1d()[:, 2]
        above_water_ref_cp_points = reference_cp_points[reference_cp_points > 0].transpose()
        index_array = np.array([np.where(all_cp_points == point)[0][0] for point in above_water_ref_cp_points])

        if isinstance(data, pd.DataFrame):
            above_water_quantities = data.iloc[index_array]
        elif isinstance(data, np.ndarray):
            if len(data.shape) == 1:
                above_water_quantities = data[index_array]
            elif len(data.shape) == 2:
                above_water_quantities = data[index_array, :]
            else:
                raise NotImplementedError
        else:
            raise NotImplementedError

        return above_water_quantities

    def extract_data_above_water_to_df(self, data):
        above_water_dfs = [pd.DataFrame(self.extract_data_above_water_by_id(data, i)) for i in range(len(self.sails))]
        merged_df_above_water = pd.concat(above_water_dfs)
        return merged_df_above_water

    def get_sail_name_for_each_element(self):
        names = []
        for i in range(len(self.sails)):
            for j in range(len(self.sails[i].panels1d)):
                names.append(self.sails[i].name)

        return pd.DataFrame({'sail_name': names})

    def get_panel_points_4debug(self):
        all_points = []
        for panel in self.panels1d:
            points = panel.get_points()
            all_points.append(points)
        return np.array(all_points)

