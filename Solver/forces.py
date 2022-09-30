import numpy as np
from Solver.vlm_solver import calc_induced_velocity

from Solver.Panel import Panel
from typing import List


def calc_moment_arm_in_shifted_csys(cp_points, v_from_old_2_new_csys):
    dx, dy, dz = v_from_old_2_new_csys
    N = len(cp_points)
    r = np.copy(cp_points)
    for i in range(N):
        r[i, 0] -= dx  # yawing moment arm
        r[i, 1] -= dy  # pitching moment arm
        if cp_points[i, 2] > 0:
            r[i, 2] -= dz  # heeling moment arm
        else:
            r[i, 2] += dz

    return r


def calc_moments(cp_points, forces):
    moments = np.cross(cp_points, forces)
    # Mx_test1 = cp_points[:, 2] * force_xyz[:, 1]
    # Mx_test2 = cp_points[:, 1] * force_xyz[:, 2]
    # Mx_test = Mx_test1 + Mx_test2
    return moments


def extract_above_water_quantities(quantities, cp_points):
    # potential pitfalls
    # 1) above_water_mask = Mx < 0  # don't do this - consider a negative camber design ;p
    # 2) total_force = np.sum(force, axis=0) / 2.  # don't do this - the z component of the force would vanish due to mirror effect
    above_water_mask = cp_points[:, 2] > 0
    above_water_quantities = quantities[above_water_mask]
    total_above_water_quantities = np.sum(above_water_quantities, axis=0)  # heeling, sway, yaw (z-axis)
    return above_water_quantities, total_above_water_quantities


def calc_V_at_cp(V_app_infw, gamma_magnitude, panels):
    """
    :param V_app_infw: apparent wind
    :param gamma_magnitude: vector with circulations
    :param panels:
    :return: Wind at cp = apparent wind + wind_induced_at_cp
    """
    panels_1d = panels.flatten()
    N = len(panels_1d)
    v_ind_coeff = np.full((N, N, 3), 0., dtype=float)

    for i in range(0, N):
        if i % 25 == 0:
            print(f"assembling v_ind_coeff matrix at cp {i}/{N}")

        cp = panels_1d[i].get_cp_position()
        for j in range(0, N):
            # velocity induced at i-th control point by j-th vortex
            v_ind_coeff[i][j] = panels_1d[j].get_induced_velocity(cp, V_app_infw[j])

    V_induced_at_cp = calc_induced_velocity(v_ind_coeff, gamma_magnitude)
    V_app_fs_at_cp = V_app_infw + V_induced_at_cp
    return V_app_fs_at_cp, V_induced_at_cp

def calc_force_LLT_xyz(V_app_fs_at_cp, gamma_magnitude, panels1d, rho):
    span_vectors = np.array([p.get_span_vector() for p in panels1d])
    for i in range(0, len(gamma_magnitude)):
        gamma = span_vectors[i] * gamma_magnitude[i]
        panels1d[i].force_xyz = rho * np.cross(V_app_fs_at_cp[i], gamma)

def calc_force_VLM_xyz(V_app_infw, gamma_magnitude, panels: List[Panel], rho):
    """
    Katz and Plotkin, p. 346 Chapter 12 / Three-Dimensional Numerical Solution
    f. Secondary Computations: Pressures, Loads, Velocities, Etc
    Eq (12.25)

    force = rho* (V_app_fw_at_cp x gamma*span)
    :param V_app_infw: apparent wind
    :param gamma_magnitude: vector with circulations
    :param panels: list
    :param rho: air density
    :return: force
    """

    V_app_fs_at_cp, V_induced_at_cp = calc_V_at_cp(V_app_infw, gamma_magnitude, panels)

    V_app_fs_at_cp_re = V_app_fs_at_cp.reshape(panels.shape[0], panels.shape[1], 3)
    gamma_re = gamma_magnitude.reshape(panels.shape)
    force_re_xyz = np.full((panels.shape[0], panels.shape[1], 3), 0., dtype=float)
    for i in range(0, panels.shape[0]):
        for j in range(0, panels.shape[1]):
            if i == 0:  # leading edge only
                gamma = panels[i, j].get_span_vector() * gamma_re[i, j]
            else:
                gamma = panels[i, j].get_span_vector() * (gamma_re[i, j]-gamma_re[i-1, j])

            force_tmp = rho * np.cross(V_app_fs_at_cp_re[i, j], gamma)
            force_re_xyz[i, j, :] = force_tmp
            panels[i, j].force_xyz = force_tmp

    return force_re_xyz, V_app_fs_at_cp, V_induced_at_cp


def get_p_from_panels(panels):
    p = np.zeros(shape=panels.shape)
    for i in range(panels.shape[0]):
        for j in range(panels.shape[1]):
            p[i, j] = panels[i, j].pressure
    return p


def get_forces_from_panels(panels):
    forces = np.full((panels.shape[0], panels.shape[1], 3), 0., dtype=float)
    for i in range(panels.shape[0]):
        for j in range(panels.shape[1]):
            forces[i, j, :] = panels[i, j].force_xyz
    return forces

def determine_vector_from_its_dot_and_cross_product(F, r_dot_F, r_cross_F):
    # https://math.stackexchange.com/questions/246594/what-is-vector-division
    # https://math.stackexchange.com/questions/1683996/determining-an-unknown-vector-from-its-cross-and-dot-product-with-known-vector
    # we know that: a x ( b x c) = b(a . c) - c(a.b)
    # so: F x ( r x F) = r(F . F) - F(F.r)
    # one shall get r, assuming that both the cross and dot product are known:
    # r = (F x (r x F) + F(F . r))/(F . F)

    F_dot_F = np.dot(F, F)
    R = (np.cross(F, r_cross_F) + F * r_dot_F) / F_dot_F
    return np.array(R)


