import numpy as np
from numba import config, njit, threading_layer, jit,  prange

# set the threading layer before any parallel target compilation
config.THREADING_LAYER = 'threadsafe'
# @jit(parallel=True, forceobj=True)


@njit(parallel=True)
def foo(a, b):
    x = np.log(2*a + b)
    x2 = np.exp(np.tanh(x))
    return x + x2



# @jit(nopython=False, forceobj=True, parallel=True)
config.THREADING_LAYER = 'omp'
# @jit(parallel=True)
def assembly_sys_of_eq_helper(V_app_infw, panels1D, ctr_p, panel_surf_normal, N):
    # v_ind_coeff_tmp = np.full((N, 3), 0., dtype=float)
    v_ind_coeff_tmp = np.zeros((N, 3), dtype=np.double)
    A_tmp = np.zeros(N, dtype=np.double)

    for j in prange(0, N):  # TODO: notice prange instead of range
        # velocity induced at i-th control point by j-th vortex
        v_ind_coeff_tmp[j] = panels1D[j].get_horse_shoe_induced_velocity(ctr_p, V_app_infw[j])
        A_tmp[j] = np.dot(v_ind_coeff_tmp[j], panel_surf_normal)

    # demonstrate the threading layer chosen
    # print("Threading layer chosen: %s" % threading_layer())
    return v_ind_coeff_tmp, A_tmp


# @jit(nopython=False, forceobj=True)
def assembly_sys_of_eq(V_app_infw, panels1D, N):


    A = np.zeros((N, N), dtype=np.double)  # Aerodynamic Influence Coefficient matrix
    RHS = np.zeros(N, dtype=np.double)
    v_ind_coeff = np.zeros((N, N, 3), dtype=np.double)

    for i in range(0, N):
        panel_surf_normal = panels1D[i].get_normal_to_panel()
        ctr_p = panels1D[i].get_ctr_point_position()
        RHS[i] = -np.dot(V_app_infw[i], panel_surf_normal)

        # for j in range(0, N):
        #         # velocity induced at i-th control point by j-th vortex
        #         v_ind_coeff[i][j] = panels1D[j].get_horse_shoe_induced_velocity(ctr_p, V_app_infw[j])
        #         A[i][j] = np.dot(v_ind_coeff[i][j], panel_surf_normal)
        v_ind_coeff_tmp, A_tmp = assembly_sys_of_eq_helper(V_app_infw, panels1D, ctr_p, panel_surf_normal, N)
        v_ind_coeff[i][:] = v_ind_coeff_tmp
        A[i][:] = A_tmp


    # x = np.arange(5000000.)
    # x = np.arange(500000000.)
    # y = x.copy()
    # foo(x, y)

    return A, RHS, v_ind_coeff  # np.array(v_ind_coeff)

# @jit(nopython=False, forceobj=True)
def calc_circulation(V_app_ifnw, panels):
    # it is assumed that the freestream velocity is V [vx,0,vz], where vx > 0

    panels1D = panels.flatten()
    A, RHS, v_ind_coeff = assembly_sys_of_eq(V_app_ifnw, panels1D, len(panels1D))
    gamma_magnitude = np.linalg.solve(A, RHS)

    return gamma_magnitude, v_ind_coeff


def calc_induced_velocity(v_ind_coeff, gamma_magnitude):
    N = len(gamma_magnitude)
    V_induced = np.full((N, 3), 0., dtype=float)
    for i in range(N):
        for j in range(N):
            V_induced[i] += v_ind_coeff[i][j] * gamma_magnitude[j]

    return V_induced


def is_no_flux_BC_satisfied(V_app_fw, panels):
    panels1D = panels.flatten()
    N = len(panels1D)
    flux_through_panel = np.zeros(shape=N)
    panels_area = np.zeros(shape=N)

    for i in range(0, N):
        panel_surf_normal = panels1D[i].get_normal_to_panel()
        panels_area[i] = panels1D[i].get_panel_area()
        flux_through_panel[i] = -np.dot(V_app_fw[i], panel_surf_normal)

    for area in panels_area:
        if np.isnan(area) or area < 1E-14:
            raise ValueError("Solution error, panel_area is suspicious")

    for flux in flux_through_panel:
        if abs(flux) > 1E-12 or np.isnan(flux):
            raise ValueError("Solution error, there shall be no flow through panel!")

    return True


def calculate_app_fs(inletConditions, v_ind_coeff, gamma_magnitude):
    V_induced = calc_induced_velocity(v_ind_coeff, gamma_magnitude)
    V_app_fs = inletConditions.V_app_infs + V_induced
    return V_induced, V_app_fs
