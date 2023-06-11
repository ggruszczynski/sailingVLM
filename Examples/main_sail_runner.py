import timeit
import shutil

from YachtGeometry.SailFactory import SailFactory
from YachtGeometry.SailSet import SailSet
from Rotations.CSYS_transformations import CSYS_transformations
from Solver.Interpolator import Interpolator
from YachtGeometry.HullGeometry import HullGeometry
from Inlet.InletConditions import InletConditions
from Inlet.Winds import ExpWindProfile, FlatWindProfile

from ResultsContainers.save_results_utils import save_results_to_file
from Solver.PanelsPlotter import display_panels_xyz_and_winds
from Solver.vlm_solver import is_no_flux_BC_satisfied

from Solver.vlm_solver import calc_circulation
from ResultsContainers.InviscidFlowResults import prepare_inviscid_flow_results_vlm
from Solver.vlm_solver import calculate_app_fs
from Utils.git_utils import get_git_branch, get_git_revision_hash

from ResultsContainers.InviscidFlowResults import InviscidFlowResults
from Solver.forces import calc_forces_on_panels_VLM_xyz
from Solver.coeff_formulas import get_CL_CD_free_wing

from InputData.rectangular_sail_vlm_case import *
# from InputData.swept_sail import *

# np.set_printoptions(precision=3, suppress=True)

start = timeit.default_timer()
interpolator = Interpolator(interpolation_type)

csys_transformations = CSYS_transformations(
    heel_deg, leeway_deg,
    v_from_original_xyz_2_reference_csys_xyz=reference_level_for_moments)

sail_factory = SailFactory(csys_transformations=csys_transformations, n_spanwise=n_spanwise, n_chordwise=n_chordwise,
                           rake_deg=rake_deg, sheer_above_waterline=sheer_above_waterline)


main_sail_geometry = sail_factory.make_main_sail(
    main_sail_luff=main_sail_luff,
    boom_above_sheer=boom_above_sheer,
    main_sail_chords=interpolator.interpolate_girths(main_sail_girths, main_sail_chords, n_spanwise + 1),
    sail_twist_deg=interpolator.interpolate_girths(main_sail_girths, main_sail_centerline_twist_deg, n_spanwise + 1),
    LLT_twist=LLT_twist,
    interpolated_camber=interpolator.interpolate_girths(main_sail_girths, main_sail_camber, n_spanwise + 1),
    interpolated_distance_from_LE=interpolator.interpolate_girths(main_sail_girths, main_sail_max_camber_distance_from_luff, n_spanwise + 1))

sail_set = SailSet([main_sail_geometry])

# jib_geometry = sail_factory.make_jib(
#     jib_luff=jib_luff,
#     foretriangle_base=foretriangle_base,
#     foretriangle_height=foretriangle_height,
#     jib_chords=interpolator.interpolate_girths(jib_girths, jib_chords, n_spanwise + 1),
#     sail_twist_deg=interpolator.interpolate_girths(jib_girths, jib_centerline_twist_deg, n_spanwise + 1),
#     mast_LOA=mast_LOA,
#     LLT_twist=LLT_twist,
#     interpolated_camber=interpolator.interpolate_girths(jib_girths, jib_sail_camber, n_spanwise + 1),
#     interpolated_distance_from_LE=interpolator.interpolate_girths(jib_girths, jib_sail_max_camber_distance_from_luff, n_spanwise + 1)
# )
# sail_set = SailSet([jib_geometry, main_sail_geometry])


wind = FlatWindProfile(alpha_true_wind_deg, tws_ref, SOG_yacht)
# wind = ExpWindProfile(
#     alpha_true_wind_deg, tws_ref, SOG_yacht,
#     exp_coeff=wind_exp_coeff,
#     reference_measurment_height=wind_reference_measurment_height,
#     reference_water_level_for_wind_profile=reference_water_level_for_wind_profile)

inlet_condition = InletConditions(wind, rho=rho, panels1D=sail_set.panels1d)

# hull = HullGeometry(sheer_above_waterline, foretriangle_base, csys_transformations, center_of_lateral_resistance_upright)

gamma_magnitude, v_ind_coeff = calc_circulation(inlet_condition.V_app_infs, sail_set.panels)
V_induced_at_ctrl_p, V_app_fs_at_ctrl_p = calculate_app_fs(inlet_condition, v_ind_coeff, gamma_magnitude)

assert is_no_flux_BC_satisfied(V_app_fs_at_ctrl_p, sail_set.panels)

inviscid_flow_results = prepare_inviscid_flow_results_vlm(gamma_magnitude, sail_set, inlet_condition, csys_transformations)
# inviscid_flow_results.estimate_heeling_moment_from_keel(hull.center_of_lateral_resistance)

df_components, df_integrals, df_inlet_IC = save_results_to_file(inviscid_flow_results, None, inlet_condition, sail_set, output_dir_name)
shutil.copy(os.path.join(case_dir, case_name), os.path.join(output_dir_name, case_name))

print(f"-------------------------------------------------------------")
print(f"Notice:\n"
      f"\tThe forces [N] and moments [Nm] are without profile drag.\n"
      f"\tThe the _COG_ CSYS is aligned in the direction of the yacht movement (course over ground).\n"
      f"\tThe the _COW_ CSYS is aligned along the centerline of the yacht (course over water).\n"
      f"\tNumber of panels (sail set with mirror): {sail_set.panels.shape}")

print(df_integrals)

q = 0.5 * rho * (np.linalg.norm(V_inf) ** 2)
S = 2*chord*half_wing_span
AR = 2 * half_wing_span / chord

sail_mirror_multiplier = 2
F_xyz_total = sail_mirror_multiplier*inviscid_flow_results.F_xyz_total
# from Rotations.geometry_calc import rotation_matrix
# A = rotation_matrix([0, 0, 1], np.deg2rad(-AoA_deg))      # if you rotate the wind instead of the lifting surface
# F = np.dot(A, F_xyz_total)    # by definiton the Lift force is the component of total aerodynamic force being perpendicular to V_inf
F = F_xyz_total
CX_vlm = F[0] / (q*S)
CY_vlm = F[1] / (q*S)
CZ_vlm = F[2] / (q*S)

a_VLM = CY_vlm / np.deg2rad(AoA_deg)
CL_analytical, CD_ind_analytical, a_analytical = get_CL_CD_free_wing(AR, AoA_deg, sweep_half_chord_deg=sweep_angle_deg)

print(f"\nAspect Ratio {AR}")
print(f"CL_vlm  {CY_vlm:.6f}    \t CL_VLM_KTH {4.897*np.deg2rad(AoA_deg):.6f}   \t CL_analytical      {abs(CL_analytical):.6f} ")
print(f"CD_vlm  {CX_vlm:.6f}    \t CD_vlm_KTH {0.0235:.6f}  \t CD_ind_analytical  {CD_ind_analytical:.6f}")
print(f"a_VLM   {abs(a_VLM):.6f}     \t a_VLM_KTH  {4.897:.6f}   \t a_analytical       {a_analytical:.6f}")

# rows_to_display = ['M_total_heeling', 'M_total_sway', 'F_sails_drag'] # select rows to print
# print(df_integrals[df_integrals['Quantity'].isin(rows_to_display)])

print(f"\nCPU time: {float(timeit.default_timer() - start):.2f} [s]")

print("Preparing visualization.")
display_panels_xyz_and_winds(sail_set.panels1d, inlet_condition, inviscid_flow_results, hull=None)

print(f"Code version\t branch: {get_git_branch()} \t commit hash: {get_git_revision_hash()}")
print("Done.")
