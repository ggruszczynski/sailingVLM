import os
import numpy as np
import time

# OUTPUT DIR
case_name = os.path.basename(__file__)  # get name of the current file
case_dir = os.path.dirname(os.path.realpath(__file__))  # get dir of the current file
time_stamp = time.strftime("%Y-%m-%d_%Hh%Mm%Ss")
output_dir_name = os.path.join("results_RC44_GENOA_VLM", time_stamp)

# SOLVER SETTINGS
n_spanwise = 25  # No of control points (above the water) per sail, recommended: 20-50
n_chordwise = 15  # No of control points (above the water) per sail, recommended: 20-50
# AoA_0lift_iterations = 15   # recommended range [10-30] note that running 1 iteration produce a shape without AoA_0lift constraint
# AoA_0lift_max_change_per_iteration_deg = 0.05  # [deg] recommended range [0.1 - 0.01]
# wind_sub_iterations = 3  # recommended value 3 - enough to converge CL and Heeling moment constrains
interpolation_type = "spline"  # either "spline" or "linear"
LLT_twist = "real_twist"  # defines how the Lifting Line discretize the sail twist.
# It can be "sheeting_angle_const" or "average_const" or "real_twist"

# SAILING CONDITIONS
leeway_deg = 0.0   # [deg]
heel_deg = 13.0    # [deg]
SOG_yacht = 3.4460  # [m/s] yacht speed - speed over ground (leeway is a separate variable)
tws_ref = 3.0860    # [m/s] true wind speed
alpha_true_wind_deg = 45.0  # [deg] true wind angle (with reference to course over ground) => Course Wind Angle to the boat track = true wind angle to centerline + Leeway
reference_water_level_for_wind_profile = -1.23  # [m] this is an attempt to mimick the deck effect
# by lowering the sheer_above_waterline
# while keeping the wind profile as in original geometry
# this shall be negative (H = sail_ctrl_point - water_level)
wind_exp_coeff = 0.1000  # [-] coefficient to determine the exponential wind profile
wind_reference_measurment_height = 22.0  # [m] reference height for exponential wind profile
rho = 1.184  # air density [kg/m3]

# GEOMETRY OF THE RIG
main_sail_luff = 17.52  # [m]
jib_luff = 17.40  # [m]
foretriangle_height = 16.98  # [m]
foretriangle_base = 5.13  # [m]
sheer_above_waterline = 1.23  # [m]
boom_above_sheer = 1.41  # [m]
rake_deg = 95.  # rake angle [deg]
mast_LOA = 0.24  # [m]

# INPUT - GEOMETRY OF THE SAIL
main_sail_girths = np.array([0.00, 1./4, 1./2, 3./4, 1.00])
main_sail_chords = np.array([5.37, 4.98, 4.45, 3.72, 2.49])
main_sail_centerline_twist_deg = 17. * main_sail_girths + 0
main_sail_camber = 1*np.array([0.028, 0.083, 0.095, 0.087, 0.013])
main_sail_max_camber_distance_from_luff = np.array([0.5, 0.5, 0.5, 0.5, 0.5])  # distance from luff (leading edge)

jib_girths = np.array([0.00, 1./4, 1./2, 3./4, 1.00])
jib_chords = np.array([7.32, 5.39, 3.55, 1.76, 0.09])
jib_centerline_twist_deg = 24. * jib_girths + 8.6
jib_sail_camber = 1*np.array([0.055, 0.100, 0.130, 0.131, 0.01])
jib_sail_max_camber_distance_from_luff = np.array([0.5, 0.5, 0.5, 0.5, 0.5])  # distance from luff (leading edge)

# OPTIMIZATION CONSTRAINTS (INEQUALITY)
# main_sail_CLmax = 1.1 * np.array([1.0, 1.0, 1.0, 1.0, 1.0])
# main_sail_CLmin = 1.0 * np.array([0.33, 0.43, 0.0, 0.0, 0.0])
# main_sail_AoA_0lift_deg_min = -11.0 * np.array([1.0, 1.0, 1.0, 1.0, 1.0])
#
# jib_CLmax = 1.0 * np.array([0.5, 1.0, 0.9, 0.54, 0.01])  # kill the Jib's circulation to avoid bump on the main sail
# jib_CLmin = 0.0 * np.array([1.0, 1.0, 1.0, 1.0, 1.0])
# jib_AoA_0lift_deg_min = -10.0 * np.array([1.0, 1.0,  1.0,  1.0,  1.0])
#
# # OPTIMIZATION CONSTRAINTS (EQUALITY)
# imposed_heeling_moment = -25200  # for the sails [Nm]

# REFERENCE CSYS
# The origin of the default CSYS is located @ waterline level and aft face of the mast
# The positive x-coord: towards stern
# The positive y-coord: towards leeward side
# The positive z-coord: above the water
# To shift the default CSYS, adjust the 'reference_level_for_moments' variable.
# Shifted CSYS = original + reference_level_for_moments
# As a results the moments will be calculated around the new origin.

# yaw_reference [m] - distance from the aft of the mast towards stern, at which the yawing moment is calculated.
# sway_reference [m] - distance from the aft of the mast towards leeward side. 0 for symmetric yachts ;)
# heeling_reference [m] - distance from the water level,  at which the heeling moment is calculated.
reference_level_for_moments = np.array([-5.099, 0, 0])  # [yaw_reference, sway_reference, heeling_reference]
# reference_level_for_moments = np.array([0, 0, 0])  # [yaw_reference, sway_reference, heeling_reference]

# GEOMETRY OF THE KEEL
# to estimate heeling moment from keel, does not influence the optimizer.
# reminder: the z coord shall be negative (under the water)
center_of_lateral_resistance_upright = np.array([0, 0, -1.25])  # [m] the coordinates for a yacht standing in upright position