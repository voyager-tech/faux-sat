# Utilized Modules
import numpy as np
from orbital_analyses import u
from copy import deepcopy as dc
from datetime import datetime, timedelta
import Requirements as Req
import constants as const

from orbital_analyses.Transform_State import gregorian2julian_date
from orbital_analyses.Transform_Coordinate import IAU_2000Reduction

# TODO: Remove this section

# Determine Inital State
# Initial Rad / Vel - R0, V0 and epoch
# TODO: Add in handling if rad and vel given as 3x1 (or handle on input script)
# r0 = (Req.Rad).magnitude  # km
r0 = np.array([1949.7385039, -5139.97811485, -4852.9027131])  # km
# v0 = (Req.Vel).magnitude  # km / sec
v0 = np.array([-6.57485096794, 0.215327040815, -3.14995010529])  # km / sec
# EpochGD_0 = Req.GD_UTC  # Gregorian Date
Epoch_GD0 = np.array([2004, 4, 6, 0, 0, 0])  # Gregorian Date
Epoch_JD0 = gregorian2julian_date(Epoch_GD0)  # Julian Date # TODO: Remove?
# Step = Req.Ssize
step_size = 1 * u.min
steps = 100
accuracy = 1e-012
max_iterations = 50

# TODO: Include time at each state into iterations and output file
# Use convert_time type code (or old propagator) and make into a function
############################################################################
# %%


def gauss_jackson_prop(r0, v0, Epoch_GD0, step_size, steps, max_iterations=50, accuracy=1e-12):
    """
    Propagate an initial orbit forward in time by a specified timespan.

    Propagator is an 8th order numerical integrator corrector that utilizes a
    f & g series with Newton-Raphson iteration to calculate the initial
    positions and velocities to propagate forward from.

    Parameters
    ----------
    r0 : numpy array (3,) - [X, Y, Z]
        - Initial radius vector components defined in km in the J2000eq frame
    v0 : numpy array (3,) - [VX, VY, VZ]
        - Initial velocity vector components defined in km in the J2000eq frame
    Epoch_GD0 : numpy array (6,) - [Yr, Mo, Day, Hr, Min, Sec]
        - Gregorian Date of the initial state
    step_size : float
        - Size of each propagator-calculated step
        - Multiplied by a time value from the package, Pint, or just in seconds
        - Note : Best practice for step_size <= 5 min (600 sec) for convergence
    steps : int
        - Total number of steps to be evauated in tandem with step_size
    max_iterations : int, optional
        - The desired number of iterations to attempt to converge a new step
    accuracy : float, optional
        - The desired accuracy for the convergence at each new step

    Returns
    -------
    location_c : tuple - (GD, rad, vel, acc, frame)
        - Returns arrays / lists of information at each step
        - rad, vel, acc arrays are in the frame specified in the frame list

    References
    ----------
    [1] Berry, Matthew M, and Liam M Healy. “Implementation of Gauss-Jackson
    Integration for Orbit Propagation.” The Journal of the Astronautical
    Sciences, vol. 52, no. 3, July 2004, pp. 331–357.
    - https://drum.lib.umd.edu/bitstream/handle/1903/2202/2004-berry-healy-jas.pdf

    [2] D. Vallado, `Fundamentals of Astrodynamics and Applications`. 4th ed.,
    Microcosm Press, 2013.
        - Alg. 8 - pg. 93, Sec. 8.6 - pg. 538
    """
    # Utilized functions for the Newton-Raphson iterator
    def guess1(Tdelt, alpha):
        # Vallado Alg. 8, Pg. 93
        X0 = (np.sqrt(const.GRAV_PARAMATER['EARTH']) * Tdelt * alpha)
        return X0

    def guess2(r0, v0, Tdelt):
        # Vallado Alg. 8, Pg. 93
        h0 = np.cross(r0, v0)
        h0_2 = np.linalg.norm(np.inner(h0, h0))
        p = (h0_2 / const.GRAV_PARAMATER['EARTH'])
        s = ((np.arccos(3 * np.sqrt(const.GRAV_PARAMATER['EARTH'] /
                                    (p ** 3)) * Tdelt)) / 2)
        w = np.arctan(np.cbrt(np.tan(s)))
        X0 = (np.sqrt(p) * 2 * np.cot(2 * w))
        return X0

    def guess3(r0, v0, Tdelt, alpha, r0_mag):
        # Vallado Alg. 8, Pg. 93
        a = (1 / alpha)
        X0 = (np.sign(Tdelt) * np.sqrt(-a) *
              np.log((-2 * const.GRAV_PARAMATER['EARTH'] * alpha * Tdelt) /
                     ((np.linalg.norm(np.inner(r0, v0))) + np.sign(Tdelt) *
                      np.sqrt(-const.GRAV_PARAMATER['EARTH'] * a) *
                      (1 - r0_mag * alpha))))
        return X0

    def c2c3(psi):
        # Vallado Alg. 1, Pg. 63
        if psi > 1e-06:
            c2 = ((1 - np.cos(np.sqrt(psi))) / psi)
            c3 = ((np.sqrt(psi) - np.sin(np.sqrt(psi))) / np.sqrt(psi ** 3))
        elif psi < -1e-06:
            c2 = ((1 - np.cosh(np.sqrt(-psi))) / psi)
            c3 = ((np.sin(np.sqrt(-psi)) - np.sqrt(-psi)) /
                  np.sqrt((-psi) ** 3))
        else:
            c2 = (1 / 2)
            c3 = (1 / 6)
        return c2, c3

    # 1. Use f & g series to calculate 8 (-4,4) rad/vel surrounding the epoch_0
    # Initialize arrays to place all calculated state vectors into
    rad = np.zeros((steps + 1, 3), dtype=float)
    vel = np.zeros((steps + 1, 3), dtype=float)
    acc = np.zeros((steps + 1, 3), dtype=float)
    GD = np.zeros((steps + 1, 6), dtype=float)
    frame = []
    # Initialize temp array for calculating convergence at each step
    rad_temp_inertial = np.zeros((9, 3), dtype=float)
    vel_temp_inertial = np.zeros((9, 3), dtype=float)
    rad_temp = np.zeros((9, 3), dtype=float)
    vel_temp = np.zeros((9, 3), dtype=float)
    acc_temp = np.zeros((9, 3), dtype=float)
    GD_temp = dc(Epoch_GD0)
    # Determine all other variables needed for loop
    r0_mag = np.linalg.norm(r0)
    v0_2 = np.inner(v0, v0)  # km^2 / sec^2
    # Set initial step value
    if type(step_size) != int and type(step_size) != float:
        h = ((step_size.to(u.sec)).magnitude)
    else:
        h = dc(step_size)

    # Newton-Raphson iteration using Universal Variabes to calculate
    # the 8 (-4, 4) rad/vel surrounding the epoch_0
    # Vallado - Alg. 8, pg. 93
    # loop to calculate rad/vel vectors of 8 other positions defined by Step
    n_range = np.linspace(-4, 4, 9, dtype=int)
    for n in n_range:
        # Compute the Rad and Vel vectors at each step around the origin
        Tdelt = n * h
        alpha0 = ((-v0_2 / const.GRAV_PARAMATER['EARTH']) + (2 / r0_mag))
        if alpha0 > 0.000001:                   # Circle or Ellipse
            Xn = guess1(Tdelt, alpha0)  # Xn = [sqrt(km)]
        elif np.absolute(alpha0) < 0.000001:    # Parabola
            Xn = guess2(r0, v0, Tdelt)  # Xn = [sqrt(km)]
        elif alpha0 < -0.000001:                # Hyperbola
            Xn = guess3(r0, v0, Tdelt, alpha0, r0_mag)  # Xn = [sqrt(km)]
        err = 1.
        while err >= 1e-06:
            Psi = ((Xn ** 2) * alpha0)
            c2, c3 = c2c3(Psi)
            r = (((Xn ** 2) * c2) + ((np.linalg.norm(np.inner(r0, v0)) /
                                     np.sqrt(const.GRAV_PARAMATER['EARTH'])) *
                                     Xn * (1 - (Psi * c3))) +
                 (r0_mag * (1 - (Psi * c2))))
            Xn1 = (Xn + ((((np.sqrt(const.GRAV_PARAMATER['EARTH']) * Tdelt) -
                           ((Xn ** 3) * c3) -
                           ((np.linalg.norm(np.inner(r0, v0)) /
                             np.sqrt(const.GRAV_PARAMATER['EARTH'])) *
                           ((Xn ** 2) * c2)) -
                           (r0_mag * Xn * (1 - (Psi * c3))))) / r))
            Xn = dc(Xn1)
            err = np.absolute(Xn - Xn1)
        # With iterated values, calculate relevant f and g functions
        f = (1 - ((Xn ** 2) / (r0_mag)) * c2)
        g = (Tdelt - ((Xn ** 3) / np.sqrt(const.GRAV_PARAMATER['EARTH'])) * c3)
        f_dot = ((np.sqrt(const.GRAV_PARAMATER['EARTH']) / (r * r0_mag)) *
                 Xn * ((Psi * c3) - 1))
        g_dot = (1 - ((Xn ** 2) / (r)) * c2)
        # Check for correctness
        # print((f * g_dot) - (f_dot * g))  # TODO: Create check here (Should = 1)
        # Calculate new rad / vel vectors
        rad_twoBod = ((f * r0) + (g * v0))
        vel_twoBod = ((f_dot * r0) + (g_dot * v0))
        # Place rad / vel vectors into matrices
        if n == 0:
            rad_temp[4, :] = r0  # TODO: Change to rad_temp_inertial
            vel_temp[4, :] = v0  # TODO: Change to vel_temp_inertial
        else:
            rad_temp[n + 4, :] = rad_twoBod  # TODO: Change to rad_temp_inertial
            vel_temp[n + 4, :] = vel_twoBod  # TODO: Change to vel_temp_inertial
        # Convert rad and vel to Fixed frame (ECEF) from their Inertial frame
        #for m in range(9):
        #    rad_temp[m, :], vel_temp[m, :] = IAU_2000Reduction(rad_temp_inert[m, :], vel_temp_inert[m, :], gd_UTC, 1)

# %% # TODO: Remove Later######################################################
    # 2. Evaluate 9 acceeration vectors for the 9 states determined above
    # Acceleration w/o any perturbations
    def evaluate_accel(rad, vel, bodies, J2=0, drag=0, three_body=0, srp=0):
        # TODO: Docstring here
        # TODO: make sure user has entered 1s or 0s for each option
        # Initialize all acceleration vectors
        acc_2body = np.zeros((len(rad), 3), dtype=float)
        acc_J2 = np.zeros(3, dtype=float)
        acc_drag = np.zeros(3, dtype=float)
        acc_3body = np.zeros(3, dtype=float)
        acc_srp = np.zeros(3, dtype=float)
        # Evaluate 2-body accelerations (Only if not evaluating 3-Body)
        if three_body == 0:
            for j in range(len(rad)):
                acc_2body[j, :] = -((rad[j, :] *
                                     const.GRAV_PARAMATER['EARTH']) /
                                    (np.linalg.norm(rad[j, :]) ** 3))
        # If selected, evaluate J2 accelerations (Planetary oblateness (grav.))
        if J2 == 1:
            0
        # If selected, evaluate drag accelerations (Atmosphereic)
        if drag == 1:
            0
        # If selected, evaluate 3-body accelerations (Includes 2-body)
        if three_body == 1:
            # D. Vallado - pg. 574
            # Collect Mu values of all requested bodies
            bodies_Mu = np.zeros(len(bodies), dtype=float)
            for j in range(len(bodies)):
                bodies_Mu[j] = const.GRAV_PARAMATER[bodies[j]]
            # D. Vallado - pg. 295-298
            # Calculate distance vectors from the satellite to the 3rd bodies
            

        # If selected, evaluate srp accelerations (Solar Radiation Pressure)
        if srp == 1:
            0
        # Sum all calculated accelerations

        return

    for j in range(9):
        acc_temp[j, :] = -((rad_temp[j, :] * const.GRAV_PARAMATER['EARTH']) /
                           (np.linalg.norm(rad_temp[j, :]) ** 3))
    # Incorporate Special Perturbations into acceleration model
    # Vallado Alg. 64, Pg. 591
    # TODO: Gather more accurate models of perturbative forces
    # https://sourceforge.net/p/gmat/git/ci/GMAT-R2018a/tree/application/data/
    # https://sourceforge.net/p/gmat/git/ci/GMAT-R2018a/tree/prototype/StateConv/
    # http://www.hayabusa.isas.jaxa.jp/kawalab/dromobile/Papers/HernandoAyuso%20Dromobile%202016.pdf
    # https://www.sciencedirect.com/science/article/pii/0898122186900258

# %% # TODO: Remove Later######################################################

    # 3. Converging the accelerations
    # Import all coefficient arrays
    # Eigth order summed adams coefficients in Ordinate Form: a(j,k), b(j,k)
    a = np.load(r'orbital_analyses\coefficient_matrices\gauss_jackson_prop\GJ_a_coeff.npy')
    b = np.load(r'orbital_analyses\coefficient_matrices\gauss_jackson_prop\GJ_b_coeff.npy')
    # Initialize all loop arrays and variables before loop
    acc_temp1 = np.zeros((9, 3), dtype=float)
    diff_acc = 1
    # 3a. Calculate s_0 and S_0
    # Calculate s_0 from C1
    C1_sum = 0
    for m in range(9):
        C1_sum += (b[4, m] * acc_temp[m, :])
    C1 = ((v0 / h) - C1_sum + (acc_temp[4, :] / 2))
    s_0 = (C1 - (acc_temp[4, :] / 2))  # equal to C'1
    # Calculate S_0 from C1 & C2
    C2_sum = 0
    for m in range(9):
        C2_sum += (a[4, m] * acc_temp[m, :])
    C2 = ((r0 / (h ** 2)) - C2_sum + C1)
    S_0 = (C2 - C1)
    # Begin while loop for acceleration convergence
    while diff_acc >= accuracy:
        # 3b
        s_n = np.zeros((9, 3), dtype=float)
        S_n = np.zeros((9, 3), dtype=float)
        s_n[4, :] = s_0
        S_n[4, :] = S_0
        for n in range(1, 5):
            # 3b.i. Calculate s_n & S_n for n = -4...4, n != 0
            # Calculate s_n for 0 < n =< 4
            s_n[n + 4, :] = (s_n[n + 3, :] + ((acc_temp[n + 3, :] +
                                               acc_temp[n + 4, :]) / 2))
            # Calculate s_n for -4 <= n < 0
            s_n[-n + 4, :] = (s_n[-n + 5, :] - ((acc_temp[-n + 5, :] +
                                                 acc_temp[-n + 4, :]) / 2))
            # Calculate S_n for 0 < n =< 4
            S_n[n + 4, :] = (S_n[n + 3, :] + s_n[n + 3, :] +
                             (acc_temp[n + 3, :] / 2))
            # Calculate S_n for -4 <= n < 0
            S_n[-n + 4, :] = (S_n[-n + 5, :] - s_n[-n + 5, :] +
                              (acc_temp[-n + 5, :] / 2))
        # 3b.ii,iii. Calculate a_sum & b_sum using the arrays, a & b
        a_sum = np.zeros((9, 3), dtype=float)
        b_sum = np.zeros((9, 3), dtype=float)
        for n in range(9):
            # Calculate b_sum
            for k in range(9):
                b_sum[n, :] += (b[n, k] * acc_temp[k, :])
            # Calculate a_sum
            for k in range(9):
                a_sum[n, :] += (a[n, k] * acc_temp[k, :])
        # 3b.iv. Calculate radius and velocity for all n
        n_range_no_0 = np.array([0, 1, 2, 3, 5, 6, 7, 8])
        for n in n_range_no_0:
            vel_temp[n, :] = (h * (s_n[n, :] + b_sum[n, :]))
            rad_temp[n, :] = ((h ** 2) * (S_n[n, :] + a_sum[n, :]))
        # 3b.v. Evaluate acceleration using new radius and velocity values
        for j in range(9):
            acc_temp1[j, :] = -((rad_temp[j, :] *
                                 const.GRAV_PARAMATER['EARTH']) /
                                (np.linalg.norm(rad_temp[j, :]) ** 3))
        # 3c. Test convergence of accelerations in acc_int
        diff_acc = np.linalg.norm(acc_temp1 - acc_temp)
        acc_temp = dc(acc_temp1)
    # Set initial step in output arrays
    rad[0, :] = rad_temp[4, :]
    vel[0, :] = vel_temp[4, :]
    acc[0, :] = acc_temp[4, :]
    GD[0, :] = Epoch_GD0
    # Function to check what the actual frame is
    frame.append('J2000eq')

    # Predict
    step_current = 0
    for s in range(1, steps + 1):
        step_current += 1
        # 4. Calculate S_step
        S_step = (S_n[8, :] + s_n[8, :] + (acc_temp[8, :] / 2))
        # 5,6. Calculate b_sum_np1, a_sum_np1 with acceleration alterations
        b_sum_np1 = 0
        for k in range(9):
            b_sum_np1 += (b[9, k] * acc_temp[k, :])
        a_sum_np1 = 0
        for k in range(9):
            a_sum_np1 += (a[9, k] * acc_temp[k, :])
        # 7. Calculate rad_np1 & vel_np1
        vel_np1 = (h * (s_n[8, :] + (acc_temp[8, :] / 2) + b_sum_np1))
        rad_np1 = ((h ** 2) * (S_step + a_sum_np1))
        # Increment step in rad_temp/vel_temp to make room for n + 1 value
        vel_temp = np.roll(vel_temp, -1, axis=0)  # First in last position
        rad_temp = np.roll(rad_temp, -1, axis=0)
        # Replace rolled vector with rad_np1 and vel_np1
        vel_temp[8, :] = vel_np1
        rad_temp[8, :] = rad_np1

        # Evaluate - Correct
        # 8. Evaluate acc_step
        # Increment current step in acc_temp to make room for n + 1 value
        acc_temp = np.roll(acc_temp, -1, axis=0)
        # Calculate new acceleration for n + 1 and replace rolled vector
        acc_temp[8, :] = -(rad_np1 * (const.GRAV_PARAMATER['EARTH'] /
                                      (np.linalg.norm(rad_np1) ** 3)))
        # 9. Increment n (Incremented by rolling the arrays above)
        # 10. Converge radius and velocity vectors
        diff_rad = 1
        diff_vel = 1
        iteration = 0
        vel_step = vel_temp[8, :]
        rad_step = rad_temp[8, :]
        # Loop while radius and velocity calcs have not converged
        while diff_rad >= accuracy and diff_vel >= accuracy and iteration <= max_iterations:
            iteration += 1
            #print("Iteration = %s" % iteration)  # TODO: Remove
            # 10.a Calculate s_step and S_step for the n value
            s_step = (s_n[8, :] + ((acc_temp[7, :] + acc_temp[8, :]) / 2))
            S_step = (S_n[8, :] + s_n[8, :] + (acc_temp[8, :] / 2))
            # 10.b calculate b_sum_4 and a_sum_4 for first iteration only
            if iteration == 1:
                b_sum_4 = 0
                for k in range(8):
                    b_sum_4 += (b[8, k] * acc_temp[k, :])
                a_sum_4 = 0
                for k in range(8):
                    a_sum_4 += (a[8, k] * acc_temp[k, :])
            # 10.c Calculate b_sum_step & a_sum_step for the current step
            b_sum_step = (b[8, 8] * acc_temp[8, :])
            a_sum_step = (a[8, 8] * acc_temp[8, :])
            # 10.d Sum results of 10.a, 10.b, 10.c to get rad and vel step
            vel_step1 = (h * (s_step + b_sum_4 + b_sum_step))
            rad_step1 = ((h ** 2) * (S_step + (a_sum_4 + a_sum_step)))
            #print("Rad Step = %s" % rad_step1)  # TODO: Remove
            # 10.e Test convergence of rad and vel vectors
            diff_vel = np.linalg.norm(vel_step - vel_step1)
            diff_rad = np.linalg.norm(rad_step - rad_step1)
            # Replace step vectors with step1 vectors
            vel_step = dc(vel_step1)
            rad_step = dc(rad_step1)
            # Re-evaluate acc_temp before next loop
            if iteration < max_iterations:
                acc_temp[8, :] = -((rad_step * const.GRAV_PARAMATER['EARTH']) /
                                   (np.linalg.norm(rad_step) ** 3))
        if iteration > max_iterations:
            raise Exception('The maximum allowable number of iterations has been reached')
        # Increment step in S_n to make room for n + 1 value
        S_n = np.roll(S_n, -1, axis=0)
        s_n = np.roll(s_n, -1, axis=0)
        # Replace rolled vector with S_step & s_step
        S_n[8, :] = S_step
        s_n[8, :] = s_step
        # If all values converge
        rad[s, :] = rad_temp[4, :]
        vel[s, :] = vel_temp[4, :]
        acc[s, :] = acc_temp[4, :]
        frame.append('J2000eq')
        #######################################################################
        # Calculate gregorian date at new step and store in output array
        t_ms = np.zeros(7, dtype=float)
        t_ms[0:5] = GD_temp[0:5]
        t_ms[5] = np.floor(GD_temp[5])  # Seconds
        t_ms[6] = (np.mod(GD_temp[5], 1) * 1e6)  # Milliseconds
        t_ms = t_ms.astype(int)
        GD_dt = datetime(t_ms[0], t_ms[1], t_ms[2], t_ms[3],
                         t_ms[4], t_ms[5], t_ms[6])  # To datetime
        # Using h, add appropriate amount of time and calculate new date
        GD_add = timedelta(seconds=h)
        GD_stepped = GD_dt + GD_add
        # Convert back to numpy array and store
        GD_np = np.datetime64(GD_stepped)
        GD_new = np.zeros(6, dtype=float)
        GD_new[0] = GD_np.astype(object).year
        GD_new[1] = GD_np.astype(object).month
        GD_new[2] = GD_np.astype(object).day
        GD_new[3] = GD_np.astype(object).hour
        GD_new[4] = GD_np.astype(object).minute
        GD_new[5] = (GD_np.astype(object).second +
                     (GD_np.astype(object).microsecond * 1e-6))
        GD[s, :] = GD_new
        GD_temp = GD[s, :]

    # Outputs - Save all dictionary entries into tuple
    location_c = (GD, rad, vel, acc, frame)
    return location_c
