# Utilized Modules
import numpy as np
from orbital_analyses import u
from copy import deepcopy as dc
import Requirements as Req

from orbital_analyses.Transform_State import Gregorian2JD

# Constants
E_Rad = 6378.137            # Earth Radius (km))
E_Mu = 398600.4418          # Earth Mu (km^3 / sec^2)

# %% Determine Inital State
# Initial Rad / Vel - R0, V0 and epoch
# r0 = (Req.Rad).magnitude  # km
r0_i = np.array([[1949.73850397223],
                 [-5139.978114855865],
                 [-4852.90271319306]])  # km
# v0 = (Req.Vel).magnitude  # km / sec
v0_i = np.array([[-6.574850967942473],
                 [0.2153270408156436],
                 [-3.149950105296086]])  # km / sec
# EpochGD_0 = Req.GD_UTC  # Gregorian Date
EpochGD_0 = np.array([[2004], [4], [6], [0], [0], [0]])  # Gregorian Date
EpochJD_0 = Gregorian2JD(EpochGD_0)  # Julian Date
# Step = Req.Ssize
step = 30 * u.sec
steps = 25

# %%###########################################################################
GD = dc(EpochGD_0)


def Gauss_Jackson_Prop(r0_i, v0_i, GD, step, steps, max_iterations=50):
    # An 8th order integrator corrector for propagation
    # 1. Use f & g series to calculate 8 (-4,4) rad/vel surrounding the epoch_0
    # Initialize vectors to place all 9 initial vectors into
    Rad_init = np.zeros((3, 9), dtype=float)
    Vel_init = np.zeros((3, 9), dtype=float)
    Acc_init = np.zeros((3, 9), dtype=float)
    # TODO: Add in handling if rad and vel given as 3x1
    r0 = (r0_i).reshape(1, 3)
    v0 = (v0_i).reshape(1, 3)
    Rad_init[:, 4] = np.squeeze(r0)
    Vel_init[:, 4] = np.squeeze(v0)
    # Determine all other variables needed for loop
    r0_mag = np.linalg.norm(r0)
    r0_2 = np.linalg.norm(np.inner(r0, r0))  # km^2
    v0_2 = np.linalg.norm(np.inner(v0, v0))  # km^2 / sec^2
    # TODO: Add coditional if step has no units
    StepSec = ((step.to(u.sec)).magnitude)
    # set initial step value
    h = 1  # TODO: confirm this isn't 0 to start

    ###########################################################################
    # TODO: Move all function definitions to top of main function or section?
    def guess1(Tdelt, alpha):
        # Vallado Alg. 8, Pg. 93
        X0 = (np.sqrt(E_Mu) * Tdelt * alpha)
        return X0

    def guess2(r0, v0, Tdelt):
        # Vallado Alg. 8, Pg. 93
        h0 = np.cross(r0, v0)
        h0_2 = np.linalg.norm(np.inner(h0, h0))
        p = (h0_2 / E_Mu)
        s = ((np.arccos(3 * np.sqrt(E_Mu / (p ** 3)) * Tdelt)) / 2)
        w = np.arctan(np.cbrt(np.tan(s)))
        X0 = (np.sqrt(p) * 2 * np.cot(2 * w))
        return X0

    def guess3(r0, v0, Tdelt, alpha, r0_mag):
        # Vallado Alg. 8, Pg. 93
        a = (1 / alpha)
        X0 = (np.sign(Tdelt) * np.sqrt(-a) *
              np.log((-2 * E_Mu * alpha * Tdelt) /
                     ((np.linalg.norm(np.inner(r0, v0))) + np.sign(Tdelt) *
                      np.sqrt(-E_Mu * a) * (1 - r0_mag * alpha))))
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

    # Newton-Raphson iteration using Universal Variabes to calculate
    # the 8 (-4,4) rad/vel surrounding the epoch_0
    # Vallado - ALg. 8, pg. 93

    # loop to calculate rad/vel vectors of 8 other positions defined by Step
    n_range = np.linspace(-4, 4, 9, dtype=int)
    for i in n_range:
        # Compute the Rad and Vel vectors at each step around the origin
        Tdelt = i * StepSec
        M_Enrgy0 = ((v0_2 / 2) - (E_Mu / r0_mag))  # TODO: Delete?
        alpha0 = ((-v0_2 / E_Mu) + (2 / r0_mag))
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
                                     np.sqrt(E_Mu)) * Xn * (1 - (Psi * c3))) +
                 (r0_mag * (1 - (Psi * c2))))
            Xn1 = (Xn + ((((np.sqrt(E_Mu) * Tdelt) -
                           ((Xn ** 3) * c3) -
                           ((np.linalg.norm(np.inner(r0, v0)) /
                             np.sqrt(E_Mu)) * ((Xn ** 2) * c2)) -
                           (r0_mag * Xn * (1 - (Psi * c3))))) / r))
            Xn = dc(Xn1)
            err = np.absolute(Xn - Xn1)
        # With iterated values, calculate relevant f and g functions
        f = (1 - ((Xn ** 2) / (r0_mag)) * c2)
        g = (Tdelt - ((Xn ** 3) / np.sqrt(E_Mu)) * c3)
        f_dot = ((np.sqrt(E_Mu) / (r * r0_mag)) * Xn * ((Psi * c3) - 1))
        g_dot = (1 - ((Xn ** 2) / (r)) * c2)
        # Check for correctness
        # print((f * g_dot) - (f_dot * g))  # TODO: Delete later (Should = 1)
        # Calculate new rad / vel vectors
        rN = ((f * r0) + (g * v0))
        vN = ((f_dot * r0) + (g_dot * v0))
        # Place rad / vel vectors into matrices
        if i == 0:
            Rad_init[:, 4] = np.squeeze(r0)
            Vel_init[:, 4] = np.squeeze(v0)
        else:
            Rad_init[:, i + 4] = np.squeeze(rN)
            Vel_init[:, i + 4] = np.squeeze(vN)

    # %% 2. Evaluate the 9 acceeration vectors for the 9 states determined above
    # Acceleration w/o any perturbations
    # 2 body propagation
    for j in range(9):
        Acc_init[:, j] = -(Rad_init[:, j] *
                           (E_Mu / (np.linalg.norm(Rad_init[:, j]) ** 3)))
    # Incorporate Special Perturbations into acceleration model
    # Vallado Alg. 64, Pg. 591
    # TODO: Gather more accurate models of perturbative forces
        # https://sourceforge.net/p/gmat/git/ci/GMAT-R2018a/tree/application/data/
        # https://sourceforge.net/p/gmat/git/ci/GMAT-R2018a/tree/prototype/StateConv/
        # http://www.hayabusa.isas.jaxa.jp/kawalab/dromobile/Papers/HernandoAyuso%20Dromobile%202016.pdf
        # https://www.sciencedirect.com/science/article/pii/0898122186900258

    # %% 3. Converging the accelerations
    # Import all coefficient arrays
    # Eigth order summed adams coefficients in Ordinate Form: a(j,k), b(j,k)
    a = np.load(r'orbital_analyses\coefficient_matrices\gauss_jackson_prop\GJ_a_coeff.npy')
    b = np.load(r'orbital_analyses\coefficient_matrices\gauss_jackson_prop\GJ_b_coeff.npy')

    # Initialize all loop arrays and variables before loop
    acc_int = dc(Acc_init)
    acc_int1 = np.zeros((3, 9), dtype=float)
    diff_acc = 1

    # 3a. Calculate s_0 and S_0
    # Calculate s_0 from C1
    C1_sum = 0
    for m in range(9):
        C1_sum += (b[4, m] * acc_int[:, m])
    C1 = ((v0 / h) - C1_sum + (acc_int[:, 4] / 2))
    s_0 = (C1 - (acc_int[:, 4] / 2))  # equal to C'1
    # Calculate S_0 from C1 & C2
    C2_sum = 0
    for m in range(9):
        C2_sum += (a[4, m] * acc_int[:, m])
    C2 = ((r0 / (h ** 2)) - C2_sum + C1)
    S_0 = (C2 - C1)
    # Begin while loop for acceleration convergence
    while diff_acc != 0:
        acc_int1 = np.zeros((3, 9), dtype=float)
        # 3b
        s_n = np.array(np.zeros([3, 9]))
        S_n = np.array(np.zeros([3, 9]))
        s_n[:, 4] = np.squeeze(s_0)
        S_n[:, 4] = np.squeeze(S_0)
        for n in range(1, 5):
            # 3b.i. Calculate s_n & S_n for n = -4...4, n != 0
            # Calculate s_n for 0 < n =< 4
            s_n[:, n + 4] = (s_n[:, n + 3] +
                             np.transpose((acc_int[:, n + 3] +
                                           acc_int[:, n + 4]) / 2))
            # Calculate s_n for -4 <= n < 0
            s_n[:, -n + 4] = (s_n[:, -n + 5] -
                              np.transpose((acc_int[:, -n + 5] +
                                            acc_int[:, -n + 4]) / 2))
            # Calculate S_n for 0 < n =< 4
            S_n[:, n + 4] = (S_n[:, n + 3] + s_n[:, n + 3] +
                             np.transpose(acc_int[:, n + 3] / 2))
            # Calculate S_n for -4 <= n < 0
            S_n[:, -n + 4] = (S_n[:, -n + 5] - s_n[:, -n + 5] +
                              np.transpose(acc_int[:, -n + 5] / 2))
        # 3b.ii,iii. Calculate a_sum & b_sum using the arrays, a & b
        a_sum = np.array(np.zeros([3, 9]))
        b_sum = np.array(np.zeros([3, 9]))
        for n in range(9):
            # Calculate b_sum
            for m in range(9):
                b_sum[:, n] += np.squeeze(np.array(b[n, m] * acc_int[:, m]))
            # Calculate a_sum
            for m in range(9):
                a_sum[:, n] += np.squeeze(np.array(a[n, m] * acc_int[:, m]))
        # 3b.iv. Calculate radius and velocity for all n
        h = 1  # h is the step (t_n = t_0 + nh) pg. 334
        rad_int = np.array(np.zeros([3, 9]))
        vel_int = np.array(np.zeros([3, 9]))
        for n in range(9):
            vel_int[:, n] = (h * (s_n[:, n] + b_sum[:, n]))
            rad_int[:, n] = ((h ** 2) * (S_n[:, n] + a_sum[:, n]))
        vel_int[:, 4] = np.transpose(Vel_init[:, 4])
        rad_int[:, 4] = np.transpose(Rad_init[:, 4])
        # 3b.v. Evaluate acceleration using new radius and velocity values
        for j in range(9):
            acc_int1[:, j] = -((rad_int[:, j]) *
                               (E_Mu / (np.linalg.norm(rad_int[:, j]) ** 3)))
        # 3c. Test convergence of accelerations in acc_int
        diff_acc = np.linalg.norm(acc_int1 - acc_int)
        acc_int = dc(acc_int1)

    # %% Predict
    # Dictionary entries to save all steps into
    rad = np.zeros((steps + 1, 3), dtype=float)
    vel = np.zeros((steps + 1, 3), dtype=float)
    acc = np.zeros((steps + 1, 3), dtype=float)
    frame = []
    # Initial step
    rad[0, :] = np.squeeze(rad_int[:, 4])
    vel[0, :] = np.squeeze(vel_int[:, 4])
    acc[0, :] = np.squeeze(acc_int[:, 4])
    frame.append('ECEF')

    # 4. Calculate S_np1 (S_5?)
    S_np1 = (S_n[:, 4] + s_n[:, 4] + (acc_int[:, 4] / 2))
    # 5,6. Calculate b_sum_4, a_sum_4 with acceleration alterations
    b_sum_4 = 0
    for m in range(9):
        b_sum_4 += (b[8, m] * acc_int[:, m])
    a_sum_4 = 0
    for m in range(9):
        a_sum_4 += (a[8, m] * acc_int[:, m])
    # 7. Calculate rad_np1 & vel_np1
    vel_np1 = (h * (s_n[:, 4] + (acc_int[:, 4] / 2) + b_sum_4))
    rad_np1 = ((h ** 2) * (S_np1 + a_sum_4))

    # %% Evaluate - Correct
    # 8. Evaluate acc_step
    acc_step = np.zeros((3, 9), dtype=float)
    # Import calculated accelerations from acc_int
    acc_step[:, 0:4] = dc(acc_int[:, 1:5])
    acc_step[:, 5:8] = dc(acc_int[:, 6:9])
    # Calculate new acceleration for n + 1 (even though it's [4] in the array)
    acc_step[:, 4] = -np.squeeze(rad_np1 *
                                 (E_Mu / (np.linalg.norm(rad_np1) ** 3)))
    # 9. Increment n
    # n = 5

    # 10.
    # While diff_rad != 0 && diff_vel != 0:
    # TODO: Max number of corrector iterations here
    diff_rad = 1
    diff_vel = 1
    iteration = 0
    max_iterations = np.inf
    # Loop while radius and velocity calcs have not converged
    while diff_rad != 0 and diff_vel != 0 and iteration <= max_iterations:
        iteration += 1
        # 10.a Calculate s_n_step
        C1_sum = 0
        for m in range(9):
            C1_sum += (b[4, m] * acc_int[:, m])
        C1 = ((v0 / h) - C1_sum + (acc_int[:, 4] / 2))
        s_n_step = (C1 - (acc_int[:, 4] / 2))  # equal to C'1

        if iteration == 1:
            # 10.b calculate b_sum and a_sum
            b_sum_4 = 0
            for m in range(9):
                b_sum_4 += (b[8, m] * acc_int[:, m])
            a_sum_4 = 0
            for m in range(9):
                a_sum_4 += (a[8, m] * acc_int[:, m])


    # TODO: Dont forget about this per each loop at the end (edited of course)
    frame.append('ECEF')
    step += 1
    rad[0, :] = np.squeeze(rad_int[:, 4])
    vel[0, :] = np.squeeze(vel_int[:, 4])
    acc[0, :] = np.squeeze(acc_int[:, 4])

    # Outputs - Save all dictionary entries into tuple
    location_c = (rad, vel, acc, frame)
    return location_c



















# %%###########################################################################

    # Initialize all loop arrays and variables before loop
    h = 1  # h is the step (t_n = t_0 + nh) pg. 334
    s_n = np.array(np.zeros([3, 8]))
    S_n = np.array(np.zeros([3, 8]))
    rad_int = np.array(np.zeros([3, 9]))
    vel_int = np.array(np.zeros([3, 9]))
    acc_int_magnitude = np.array(np.zeros([1, 9]))
    acc_int = dc(Acc_init)
    diff = 1
    acc_error = 1

    # Loop to converge accelerations
    while diff != 0:
        acc_error1 = dc(acc_error)
        # Calculate s_0 from C1
        C1 = ((v0_i / h) +
              np.transpose(np.transpose(acc_int[:, 4] / 2) -
                           np.dot(b[4, :], np.transpose(acc_int))))
        s_n[:, 4] = np.transpose(C1 - (acc_int[:, 4] / 2))  # equal to C'1
        # Calculate s_n for -4 <= n <= 4
        for n in range(1, 5):
            # Positive n value
            s_n[:, n + 4] = (s_n[:, n + 3] +
                             np.transpose((acc_int[:, n + 3] +
                                           acc_int[:, n + 4]) / 2))
            # Negative n value
            s_n[:, -n + 4] = (s_n[:, -n + 5] -
                              np.transpose((acc_int[:, -n + 5] +
                                            acc_int[:, -n + 4]) / 2))
        # Calculate S_0 from C2 and C1
        C2 = ((r0_i / (h ** 2)) + (C1) -
              np.transpose(a[4, :] * np.transpose(acc_int)))
        S_n[:, 4] = np.transpose(C2 - C1)
        # Calculate S_n for -4 <= n <= 4
        for n in range(1, 5):
            # Positive n value
            S_n[:, n + 4] = (S_n[:, n + 3] + s_n[:, n + 3] +
                             np.transpose(acc_int[:, n + 3] / 2))
            # Negative n value
            S_n[:, -n + 4] = (S_n[:, -n + 5] - s_n[:, -n + 5] +
                              np.transpose(acc_int[:, -n + 5] / 2))
        # Calculate redius and velocity again using new s0 & S0 values
        # This is used to test the convergence of the accelerations
        for n in range(9):
            vel_int[:, n] = (h * (s_n[:, n] +
                                  (b[n, :] * np.transpose(acc_int))))
            rad_int[:, n] = ((h ** 2) * (S_n[:, n] +
                                         (a[n, :] * np.transpose(acc_int))))
        # Evaluate the updated acceleration using appropriate force models
        for j in range(9):
            acc_int[:, j] = -np.reshape(rad_int[:, j] *
                                        (E_Mu /
                                         (np.linalg.norm(rad_int[:, j]) **
                                          3)), (3, -1))
        # Determine error from acc_0
        for i in range(9):
            acc_int_magnitude[0, i] = np.linalg.norm(acc_int[:, i])
        acc_error = np.linalg.norm(acc_int_magnitude - acc_int_magnitude[0, 4])
        diff = acc_error1 - acc_error

    # %% 4. Predictor
    # Calculate S_n+1 with new, converged accelerations (S_5?) n=5?
    S_5 = np.reshape(S_n[:, 8] + s_n[:, 8] +
                     np.transpose(acc_int[:, 8] / 2), (3, -1))
    vel_5 = np.reshape(h * (s_n[:, 8] + np.transpose(acc_int[:, 8] / 2) +
                            (b[9, :] * np.transpose(acc_int))), (3, -1))
    rad_5 = np.reshape((h ** 2) * (S_5 + np.transpose(a[9, :] *
                                   np.transpose(acc_int))), (3, -1))





























