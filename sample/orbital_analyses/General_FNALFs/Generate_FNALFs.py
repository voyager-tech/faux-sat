# Utilized modules
import numpy as np
from sympy import Symbol, diff, S
import pickle

eg_n = np.linspace(2.0, 200, num=199, dtype=np.int)
eg_fun = np.asmatrix(np.zeros((20298, 1), dtype=np.float))
eg_dif = np.asmatrix(np.zeros((20298, 1), dtype=np.float))

eg_fun_l = []
eg_dif_l = []

# Symbolic variable defintion
u = Symbol('u')

# Initial Values for recursion
# 0, 0
eg_fun_l.append(S(1))
eg_dif_l.append(S(0))

# 1, 0
eg_fun_l.append(S(np.sqrt(3)))  # Confirm Value
eg_dif_l.append(S(0))

# 1, 1
eg_fun_l.append(S(np.sqrt(3) * u))
eg_dif_l.append(S(np.sqrt(3)))

# Holmes & Featherstone
# http://mitgcm.org/~mlosch/geoidcookbook/node11.html
# Legrende Function for each n, m combination
for n in eg_n:
    print("n=%s" % n)
    for m in range(n + 1):

        if n > m:  # Non-sectorial
            # Other values to calcualate
            eg_Anm = np.sqrt((((2 * n) - 1) * ((2 * n) + 1)) /
                             ((n - m) * (n + m)))
            eg_Bnm = np.sqrt((((2 * n) + 1) * (n + m - 1) *
                              (n - m - 1)) /
                             ((n - m) * (n + m) * ((2 * n) - 3)))
            if (n - m) == 1:
                n1m0 = (((((n - 1) ** 2) + (n - 1)) / 2) + m)
                P_n1m0 = eg_fun_l[int(n1m0)]
                P_n2m0 = 0
            else:
                n1m0 = (((((n - 1) ** 2) + (n - 1)) / 2) + m)
                P_n1m0 = eg_fun_l[int(n1m0)]
                n2m0 = (((((n - 2) ** 2) + (n - 2)) / 2) + m)
                P_n2m0 = eg_fun_l[int(n2m0)]
            P_nm = ((eg_Anm * P_n1m0) - (eg_Bnm * P_n2m0))
#            print("P_nm=%s" % P_nm)

        if n == m:  # Sectorial
            mult = np.asmatrix(np.zeros((1, m - 1), dtype=np.float64))
            for j in range(2, m + 1):
                mult[0, (j - 2)] = np.sqrt(((2 * j) + 1) / (2 * j))
            P_nm = ((u ** m) * (np.sqrt(3)) * np.prod(mult, axis=1))
#            print("P_nm=%s" % P_nm)

        eg_fun_l.append(S(P_nm))
        P_diff = diff(P_nm, u)
#        print("P_diff=%s" % P_diff)
        eg_dif_l.append(P_diff)

###############################################################################

    # Import Export statement examples
    with open('General_FNALFs/DO_200', 'wb') as fp:
        pickle.dump(eg_fun_l, fp)

    with open('General_FNALFs/DO_200_dif', 'wb') as fp:
        pickle.dump(eg_dif_l, fp)

    with open('General_FNALFs/DO_200', 'rb') as fp:
        DO_200 = pickle.load(fp)







###############################################################################
###############################################################################































## slideheaven.com_the-integral-formulas-of-the-associated-legendre-f.pdf
## Legrende Function for each n, m combination
#    for n in eg_n:
#        print(n)
#        for m in range(n + 1):
#            print(m)
#            for i in range(eg_steps):
#                # New Constants for each change in m
#                x = Symbol('x', real=True)
#                T = Symbol('T', positive=True)  # Theta
#                L = Symbol('L')  # Lambda
#
#                eg_s = np.sin(eg_RLL[1, i] * (np.pi / 180))  # geocentric lat
#
#                # Determine how far along degree/order you are
#                DO = ((((n ** 2) + n) / 2) + m - 3)  # Refers to 0 index
#
#                # Other values to calcualate
#                eg_rho = sqrt(x**2 + (-x**2 + 1)*sin(T)**2)
#                eg_alpha = atan(sqrt(-x**2 + 1)*sin(T)/x)
#
##                eg_rho = (sqrt((sin(x) ** 2) + (1 - (sin(x) ** 2)) * (sin(T) ** 2)))
##                eg_alpha = atan((sqrt(1 - (sin(x) ** 2)) / (sin(x))) * sin(T))
##                eg_fun1 = ((eg_rho ** n) * cos(n * eg_alpha))
#
#                # Calculate Pn
#                eg_fun1 = simplify((eg_rho ** n) * cos(n * eg_alpha))
##                eg_fun1 = eg_fun1.subs(x, eg_s)
#                eg_Pn = simplify((1 / pi) * (integrate(eg_fun1, (T, 0, pi))))
#
#                # Calculate P_nm (0)
#                if np.mod((n - m), 2) != 0:  # if n - m is odd, P_nm(0) = 0
#                    eg_P_nm0 = 0
#                if np.mod((n - m), 2) == 0:  # if n - m is even, find P_nm(0)
#                    k = ((n - m) / 2).astype(int)
#                    if k == 0:
#                        Yk = 1
#                    if k == 1:
#                        Yk = (np.sqrt(2) / 2)
#                    if k >= 2:
#                        Yk = (np.sqrt(2) / 2)
#                        for j in range(2, k + 1):
#                            Yk_1 = (np.sqrt((2 * j + 1) / (2 * j + 2)) * Yk)
#                            Yk = Yk_1
#                    if (k + m) == 1:
#                        Ykm = (np.sqrt(2) / 2)
#                    if (k + m) >= 2:
#                        Ykm = (np.sqrt(2) / 2)
#                        for j in range(2, k + m + 1):
#                            Ykm_1 = (np.sqrt((2 * j + 1) / (2 * j + 2)) * Ykm)
#                            Ykm = Ykm_1
#
#                eg_P_nm0 = ((-1 ** k) * np.sqrt((2 * n + 1) / 2) * Yk * Ykm)
#
#                # Calculate P_nm (x) for n-m = even number
#                P_nmX1 = (eg_Pn * (sqrt(1 - (x ** 2)) * cos(L)) * cos(m * L))
##                P_nmX1 = P_nmX1.subs(x, eg_s)
#                P_nmX2 = (integrate(P_nmX1, (L, 0, (pi / 2))))
#                eg_P_nmX = ((1 / pi) * ((2 * n + 1) / (eg_P_nm0)) * P_nmX2)
#                print(eg_P_nmX)
#                eg_Lege_Fun.append(eg_P_nmX)
#
#                eg_P_nmX1 = simplify(eg_P_nmX)
#                eg_P_nmD = (integrate(eg_P_nmX1, x))
#                print(eg_P_nmD)
#                eg_Lege_dif.append(eg_P_nmD)
#
#                # Calculate P_nm (x) for n-m = odd
#
#
#                # Calculate Normalized Associated Legendre Function
#                if m == 0:
#                    k = 1
#                if m != 0:
#                    k = 0
#
#                eg_Lege_Fun.append(eg_mid_Fun)
#                eg_Lege_dif.append(eg_mid_dif)
#
#    eg_fun = np.zeros((2333877, 1), dtype=np.float)
#    eg_dif = np.zeros((2333877, 1), dtype=np.float)
#
################################################################################
#
#    # Legrende Function for each n, m combination
#    for n in eg_n:
#        print(n)
#        for m in range(n + 1):
#            print(m)
#            for i in range(eg_steps):
#                # New Constants for each change in m
#                s = Symbol('s')
#
#                # Determine how far along degree/order you are
#                DO = ((((n ** 2) + n) / 2) + m - 3)  # Refers to 0 index
#
#                # Calculate Legendre Polynomial defined by n
#                eg_mid_n = diff((((s ** 2) - 1) ** n), s, n)
#                eg_Lege_n = ((1 / ((2 ** n) * fac(n))) * eg_mid_n)
#
#                # Calculate Legendre Polynomial defined by m
#                eg_mid_nm = diff(eg_Lege_n, s, m)
#                eg_Lege_nm = (((1 - (s ** 2)) ** (m / 2)) * eg_mid_nm)
#
#                # Calculate Normalized Associated Legendre Function
#                if m == 0:
#                    k = 1
#                if m != 0:
#                    k = 0
#                eg_mid_Fun = ((np.sqrt(((2 * n) + 1) / (2 * (1 + k)))) *
#                              (np.sqrt((fac(n - m)) / (fac(n + m)))) *
#                              eg_Lege_nm)
#                eg_mid_dif = diff(eg_mid_Fun, s)
#                eg_Lege_Fun.append(eg_mid_Fun)
#                eg_Lege_dif.append(eg_mid_dif)
#
#    eg_fun = np.zeros((2333877, 1), dtype=np.float)
#    eg_dif = np.zeros((2333877, 1), dtype=np.float)
#
#    for j in range(eg_steps):
#        # Define symbols
#        s = Symbol('s')
#        c = Symbol('c')
#
#        # Determine sin & cos
#        eg_s = np.sin(eg_RLL[1, i] * (np.pi / 180))  # geocentric latitude
#
#        eg_fun1 = eg_Lege_Fun[j]
#        eg_dif1 = eg_Lege_dif[j]
#
#        eg_fun_s = eg_fun1.sub(s, eg_s)
#        eg_dif_s = eg_dif1.sub(s, eg_s)
#
#        eg_fun[j, 0] = np.float(eg_fun_s)
#        eg_dif[j, 0] = np.float(eg_dif_s)
#
################################################################################
#
#    # Legrende Function for each n, m combination
#    for n in eg_n:
#        for m in range(n + 1):
#            for i in range(eg_steps):
#                # New Constants for each change in m
#                s = Symbol('s')
#
#                eg_s = np.sin(eg_RLL[1, i] * (np.pi / 180))  # (geocentric latitude) in rads
#
#                # Determine how far along degree/order you are
#                DO = ((((n ** 2) + n) / 2) + m - 3)  # Refers to 0 index
#
#                # Calculate Legendre Polynomial defined by n
#                eg_mid_n = diff((((s ** 2) - 1) ** n), s, n)
#                eg_Lege_n = ((1 / ((2 ** n) * fac(n))) * eg_mid_n)
#
#                # Calculate Legendre Polynomial defined by m
#                eg_mid_nm = diff(eg_Lege_n, s, m)
#                eg_Lege_nm = ((c ** m) * eg_mid_nm)
#
#                # Calculate Normalized Associated Legendre Function
#                if m == 0:
#                    k = 1
#                if m != 0:
#                    k = 2
#                eg_mid_Fun = ((np.sqrt((fac(n - m) * ((2 * n) + 1) * k) /
#                                       (fac(n + m)))) * eg_Lege_nm)
#                eg_mid_dif = diff(eg_mid_Fun, s)
#
#                eg_mid_Fun1 = (eg_mid_Fun.subs(c, eg_c))
#                eg_mid_dif1 = (eg_mid_dif.subs(c, eg_c))
#                eg_mid_Fun2 = (eg_mid_Fun1.subs(s, eg_s))
#                eg_mid_dif2 = (eg_mid_dif1.subs(s, eg_s))
#                eg_Lege_Fun[DO, i] = (eg_mid_Fun2)
#                eg_Lege_dif[DO, i] = (eg_mid_dif2)

    








c = 0
for a in eg_n:
    for b in range(a + 1):
        c = c + 1
print(c)
        







