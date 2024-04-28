import time

## evaluation of double line divisor function
## arithmetics on fields, not on multiplicative group
def double_line(line_point, eval_point, E, phi, reverse = False):
    ######################## arithemtic on finite field of line_point
    ## lambda = 3x^2 / 2y
    (x_L, y_L) = (list(line_point)[0], list(line_point)[1])
    (x_E, y_E) = (list(eval_point)[0], list(eval_point)[1])
    alpha = (3 * x_L^2) / (2 * y_L)
    x_2L = alpha * alpha - 2 * x_L
    y_2L = -y_L - alpha * (x_2L - x_L)
    
    ######################## arithmetic on mixed finite field
    ## x_E, y_E \in F2
    ## y_L, x_L, alpha, x_2L \in F1
    if reverse:
        ## evaluation of slop line l_2T
        e_1 = phi(y_E) - y_L - alpha * (phi(x_E) - x_L)
        ## evaluation of vertical line v_2T
        e_2 = phi(x_E) - x_2L
    else:
        ## evaluation of slop line l_2T
        e_1 = y_E - phi(y_L) - phi(alpha) * (x_E - phi(x_L))
        ## evaluation of vertical line v_2T
        e_2 = x_E - phi(x_2L)

    return E(x_2L, y_2L), e_1, e_2

## evaluation of add line divisor function
## arithmetics on fields, not on multiplicative group
def add_line(line_left_point, line_right_point, eval_point, E, phi, reverse = False):
    ######################## arithemtic on finite field of line_point
    ## lambda = (y2 - y1) / (x2 - x1)
    (x_L, y_L) = (list(line_left_point)[0], list(line_left_point)[1])
    (x_R, y_R) = (list(line_right_point)[0], list(line_right_point)[1])
    (x_E, y_E) = (list(eval_point)[0], list(eval_point)[1])
    alpha = (y_L - y_R) / (x_L - x_R)
    x_LR = alpha * alpha - x_L - x_R
    y_LR = -y_L - alpha * (x_LR - x_L)
    
    ######################## arithmetic on mixed finite field
    ## x_E, y_E \in F2
    ## y_L, x_L, alpha, x_LR \in F1
    if reverse:
        ## evaluation of slop line l_{T + P}
        e_1 = phi(y_E) - y_L - alpha * (phi(x_E) - x_L)
        ## evaluation of vertical line v_{T + P}
        e_2 = phi(x_E) - x_LR
    else:
        ## evaluation of slop line l_{T + P}
        e_1 = y_E - phi(y_L) - phi(alpha) * (x_E - phi(x_L))
        ## evaluation of vertical line v_{T + P}
        e_2 = x_E - phi(x_LR)
    
    return E(x_LR, y_LR), e_1, e_2

## General Miller Loop Entry
def MillerLoop(P, Q, G, q, phi, reverse = False):
    T = P
    f1 = 1
    f2 = 1
    e_bits = [int(i) for i in bin(q)[2:]]
    ## last bit cannot be evaluated, since the slope would be a vertical line
    for i in range(1, len(e_bits)):
        if (i == len(e_bits) - 1) and (e_bits[i] == 0):
            f1 = f1 * (list(Q)[0] - list(T)[0])
            T = 2 * T
            break
        T, e_1, e_2 = double_line(T, Q, G, phi, reverse)
        f1, f2 = (f1 * f1 * e_1, f2 * f2 * e_2)
        if (i == len(e_bits) - 1) and (e_bits[i] == 1):
            f1 = f1 * (list(Q)[0] - list(T)[0])
            T = T + P
            break
        if e_bits[i] == 1:
            T, e_1, e_2 = add_line(T, P, Q, G, phi, reverse)
            f1, f2 = (f1 * e_1, f2 * e_2)
    assert(T == G(0))
    
    return f1 / f2

## Weil Pairing Entry
def WeilPairing(P, Qx, G1, G12, q, phi):
    t0 = time.perf_counter()
    f_rP_Q = MillerLoop(P, Qx, G1, q, phi, False)
    t1 = time.perf_counter()
    f_rQ_P = MillerLoop(Qx, P, G12, q, phi, True)
    t2 = time.perf_counter()
    mu_r = ((-1) ** q) * (f_rP_Q / f_rQ_P)
    print('\n ##[Weil Pairing] Time consuming: t[f(P, Qx)] = {:.3f},  t[f(Qx, P)] = {:.3f}'.format(t1 - t0, t2 - t1))
    
    return mu_r


G1, G2_t, G12, G12_t = (Efp, Efp2_t, Efp12, Efp12_t)
C1, C2 = (cofactor_E1, cofactor_E2_t)

## make sure they are in G1 and G2_t repectively
P, Q = (C1 * G1.random_element(), C2 * G2_t.random_element())
assert(q * P == G1(0))
assert(q * Q == G2_t(0))

## untwist from E2_t to E12: Q -> Qx
Qx = into_E12(Q, beta, Fp, w, beta_t_x, beta_t_y, G12)
assert(q * Qx == G12(0))
assert(trace_map(Qx, 12, p, G12) == G12(0))

####################################### Weil Pairing Testation 
## P is defined over E(Fp), Qx is defined over E(Fpk)
## phi maps Fp to Fp12
phi = Hom(Fp, Fp12)(Fp.gen().minpoly().roots(Fp12)[0][0])
assert(P.curve() is not Qx.curve())
mu_r_weil = WeilPairing(P, Qx, G1, G12, q, phi)
## make sure pairing result is in q-torsion subgroup
assert(mu_r_weil ** q == Fp12(1))
#######################################