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
        # ## evaluation of vertical line v_2T
        # e_2 = phi(x_E) - x_2L
    else:
        ## evaluation of slop line l_2T
        e_1 = y_E - phi(y_L) - phi(alpha) * (x_E - phi(x_L))
        # ## evaluation of vertical line v_2T
        # e_2 = x_E - phi(x_2L)

    return E(x_2L, y_2L), e_1

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
        # ## evaluation of vertical line v_{T + P}
        # e_2 = phi(x_E) - x_LR
    else:
        ## evaluation of slop line l_{T + P}
        e_1 = y_E - phi(y_L) - phi(alpha) * (x_E - phi(x_L))
        # ## evaluation of vertical line v_{T + P}
        # e_2 = x_E - phi(x_LR)
    
    return E(x_LR, y_LR), e_1

## General Miller Loop Entry
def MillerLoop(P, Q, G, q, phi, reverse = False):
    ## if power q is negative or not
    P = P if q > 0 else -P
    q = q if q > 0 else -q
    
    T = P
    f1 = 1
    e_bits = [int(i) for i in bin(q)[2:]]

    print('Miller Loop Length: {}'.format(len(e_bits)))
    
    for i in range(1, len(e_bits)):
        ##### strip this specific manner used in Tate Pairing
        # if (i == len(e_bits) - 1) and (e_bits[i] == 0):
        #     f1 = f1 * (list(Q)[0] - list(T)[0])
        #     T = 2 * T
        #     break
        T, e_1 = double_line(T, Q, G, phi, reverse)
        f1 = f1 * f1 * e_1
        ##### strip this specific manner used in Tate Pairing
        # if (i == len(e_bits) - 1) and (e_bits[i] == 1):
        #     f1 = f1 * (list(Q)[0] - list(T)[0])
        #     T = T + P
        #     break
        if e_bits[i] == 1:
            T, e_1 = add_line(T, P, Q, G, phi, reverse)
            f1 = f1 * e_1
    
    return f1

## trival implementation of easy part, Frobenius not used here actually
## exp = (p^6 - 1) * (p^2 + 1)
## 2 * Frobenius + 2 * Mul + 1 * Inv
def easy_part(f):
    ff = f
    ## 1 * Frobenius
    t0 = f ** (p ** 6)
    ## 1 * Inv
    t1 = 1 / f
    ## 1 * Mul
    f = t0 * t1
    ## 1 * Frobenius
    t0 = f ** (p ** 2)
    ## 1 * Mul
    f = t0 * f
    
    actual = ff ** (((p ** 6) - 1) * ((p ** 2) + 1))
    assert(actual == f)
    
    return f

## reference from Algorithm 1 of "On the Computation of the Optimal Ate Pairing at the 192-bit Security Level"
## trival implementation of hard part, Frobenius not used here actually
## exp = (p^4 - p^2 + 1) / r
def hard_part(f, u, p, q):
    ## 1 * Sqr + 1 * Inv
    t0 = 1 / (f * f)
    ## 1 * Pow
    t5 = f ** u
    ## 1 * Sqr
    t1 = t5 * t5
    ## 1 * Mul
    t3 = t0 * t5
    
    ## 1 * Pow
    t0 = t3 ** u
    
    ## 1 * Pow
    t2 = t0 ** u
    
    ## 1 * Pow
    t4 = t2 ** u
    
    ## 1 * Mul
    t4 = t1 * t4
    ## 1 * Pow
    t1 = t4 ** u
    ## 1 * Inv
    t3 = 1 / t3
    ## 1 * Mul
    t1 = t3 * t1
    ## 1 * Mul
    t1 = t1 * f # f^\lambda_0
    
    # 1 * Inv
    t3 = 1 / f
    ## 1 * Mul
    t0 = t0 * f
    ## 1 * Frobenius
    t0 = t0 ** (p ** 3) # f^\lambda_3
    
    ## 1 * Mul
    t4 = t3 * t4
    ## 1 * Frobenius
    t4 = t4 ** p # f^\lambda_1
    
    ## 1 * Mul
    t5 = t2 * t5
    ## 1 * Frobenius
    t5 = t5 ** (p ** 2) # f^\lambda_2
    
    ## 3 * Mul
    t5 = t5 * t0
    t5 = t5 * t4
    t5 = t5 * t1
    
    ## third power of actual pairing result
    actual = f ** (((p ** 4) - (p ** 2) + 1) // q)
    assert(t5 == actual ** 3)
    assert(t5 ** q == 1)
    
    return t5

## Final Exponentiation Entry
def FinalExponentiation(f, p, k, q, u, trivial = True):
    if trivial:
        mu_r = f ** (((p ** k) - 1) // q)
    else:
        t0 = time.perf_counter()
        f = easy_part(f)
        t1 = time.perf_counter()
        mu_r = hard_part(f, u, p, q)
        t2 = time.perf_counter()
        print('\n     ##[Hard Part of Ate Pairing] Time consuming: t[easy] = {:.3f},  t[hard] = {:.3f}'.format(t1 - t0, t2 - t1))
    return mu_r

## Ate Pairing Entry
def AtePairing(P, Qx, G1, q, phi, p, k, u, T, trivial = True):
    t0 = time.perf_counter()
    f = MillerLoop(P, Qx, G1, T, phi, False)
    t1 = time.perf_counter()
    mu_r = FinalExponentiation(f, p, k, q, u, trivial)
    t2 = time.perf_counter()
    print('\n ##[Ate Pairing] Time consuming: t[f(P, Qx)] = {:.3f},  t[exp] = {:.3f}'.format(t1 - t0, t2 - t1))
    
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

## parameter for p(x), q(x), and t(x)
x = -15132376222941642752
t = x + 1
## p = ((x - 1)^2 * (x^4 - x^2 + 1)) / 3 + x
assert((pow((x - 1), 2) * (pow(x, 4) - pow(x, 2) + 1)) // 3 + x == p)
## q = x^4 - x^2 + 1
assert(pow(x, 4) - pow(x, 2) + 1 == q)
## t = x + 1
assert(abs(p + 1 - t) == Efp.order())

## p \equiv T \mod q
T = t - 1
####################################### Ate Pairing Testation
mu_r_ate = AtePairing(P, Qx, G1, q, phi, p, k, x, T, False)
assert(mu_r_ate ** q == Fp12(1))