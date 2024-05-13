import os

cur_path = os.getcwd()
load('{}/optimal_ate.sage'.format(cur_path))

def line_double(T):
    T = T.force_affine()
    assert(T.z == T.one_element())
    x, y = T.x, T.y
    ## slope: alpha = 3 * x^2 / 2 * y
    alpha = x.square().mul_scalar(Fp(3)).mul(y.mul_scalar(Fp(2)).inverse())
    bias = y - alpha * x

    ## projective coordinate
    # x, y, z = T.x, T.y, T.z
    # alpha = x.square().mul_scalar(Fp(3)).mul(y.mul(z).double().inverse())
    # bias = y.sub(alpha.mul(x).mul(z)).mul(z.square().mul(z).inverse())

    return alpha, bias

def line_add(T, P):
    T = T.force_affine()
    P = P.force_affine()
    assert(T.z == T.one_element())
    assert(P.z == P.one_element())
    x1, y1 = T.x, T.y
    x2, y2 = P.x, P.y
    ## slope: alpha = (y2 - y1) / (x2 - x1)
    alpha = y2.sub(y1).mul((x2.sub(x1)).inverse())
    ## bias: b = y1 - alpha * x1
    bias = y1.sub(alpha.mul(x1))
    # bias = y1 - alpha * x1

    return alpha, bias

################################################## cache line parameters for [6x + 2 + p - p^2 + p^3]Q
def line_function(Q, e, lamb):
    ############################################## double-add part, 6x + 2
    naf_digits = list(reversed(to_naf(e)))[1:]
    L = []
    T = Q
    for i, digit in enumerate(naf_digits):
        alpha, bias = line_double(T)
        T = T.double()
        L.append((alpha, bias))
        if digit^2 == 1:
            Qt = Q if digit == 1 else Q.negate()
            alpha, bias = line_add(T, Qt)
            T = T.add(Qt)
            L.append((alpha, bias))

    assert(T == Q.scalar_mul(e))

    ############################################# frobenius map part, p - p^2 + p^3
    # Q1 = pi(Q)
    # x = x' * beta^(2 * (p - 1) // 6)
    # y = y' * beta^(3 * (p - 1) // 6))
    pi_1_Q = G2(
        Q.x.conjugate_of().mul(Fp12.beta_pi_1[1]),
        Q.y.conjugate_of().mul(Fp12.beta_pi_1[2]),
        Fp2.ONE())
    assert(pi_1_Q.is_on_curve() == True)
    assert(pi_1_Q == Q.scalar_mul(px(x)))

    # Q2 = pi2(Q)
    # x = x * beta * (2 * (p^2 - 1) // 6)
    # y = y * beta * (3 * (p^2 - 1) // 6) = -y
    pi_2_Q = G2(
        Q.x.mul(Fp12.beta_pi_2[1]),
        Q.y.mul(Fp12.beta_pi_2[2]),
        Fp2.ONE())
    assert(pi_2_Q.is_on_curve() == True)
    assert(pi_2_Q == Q.scalar_mul(px(x) ** 2))

    # Q3 = pi3(Q)
    # x = x' * beta * (2 * (p^3 - 1) // 6)
    # y = y' * beta * (3 * (p^3 - 1) // 6)
    pi_3_Q = G2(
        Q.x.conjugate_of().mul(Fp12.beta_pi_3[1]),
        Q.y.conjugate_of().mul(Fp12.beta_pi_3[2]),
        Fp2.ONE())
    assert(pi_3_Q.is_on_curve() == True)
    assert(pi_3_Q == Q.scalar_mul(px(x) ** 3))

    alpha, bias = line_add(T, pi_1_Q)
    T = T.add(pi_1_Q)
    L.append((alpha, bias))

    assert(T == Q.scalar_mul(e + px(x)))

    alpha, bias = line_add(T, pi_2_Q.negate())
    T = T.add(pi_2_Q.negate())
    L.append((alpha, bias))

    k = e + px(x) - px(x) ** 2
    assert(T == Q.scalar_mul(k if k > 0 else rx(x) - (-k % rx(x))))

    # alpha, bias = line_add(T, pi_3_Q)
    alpha, bias = Fp2.ZERO(), T.x.mul(T.z.inverse().square())
    T = T.add(pi_3_Q)
    L.append((alpha, bias))

    assert(T == Q.scalar_mul(lamb))
    assert(T.is_infinite() == True)
    assert(alpha.is_zero() == True)

    return L

def test_line_precomputation():
    ####################################################
    ## assume we want to prove e(P1, Q1) = e(P2, Q2), namely e(P1, Q1) * e(P2, -Q2) = 1
    ## fixed point Q in G2, public known to verifier
    Q1 = g2.scalar_mul(1).force_affine()
    Q2 = g2.scalar_mul(3).force_affine()
    
    ## point P sent from prover
    P1 = g1.scalar_mul(3).force_affine()
    P2 = g1.scalar_mul(1).force_affine()
    
    ## indexer (oracle) for fixed point Q
    e = 6 * x + 2
    lamb = lambdax(x)
    print('parameter x = {}\n'.format(x))
    L1 = line_function(Q1.force_affine(), e, lamb)
    L2 = line_function(Q2.negate().force_affine(), e, lamb)
    print('\n [Oracle] line function for Q1, and Q2 are both precomputated. \n')

# test_line_precomputation()