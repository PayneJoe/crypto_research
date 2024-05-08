import os
import copy

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

    return alpha, bias

def line_function(Q, e):
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
    return L

####################################################
## assume we want to prove e(P1, Q1) = e(P2, Q2), namely e(P1, Q1) * e(P2, Q2) = 1
## fixed point Q \in G2, public known to verifier
Q1 = g2.scalar_mul(1)
Q2 = g2.scalar_mul(3)

## point P sent from prover
P1 = g1.scalar_mul(3)
P2 = g1.scalar_mul(1)

## indexer (oracle) for fixed point Q
e = 6 * x + 2
print('parameter x = {}\n'.format(x))
L1 = line_function(Q1, e)
L2 = line_function(Q2, e)
print('\n [Oracle] line function for Q1, and Q2 are both precomputated. \n')