import os

cur_path = os.getcwd()
load('{}/line_precomputation.sage'.format(cur_path))

## (-b) + y_P * w^3 + (-alpha * x_P) * w^2 where w \in Fp12
def line_evaluation(alpha, bias, P):
    P = P.force_affine()
    assert(P.z == P.one_element())
    assert(type(P) == G1)
    assert(type(alpha) == Fp2)
    assert(type(bias) == Fp2)

    return Fp6(
        alpha.negative_of().mul_scalar(P.x),
        Fp2(Fp.ZERO(), P.y),
        bias.negative_of()
    )

def miller_loop_verify(P, L, e):
    f = Fp12.ONE()
    lc = 0
    naf_digits = list(reversed(to_naf(e)))[1:]
    for i, digit in enumerate(naf_digits):
        alpha, bias = L[lc]
        le = line_evaluation(alpha, bias, P)
        f = f.square()
        mul_line(f, le.x, le.y, le.z)
        if digit^2 == 1:
            lc = lc + 1
            alpha, bias = L[lc]
            le = line_evaluation(alpha, bias, P)
            mul_line(f, le.x, le.y, le.z)
        lc = lc + 1
    assert(lc == len(L))
    return f

############################################ miller loop for two pairings, e(Q1, P1) and e(Q2, P2), with precomputed line functions
f1 = miller_loop_verify(P1, L1, e) 
f2 = miller_loop_verify(P2, L2, e)

print('\n e(Q1, P1) = {}, e(Q2, P2) = {} \n'.format(f1, f2))
print('Miller loop verification with precomputated line functions done.')