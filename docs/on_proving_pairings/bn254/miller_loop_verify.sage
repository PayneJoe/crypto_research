import os

cur_path = os.getcwd()
load('{}/line_precomputation.sage'.format(cur_path))

## we use affine coordinate to verify the line evaluation
## (-b) + y_P * w^3 + (-alpha * x_P) * w^2 where w \in Fp12
def line_evaluation(alpha, bias, P):
    P = P.force_affine()
    assert(P.z == P.one_element())
    assert(type(P) == G1)
    assert(type(alpha) == Fp2)
    assert(type(bias) == Fp2)

    return (
        bias.negative_of(),
        alpha.negative_of().mul_scalar(P.x),
        Fp2(Fp.ZERO(), P.y)
    )

def miller_loop_verify(Q, P, L, e, f_baseline):
    f = Fp12.ONE()
    f_check = Fp12.ONE()
    lc = 0
    naf_digits = list(reversed(to_naf(e)))[1:]
    T = Q
    for i, digit in enumerate(naf_digits):
        alpha, bias = L[lc]
        ####################### we only need to verify the line evaluation, and accumulating the line evaluation
        le = line_evaluation(alpha, bias, P)
        f = f.square()
        f = mul_line(f, le[0], le[1], le[2])
        ######################

        ######################  the baseline use the projective coordinate to do the line evaluation(avoiding inversion), verification use affine coordinate
        aux = T.z.square().mul(T.z).mul(T.y).mul_scalar(Fp(4))
        T = T.double()
        f_check = f_check.square()
        f_check = mul_line(f_check, le[0].mul(aux), le[1].mul(aux), le[2].mul(aux))
        assert(f_check == f_baseline[lc])
        #####################

        if digit^2 == 1:
            lc = lc + 1
            alpha, bias = L[lc]
            le = line_evaluation(alpha, bias, P)
            f = mul_line(f, le[0], le[1], le[2])

            ####### just for testation
            aux = T.z.square().mul(Q.x).sub(T.x).mul(T.z).mul_scalar(Fp(4))
            T = T.add(Q) if digit == 1 else T.add(Q.negate())
            f_check = mul_line(f_check, le[0].mul(aux), le[1].mul(aux), le[2].mul(aux))
            assert(f_check == f_baseline[lc])
            ########

        lc = lc + 1
    assert(lc == len(L))
    return f

############################################ miller loop for two pairings, e(Q1, P1) and e(Q2, P2), with precomputed line functions
# f1 = miller_loop_verify(P1, L1, e) 
# assert(f1 == miller(Q1, P1)[1])
# f2 = miller_loop_verify(P2, L2, e)

# print('\n e(Q1, P1) = {}, e(Q2, P2) = {} \n'.format(f1, f2))
# print('Miller loop verification with precomputated line functions done.')

_, f1_check = miller(Q1, P1)
miller_loop_verify(Q1, P1, L1, e, f1_check) 