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

## make sure the line evaluation result on verifier-side is consistant with the one on prover-side
def miller_loop_check(P, L, e, f_check):
    f = Fp12.ONE()
    lc = 0
    naf_digits = list(reversed(to_naf(e)))[1:]
    for i, digit in enumerate(naf_digits):
        alpha, bias = L[lc]
        ####################### we only need to verify the line evaluation, and accumulating the line evaluation
        le = line_evaluation(alpha, bias, P)
        f = f.square()
        f = mul_line(f, le[0], le[1], le[2])
        assert(f == f1_check[lc])
        ######################

        ######################  the baseline use the projective coordinate to do the line evaluation(avoiding inversion), verification use affine coordinate
        # aux = T.z.square().mul(T.z).mul(T.y).mul_scalar(Fp(4))
        # T = T.double()
        # f_check = f_check.square()
        # f_check = mul_line(f_check, le[0].mul(aux), le[1].mul(aux), le[2].mul(aux))
        # assert(f_check == f_baseline[lc])
        #####################

        if digit^2 == 1:
            lc = lc + 1
            alpha, bias = L[lc]
            le = line_evaluation(alpha, bias, P)
            f = mul_line(f, le[0], le[1], le[2])
            assert(f == f1_check[lc])

            ####### just for checking, ignored while verification
            # aux = T.z.square().mul(Q.x).sub(T.x).mul(T.z).mul_scalar(Fp(4))
            # T = T.add(Q) if digit == 1 else T.add(Q.negate())
            # f_check = mul_line(f_check, le[0].mul(aux), le[1].mul(aux), le[2].mul(aux))
            # assert(f_check == f_baseline[lc])
            ########

        lc = lc + 1
    assert(lc == len(L) - 3)
    return f

## miller loop for multiple pairings, especially for \prod_i^n e(Qi, Pi) = 1
def multi_miller_loop(eval_points, lines, e):
    assert(len(eval_points) == len(lines))

    f_list = []
    # f_check = [Fp12.ONE() for _ in range(len(lines))]

    lc = 0
    f = Fp12.ONE()
    naf_digits = list(reversed(to_naf(e)))[1:]
    ## double-add part, 6x + 2 
    for i, digit in enumerate(naf_digits):
        f = f.square()
        # f_check = [e.square() for e in f_check]

        for j, (P, L) in enumerate(zip(eval_points, lines)):
            alpha, bias = L[lc]
            le = line_evaluation(alpha, bias, P)
            f = mul_line_base(f, le[0], le[1], le[2])
            # f = mul_line(f, le[0], le[1], le[2])
            f_list.append(f)

            # f_check[j] = mul_line_base(f_check[j], le[0], le[1], le[2])

            if digit^2 == 1:
                alpha, bias = L[lc + 1]
                le = line_evaluation(alpha, bias, P)
                f = mul_line_base(f, le[0], le[1], le[2])
                # f = mul_line(f, le[0], le[1], le[2])
                f_list.append(f)

                # f_check[j] = mul_line_base(f_check[j], le[0], le[1], le[2])

        lc = (lc + 1) if digit == 0 else (lc + 2)
    
    ## frobenius map part, p - p^2 + p^3
    for j, (P, L) in enumerate(zip(eval_points, lines)):
        for k in range(3):
            alpha, bias = L[lc + k]
            if k == 2:
                    eval = Fp12(
                        Fp6.ZERO(),
                        Fp6(Fp2.ZERO(), bias.negative_of(), Fp2(Fp.ZERO(), P.x))
                    )
                    f = f.mul(eval)
            else:
                le = line_evaluation(alpha, bias, P)
                f = mul_line_base(f, le[0], le[1], le[2])
            # f = mul_line(f, le[0], le[1], le[2])
            f_list.append(f)

            # f_check[j] = mul_line_base(f_check[j], le[0], le[1], le[2])
    lc = lc + 3

    # assert(f_check[0].mul(f_check[1]) == f)

    assert(lc == len(lines[0]))
    return f, f_list

def test_precomputed_multi_miller_loop():
    ############################### Testation for miller loop,  with precomputed line functions
    # _, f1_check = miller(Q1, P1)
    # miller_loop_check(P1, L1, e, f1_check) 
    # print('[Test] miller loop with precomputed lines.\n\n')

    ####################################################
    ## assume we want to prove e(P1, Q1) = e(P2, Q2), namely e(P1, Q1) * e(P2, -Q2) = 1
    ## fixed point Q in G2, public known to verifier
    Q1 = g2.scalar_mul(1).force_affine()
    Q2 = g2.scalar_mul(3).force_affine()

    ## point P sent from prover
    P1 = g1.scalar_mul(3).force_affine()
    P2 = g1.scalar_mul(1).force_affine()

    e = 6 * x + 2
    lamb = lambdax(x)
    print('parameter x = {}\n'.format(x))
    L1 = line_function(Q1.force_affine(), e, lamb)
    L2 = line_function(Q2.negate().force_affine(), e, lamb)
    
    print('==================== Test Module of Miller Loop ================= \n')
    f1_a, f1_a_list = miller(Q1, P1)
    f1_b, f1_b_list = multi_miller_loop([P1], [L1], e)
    assert(len(f1_a_list) == len(f1_b_list))
    t1 = (f1_a == f1_b)
    print('[Test] multi_miller_loop(P1, L1, e)? {}\n'.format(t1))
    
    f2_a, f2_a_list = miller(Q2.negate().force_affine(), P2)
    f2_b, f2_b_list = multi_miller_loop([P2], [L2], e)
    assert(len(f2_a_list) == len(f2_b_list))
    assert(f2_a == f2_b)
    t2 = (f1_a.mul(f2_a) == f1_b.mul(f2_b))
    print('[Test] multi_miller_loop(P2, L2, e)? {}\n'.format(t2))
    
    f, _ = multi_miller_loop([P1, P2], [L1, L2], e)
    t3 = (f == f1_a.mul(f2_a))
    print('[Test] multi_miller_loop([P1, P2], [L1, L2], e)? {}\n'.format(t3))

    t4 = (final_exp(f) == Fp12.ONE())
    print('[Test] multi_miller_loop([P1, P2], [L1, L2], e)^h == 1? {}\n'.format(t4))
    print('==================== End of Test Module of Miller Loop ================= \n')
    ###############################

# test_precomputed_multi_miller_loop()