import os

cur_path = os.getcwd()
load('{}/miller_loop_verify.sage'.format(cur_path))
load('{}/final_exp_verify.sage'.format(cur_path))

################################
# verify c^lambda = f * wi, namely c_inv^lambda * f * wi = 1
################################
def verify_pairings(eval_points, lines, e, c, c_inv, wi):
    assert(len(eval_points) == len(lines))
    assert(c.mul(c_inv) == Fp12.ONE())

    lc = 0
    # f = Fp12.ONE()
    f = c_inv
    naf_digits = list(reversed(to_naf(e)))[1:]
    ## double-add part, 6x + 2 
    for i, digit in enumerate(naf_digits):
        f = f.square()
        ## update c^lambda
        if digit^2 == 1:
            f = f.mul(c_inv) if digit == 1 else f.mul(c)

        for j, (P, L) in enumerate(zip(eval_points, lines)):
            alpha, bias = L[lc]
            le = line_evaluation(alpha, bias, P)
            f = mul_line_base(f, le[0], le[1], le[2])

            if digit^2 == 1:
                alpha, bias = L[lc + 1]
                le = line_evaluation(alpha, bias, P)
                f = mul_line_base(f, le[0], le[1], le[2])
        lc = (lc + 1) if digit == 0 else (lc + 2)
    
    ## update c^lambda
    f = f.mul(c_inv.frobenius()).mul(c.frobenius_p2()).mul(c_inv.frobenius_p3())
    ## update the scalar
    f = f.mul(wi)

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
    lc = lc + 3

    assert(lc == len(lines[0]))
    assert(f == Fp12.ONE())

    return f

##################################
# obtain c and wi, satisfying c^lambda = f * wi
##################################
def compute_final_exp_witness(f):
    ## constant parameters
    p, r, h, lamb, m = px(x), rx(x), hx(x), lambdax(x), mx(x)
    d = gcd(m, h)[0]
    mm = m // d
    s = 3
    t = (p^12 - 1) // 3^s
    k = (t + 1) // 3
    assert(r * h == p^12 - 1)
    assert(m * r == lamb)
    assert(d * mm == m)
    assert(d == 3)
    f_copy = f

    ########################################################## tower field 
    ## Fp2 = Fp[u] / X^2 - alpha, where alpha = -1
    ## Fp6 = Fp2[v] /X^3 - beta, where beta = alpha + 9
    ## Fp12 = Fp6[w] / X^2 - gamma, where gamma = v
    ## full extension field, Fp12 = Fp[w] / X^12 - 18 * X^6 + 82 
    ## we do not implement the random traits for Fp12, so we use Galois Field of sagemath instead
    Fp_local = GF(p)
    X = Fp_local['X'].gen()
    pol12 = X^12 - 18 * X^6 + 82 
    Fp12_local = Fp_local.extension(pol12, 'y')
    ###########################################################
    # assert(final_exp(f) == Fp12.ONE())
    # print(f.exp(3^(s - 1) * t))
    # assert(f.exp(3^(s - 1) * t) != Fp12.ONE())

    f = Fp12_local(
        [
        (f.y.z.y - Fp(9).mul(f.y.z.x)).value(),
        (f.x.z.y - Fp(9).mul(f.x.z.x)).value(),
        (f.y.y.y - Fp(9).mul(f.y.y.x)).value(),
        (f.x.y.y - Fp(9).mul(f.x.y.x)).value(),
        (f.y.x.y - Fp(9).mul(f.y.x.x)).value(),
        (f.x.x.y - Fp(9).mul(f.x.x.x)).value(),
        (f.y.z.x).value(),
        (f.x.z.x).value(),
        (f.y.y.x).value(),
        (f.x.y.x).value(),
        (f.y.x.x).value(),
        (f.x.x.x).value(),
        ]
    )

    ## make sure f is r-th residue, but not d-th (or 3-th) residue
    assert(f ** h == Fp12_local(1))
    assert(f ** (3^(s - 1) * t) != Fp12_local(1))

    g = gcd(d, mm * r * h)
    assert(g[0] == d)

    ########################################################### choose a proper scalar wi for f
    w, z = Fp12_local(1), Fp12_local(1)
    while w == Fp12_local(1):
        ## choose z which is 3-th non-residue
        legendre = Fp12_local(1)
        while legendre == Fp12_local(1):
            z = Fp12_local.random_element()
            legendre = z ** (3^(s - 1) * t)
        ## obtain w which is t-th power of z
        w = z ** t
    ## make sure 27-th root w, is 3-th non-residue and r-th residue
    assert(w ** (3^(s - 1) * t) != Fp12_local(1))
    assert(w ** h == Fp12_local(1))
    ## just two option, w and w^2, since w^3 must be cubic residue, leading f*w^3 must not be cubic residue
    wi = w
    if (f * w) ** (3^(s - 1) * t) != Fp12_local(1):
        assert((f * w^2) ** (3^(s - 1) * t) == Fp12_local(1))
        wi = w^2
    assert(wi ** h == Fp12_local(1))
    f1 = f * wi
    ###########################################################

    assert(f1 ** h == Fp12_local(1))
    f2 = rth_root(r, h, f1, Fp12_local)
    assert(f2 ** r == f1)
    f3 = rth_root(mm, r * h, f2, Fp12_local)
    assert(f3 ** (mm * r) == f1)
    c = tonelli_shanks3(w, s, t, f3, k, Fp12_local)
    assert(d * mm * r == lamb)

    ## make sure the witness c and wi is what we want
    assert(c ** lamb == f * wi)

    ## convert into self contained Fp12 object
    c_list = c.polynomial().list()
    c = Fp12(
        Fp6(
            Fp2(int(c_list[11]), int(9 * c_list[11] + c_list[5])),
            Fp2(int(c_list[9]), int(9 * c_list[9] + c_list[3])),
            Fp2(int(c_list[7]), int(9 * c_list[7] + c_list[1])),
        ),
        Fp6(
            Fp2(int(c_list[10]), int(9 * c_list[10] + c_list[4])),
            Fp2(int(c_list[8]) , int(9 * c_list[8] + c_list[2])),
            Fp2(int(c_list[6]) , int(9 * c_list[6] + c_list[0])),
        ),
    )

    wi_list = wi.polynomial().list()
    wi_list = wi_list + [Fp_local(0)] * (12 - len(wi_list)) if len(wi_list) < 12 else wi_list
    wi = Fp12(
        Fp6(
            Fp2(int(wi_list[11]), int(9 * wi_list[11] + wi_list[5])),
            Fp2(int( wi_list[9]), int(9 * wi_list[9] + wi_list[3])),
            Fp2(int( wi_list[7]), int(9 * wi_list[7] + wi_list[1])),
        ),
        Fp6(
            Fp2(int(wi_list[10]), int(9 * wi_list[10] + wi_list[4])),
            Fp2(int( wi_list[8]), int(9 * wi_list[8] + wi_list[2])),
            Fp2(int( wi_list[6]), int(9 * wi_list[6] + wi_list[0])),
        ),
    )

    assert(c.exp(lamb) == f_copy.mul(wi))

    return c, wi

####################################################
## assume we want to prove e(P1, Q1) = e(P2, Q2), namely e(P1, Q1) * e(P2, -Q2) = 1

## fixed point Q in G2, public known to verifier
Q1 = g2.scalar_mul(1).force_affine()
Q2 = g2.scalar_mul(3).force_affine()

## point P being evaluated sent from prover
P1 = g1.scalar_mul(3).force_affine()
P2 = g1.scalar_mul(1).force_affine()

## witness for miller loop sent from prover
f1, _ = miller(Q1, P1)
f2, _ = miller(Q2.negate(), P2)
c, wi = compute_final_exp_witness(f1.mul(f2))
c_inv = c.inverse()
print('[## Prover-Side] Witness for miller loop has be obtained.\n\n')

## verify pairing with provided witness c, wi and precomputed lines  
e = 6 * x + 2
lamb = lambdax(x)
L1 = line_function(Q1.force_affine(), e, lamb)
L2 = line_function(Q2.negate().force_affine(), e, lamb)
verify_pairings([P1, P2], [L1, L2], e, c, c_inv, wi)
print('[## Verifier-Side] Verify pairings with provided mill loop witness and procomputed lines done.\n\n')
