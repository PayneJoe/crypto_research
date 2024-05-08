import os

cur_path = os.getcwd()
load('{}/curves.sage'.format(cur_path))

##################################################################### Optimal Ate Pairing
def line_func_add(r, p, q, r2):
    assert type(r) == G2
    assert type(p) == G2
    assert type(q) == G1
    assert type(r2) == Fp2

    r_t = r.z.square()
    B = p.x * r_t
    D = p.y + r.z
    D = D.square()
    D -= r2
    D -= r_t
    D *= r_t

    H = B - r.x
    I = H.square()

    E = I.double().double()

    J = H * E
    L1 = D - r.y
    L1 -= r.y

    V = r.x * E

    r_x = L1.square()
    r_x -= J
    r_x -= V.double()

    r_z = r.z + H
    r_z = r_z.square()
    r_z -= r_t
    r_z -= I

    t = V - r_x
    t *= L1
    t2 = r.y * J
    t2 = t2.double()
    r_y = t - t2

    r_out = G2(r_x, r_y, r_z)

    t = p.y + r_z
    t = t.square()
    t = t - r2
    t = t - (r_z.square())

    t2 = L1 * p.x
    t2 = t2.double()
    a = t2 - t

    c = r_z.mul_scalar(q.y).double()

    b = L1.negative_of()
    b = b.mul_scalar(q.x).double()

    return (a, b, c, r_out)

def line_func_double(r, q):
    assert type(r) == G2
    assert type(q) == G1

    # cache this?
    r_t = r.z.square()

    A = r.x.square()
    B = r.y.square()
    C = B.square()

    D = r.x + B
    D = D.square()
    D -= A
    D -= C
    D = D.double()

    E = A.double() + A
    F = E.square()

    C8 = C.double().double().double() # C*8

    r_x = F - D.double()
    r_y = E * (D - r_x) - C8

    # (y+z)*(y+z) - (y*y) - (z*z) = 2*y*z
    r_z = (r.y + r.z).square() - B - r_t

    assert r_z == r.y*r.z.double()

    r_out = G2(r_x, r_y,r_z)
    #assert r_out.is_on_curve()

    a = r.x + E
    a = a.square()
    a -= (A + F + B.double().double())

    t = E * r_t
    t = t.double()
    b = t.negative_of()
    b = b.mul_scalar(q.x)

    c = r_z * r_t
    c = c.double().mul_scalar(q.y)

    return (a,b,c,r_out)

def mul_line(r, a, b, c):
    assert type(r) == Fp12
    assert type(a) == Fp2
    assert type(b) == Fp2
    assert type(c) == Fp2

    # See function fp12e_mul_line in dclxvi

    t1 = Fp6(Fp2.ZERO(), a, b)
    t2 = Fp6(Fp2.ZERO(), a, b.add(c))

    t1 = t1.mul(r.x)
    t3 = r.y.mul_scalar(c)

    x = r.x.add(r.y)
    y = t3
    x = x.mul(t2)
    x = x.sub(t1)
    x = x.sub(y)
    y = y.add(t1.mul_tau())

    return Fp12(x, y) 

def miller(q, p):

    import copy

    assert type(q) == G2
    assert type(p) == G1

    Q = copy.deepcopy(q)
    Q = Q.force_affine()

    P = copy.deepcopy(p)
    P = P.force_affine()

    mQ = copy.deepcopy(Q)
    mQ = mQ.negate()

    f = Fp12(Fp6.ZERO(), Fp6.ONE())
    T = Q

    Qp = Q.y.square()
    
    # 6x + 2 in NAF
    naf_6xp2 = list(reversed(to_naf(6 * x + 2)))[1:]

    f_list = []
    for i, naf_i in enumerate(naf_6xp2):
        # Skip on first iteration?
        f = f.square()

        a, b, c, T = line_func_double(T, P)
        f = mul_line(f, a, b, c)
        f_list.append(f)

        if naf_i == 1:
            a, b, c, T = line_func_add(T, Q, P, Qp)
            f = mul_line(f, a, b, c)
            f_list.append(f)
        elif naf_i == -1:
            a, b, c, T = line_func_add(T, mQ, P, Qp)
            f = mul_line(f, a, b, c)
            f_list.append(f)

    # Q1 = pi(Q)
    Q1 = G2(
        Q.x.conjugate_of().mul(Fp12.beta_pi_1[1]),
        Q.y.conjugate_of().mul(Fp12.beta_pi_1[2]),
        Fp2.ONE())

    # Q2 = pi2(Q)
    Q2 = G2(
        Q.x.mul_scalar(Fp12.beta_pi_2[1].y),
        Q.y,
        Fp2.ONE())

    Qp = Q1.y.square()
    a, b, c, T = line_func_add(T, Q1, P, Qp)
    mul_line(f, a, b, c)

    Qp = Q2.y.square()
    a, b, c, T = line_func_add(T, Q2, P, Qp)
    mul_line(f, a, b, c)

    return f, f_list

def final_exp(inp):
    assert type(inp) == Fp12

    # Algorithm 31 from https://eprint.iacr.org/2010/354.pdf

    t1 = inp.conjugate_of()
    inv = inp.inverse()

    t1 = t1.mul(inv)
    # Now t1 = inp^(p**6-1)

    t2 = t1.frobenius_p2()
    t1 = t1.mul(t2)

    fp1 = t1.frobenius()
    fp2 = t1.frobenius_p2()
    fp3 = fp2.frobenius()

    fu1 = t1.exp(x)
    fu2 = fu1.exp(x)
    fu3 = fu2.exp(x)

    y3 = fu1.frobenius()
    fu2p = fu2.frobenius()
    fu3p = fu3.frobenius()
    y2 = fu2.frobenius_p2()

    y0 = fp1.mul(fp2)
    y0 = y0.mul(fp3)

    y1 = t1.conjugate_of()
    y5 = fu2.conjugate_of()
    y3 = y3.conjugate_of()
    y4 = fu1.mul(fu2p)
    y4 = y4.conjugate_of()

    y6 = fu3.mul(fu3p)
    y6 = y6.conjugate_of()

    t0 = y6.square()
    t0 = t0.mul(y4)
    t0 = t0.mul(y5)

    t1 = y3.mul(y5)
    t1 = t1.mul(t0)
    t0 = t0.mul(y2)
    t1 = t1.square()
    t1 = t1.mul(t0)
    t1 = t1.square()
    t0 = t1.mul(y1)
    t1 = t1.mul(y0)
    t0 = t0.square()
    t0 = t0.mul(t1)

    return t0

def optimal_ate(a, b):
    assert type(a) == G2
    assert type(b) == G1

    e, f_m = miller(a, b)
    mu = final_exp(e)

    if a.is_infinite() or b.is_infinite():
        return Fp12(Fp6.ZERO(), Fp6.ONE())

    return mu

################################################################# Testation
print('\n================ Test Module of Optimal Ate =====================\n')
# Q, P = g2.scalar_mul(3), g1.scalar_mul(4)
# mu_r = optimal_ate(Q, P)
# t1 = mu_r.exp(rx(x)) == Fp12.ONE()
# print('[Test] optimal_ate(Q, P) ** r == Fp12.ONE()? {}\n'.format(t1))
print('\n============== End of Test Module of Optimal Ate ================\n')