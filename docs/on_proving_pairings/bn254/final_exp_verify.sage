import os
import copy

cur_path = os.getcwd()
load('{}/utils.sage'.format(cur_path))

## refer from https://eprint.iacr.org/2009/457.pdf
def tonelli_shanks3(c, s, t, a, k, F):
    r = a^t
    ## compute cubic root of (a^t)^-1, h
    h, cc, c = 1, c ** (3^(s - 1)), 1 / c
    for i in range(1, s):
        delta = s - i - 1
        if delta < 0:
            d = r ** ((3^s * t) // 3^(-delta))
        else:
            d = r ** 3^delta
        if d == cc:
            h, r = h * c, r * c^3
        elif d == cc^2:
            h, r = h * c^2, r * (c^3)^2
        c = c^3
    ## recover cubic root of a
    r = a^k * h
    if t == 3 * k + 1:
        r = 1 / r
        
    assert(r ** 3 == a)
    return r

def rth_root(r, r_co, f, F):
    assert(f ** r_co == F(1))
    g = gcd(r, r_co)
    assert(g[0] == 1)
    
    ## r' r = 1 mod h
    ## when sign > 0, v * r - u * M = 1
    ## when sign < 0, u * M - v * r = 1
    if g[2]:
        r_inv = r_co - g[1]
    else:
        r_inv = g[1] 
    assert((r * r_inv) % r_co == 1)
    
    root = f ** r_inv
    assert(root ** r == f)
    
    return root

def test_final_exp_witness():
    ## constant parameters
    # x = 4965661367192848881
    p, r, h, lamb, m = px(x), rx(x), hx(x), lambdax(x), mx(x)
    d = gcd(m, h)[0]
    mm = m // d

    assert(r * h == p^12 - 1)
    assert(m * r == lamb)
    assert(d * mm == m)
    assert(d == 3)

    ########################################################## tower field
    ## Fp2 = Fp[u] / X^2 - alpha, where alpha = -1
    ## Fp6 = Fp2[v] /X^3 - beta, where beta = alpha + 9
    ## Fp12 = Fp6[w] / X^2 - gamma, where gamma = v
    ## full extension field, Fp12 = Fp[w] / X^12 - 18 * X^6 + 82 
    Fp = GF(p)
    X = Fp['X'].gen()
    # pol12 = X^12 - 6 * X^6 + 10 
    pol12 = X^12 - 18 * X^6 + 82
    Fp12 = Fp.extension(pol12, 'y')
    ###########################################################

    s = 3
    t = (p^12 - 1) // 3^s
    k = (t + 1) // 3

    ################################################### step 1: choose f which is a r-th residue, but not a 3-th residue
    a = Fp12.random_element()
    f = a ** r
    while f ** (3^(s - 1) * t) == Fp12(1):
        a = Fp12.random_element()
        f = a ** r
    assert(f ** h == Fp12(1))
    print('[## 1]: Choosing a proper f representing miller loop result done!\n')
    ###################################################

    ################################################### step 2: given f, resolve r-th root of f
    f1 = rth_root(r, h, f, Fp12)
    assert(f1 ** r, f)
    print('[##2]: Resolving r-th root of f done!\n')
    ##################################################

    ################################################## step 3: given f1, resolve m'-th root of f1
    f2 = rth_root(mm, r * h, f1, Fp12)
    assert(f2 ** mm, f1)
    print('[##3]: Resolving m\'-th root of f1 done!\n')
    ##################################################

    ################################################## step 4: give f2, resolve d-th(d = 3) root of f2
    ## failed since f2 is not d-th residue
    g = gcd(d, mm * r * h)
    assert(g[0] == d)
    print('[##4]: Resolving d-th root of f2 failed, f need to be scaled... \n')
    ##################################################

    ################################################# step 5: choose a proper scalar wi for f
    w, z = Fp12(1), Fp12(1)
    while w == Fp12(1):
        ## choose z which is 3-th non-residue
        legendre = Fp12(1)
        while legendre == Fp12(1):
            z = Fp12.random_element()
            legendre = z ** (3^(s - 1) * t)
        ## obtain w which is t-th power of z
        w = z ** t
    ## make sure 27-th root w, is 3-th non-residue and r-th residue
    assert(w ** (3^(s - 1) * t) != Fp12(1))
    assert(w ** h == Fp12(1))
    ## just two option, w and w^2, since w^3 must be cubic residue, leading f*w^3 must not be cubic residue
    wi = w
    if (f * w) ** (3^(s - 1) * t) != Fp12(1):
        assert((f * w^2) ** (3^(s - 1) * t) == Fp12(1))
        wi = w^2
    assert(wi ** h == Fp12(1))
    f1 = f * wi
    print('[##5]: Choosing proper z(3-th non-residue), w(t-th power of z) and scalar wi(w or w^2) for f done!\n')
    #################################################

    ################################################# step 6: repeat above steps of r-th root, resolve r-th, mm-th, and d-th root
    assert(f1 ** h == Fp12(1))
    f2 = rth_root(r, h, f1, Fp12)
    assert(f2 ** r == f1)
    f3 = rth_root(mm, r * h, f2, Fp12)
    assert(f3 ** (mm * r) == f1)
    c = tonelli_shanks3(w, s, t, f3, k, Fp12)
    assert(d * mm * r == lamb)
    print('[##6]: Resolving r-th root of scaled f, which is f1, done, output is f2; \
          then resolving m\'-th root of f2 done, output is f3; \
          resolving d-th root of f3 done, output is c')
    #################################################

    print('\n=======================================================================================\n')
    print('witness for verification of final exponentiation: \n\n c = {} \n\n wi = {}\n'.format(c, wi))
    print('\n=======================================================================================\n')

    ################################################# verification of final exponentiation
    assert(c ** lamb == f * wi)
    #################################################
    print('verification for final exponentiation done!')

test_final_exp_witness()