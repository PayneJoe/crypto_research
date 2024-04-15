tag = "bls12_scalar"

def inv_mod(x, N):
    assert(x < N)
    
    (A, B) = (N, x)
    (Ua, Ub) = (0, 1)
    (Va, Vb) = (1, 0)
    i = 0
    while B:
        q = A // B
        (A, B) = (B, A % B)
        (Ua, Ub) = (Ub, Ua + q * Ub)
        (Va, Vb) = (Vb, Va + q * Vb)
        i += 1
    r, u, v = A, Ua, Va
    
    assert(r == 1)
    
    return u, v, i % 2 == 0

def odd_factor(M):
    t = M - 1
    e = 0
    while t % 2 == 0:
        t = t >> 1
        e += 1
    return e, t

def sample_nqr(M):
    Fp = GF(M)
    t = Fp(2)
    while True:
        if t ** ((M - 1) // 2) == Fp(-1):
            break
        t += 1
    #print(Fp(5).sqrt())
    return t

NUM_LIMBS = {
    'bls12_base': 6,
    'bls12_scalar': 4,
}
MODULUS = {
    'bls12_base': 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787,
    'bls12_scalar': 52435875175126190479447740508185965837690552500527637822603658699938581184513,
}

W = 2 ** 64
r = W ** NUM_LIMBS[tag]

## params for Montgomery Reduction
R = r % MODULUS[tag]
R2 = (r * r) % MODULUS[tag]
R3 = (r * r * r) % MODULUS[tag]
u, v, sign = inv_mod(MODULUS[tag], r)
M0 = None
## when sign > 0, v * r - u * M = 1, M0 = ((-M^{-1}) % r) % W = u % W
## when sign < 0, u * M - v * r = 1, M0 = ((-M^{-1}) % r) % W = (-u) % W
if sign:
    assert(u * MODULUS[tag] == (v * r - 1))
    M0 = u % W
else:
    assert(u * MODULUS[tag] == (v * r + 1))
    M0 = (r - u) % W

## params for Sqrt 
## MODULUS = 2^E * RODD + 1
E, RODD = odd_factor(MODULUS[tag])
N = sample_nqr(MODULUS[tag])
print('\n M = {}, \n r = {}, \n R = {}, \n R2 = {}, \n R3 = {}, \n M0 = {}, \n E = {}, \n RODD = {}, \n N = {} \n'.format(
    MODULUS[tag], r, R, R2, R3, M0, E, RODD, N
))
print(N.sqrt())