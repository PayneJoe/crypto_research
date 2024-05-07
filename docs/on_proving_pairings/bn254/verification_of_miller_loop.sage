def gcd(x, N):
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
    
    return r, u, i % 2 == 0

########################################################### parameter polynomials
x = ZZ['x'].gen()
## prime field, p
px = 36 * x^4 + 36 * x^3 + 24 * x^2 + 6 * x + 1
## largest prime facotr, r
rx = 36 * x^4 + 36 * x^3 + 18 * x^2 + 6 * x + 1
## trace of frobenius, t
tx = 6 * x^2 + 1
## cyclotomic group of 12-extension degree, phi_12
phi_12 = px^4 - px^2 + 1
## power final exponentiation, h
hx = (px^12 - 1) // rx
## optimal lambda in miller loop, lambda
lambdax = 6 * x + 2 + px - px^2 + px^3
## multiples of r, m
mx = lambdax // rx
#########################################################

## constant parameters
# x = 4965661367192848881
x = 6518589491078791937
p, r, h, lamb, m = px(x), rx(x), hx(x), lambdax(x), mx(x)
d = gcd(m, h)[0]
mm = m // d

########################################################## tower field
## Fp2 = Fp[u] / X^2 - alpha, where alpha = -1
## Fp6 = Fp2[v] /X^3 - beta, where beta = u + 3
## Fp12 = Fp6[w] / X^2 - gamma, where gamma = v
## full extension field, Fp12 = Fp[w] / X^12 - 6 * X^6 + 10

Fp = GF(p)
X = Fp['X'].gen()
alpha = Fp(-1)

pol2 = X^2 - alpha
Fp2 = Fp.extension(pol2, 'u')
u = Fp2.gen()
beta = u + 3

pol6 = X^3 - beta
Fp6 = Fp2.extension(pol6, 'v')

pol12 = X^12 - 6 * X^6 + 10
Fp12 = Fp.extension(pol12, 'y')
###########################################################

########################################################### G1
a, b = Fp(0), Fp(3)
G1 = EllipticCurve(Fp, [a, b])
assert(G1.order() % r == 0)
co_g1 = G1.order() // r
g1 = G1.random_element() * co_g1
assert(g1 * r == G1(0))
###########################################################
print('[##G1]: cofactor = {}, \n generator = {}\n'.format(co_g1, g1))

########################################################### G2
beta_t = beta
if EllipticCurve(Fp2, [a, b * beta_t]).order() % r != 0:
    # beta_t = beta ** 5
    # print('[!] beta_t = beta^5 \n')
    beta_t = 1 / beta
G2 = EllipticCurve(Fp2, [a, b * beta_t])
assert(G2.order() % r == 0)
co_g2 = G2.order() // r
g2 = G2.random_element() * co_g2
assert(g2 * r == G2(0))
###########################################################
print('[##G2]: cofactor = {}, \n generator = {}\n'.format(co_g2, g2))