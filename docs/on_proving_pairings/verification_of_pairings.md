## prelimilaries

### equivalent class

If two pairing result are equal, their miller loop results must lie on the same **equivalent class**.

$$
\def\arraystretch{1.5}
   \begin{array}{c:c:c:c}
    \hline
   \color{blue}{1} & \color{red}{f^3} & \color{red}{f^6} & \color{red}{f^9} \\ \hdashline
   \color{green}{f^4} & f^1 & f^7 & f^{10} \\ \hdashline
   \color{green}{f^8} & f^2 & f^5 & f^{11} \\ \hline
\end{array}
$$

For instance, $f^1$ and $f^7$ lie on the $f^4$-represented **equivalent class**.

<br />

Assume $r = 3$, set of multiplicative group is $F^{\times} = \{f^{i}\}, i \in [0, 12]$.

Then the $r$-torsion subgroup of $F^{\times}$ is $F^{\times}[r] = \{1, f^4, f^8\}$, and $r F^{\times} = \{1, f^3, f^6, f^9\}$, finally the quotient group:
$$
F^{\times} / r F^{\times} = \{\{1, f^3, f^6, f^9\}, \{f^4, f^1, f^7, f^{10}\}, \{f^8, f^2, f^5, f^{11}\}\}
$$

<br />

If exist two **Miller Loop** results $\mu_{r_a}$ and $\mu_{r_b}$, they have the same equivalence class, namely they are equal after **Final Exponentiation**:
$$
\mu_{r_a}^{\frac{p^k - 1}{r}} \equiv \mu_{r_b}^{\frac{p^k - 1}{r}}
$$

According to the property of coset, they must have some specific relationship:
$$
\mu_{r_a} \cdot c^r \equiv \mu_{r_b}
$$
where $c \in F^{\times}$, then $c^r$ must lies on on set $r F^{\times}$.

<br />

### Tonelli-Shanks for cubic root

referring https://eprint.iacr.org/2009/457.pdf

<br />

## verification of final exponentiation

prove:
$$
e(P_1, Q_1) = e(P_2, Q_2)
$$
where the right operators $Q_1$ and $Q_2$ are usually fixed, especially in **KZG** commitment.

<br />

take it as a product of two pairings:
$$
e(P_1, Q_1) \cdot e(-P_2, Q_2) \equiv 1
$$

if there are multiple pairings, it would be like this:
$$
\prod_{i = 0}^n e(P_i, Q_i) \equiv 1
$$

<br />

### $f$ is $r$-th residue

Assume $n = 2$, and $e(P_1, Q_1) = \mu_{r_a}, e(-P_2, Q_2) = \mu_{r_b}$, and miller loop results of these two pairings are $f_a$ and $f_b$ respectively:
$$
\begin{aligned}
f_a^h = \mu_{r_a} \\
f_b^h = \mu_{r_b} \\
\end{aligned}
$$

so we need to prove:
$$
(f_a \cdot f_b)^h \equiv 1
$$
where $h = \frac{p^k - 1}{r}$, namely $(f_a \cdot f_b)$ is a $r$-th residue, let $f_a \cdot f_b = \xi_r$.

<br />

### $f$ is $m'$-th residue

<br />

### $f$ is not $d$-th residue

<br />

### scaled $f$ is $d$-th residue

<br />

### witness for verification of final exponentiation

```python=
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

x = ZZ['x'].gen()

## parameter polynomials
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

## constant parameters
x = 4965661367192848881
p, r, h, lamb, m = px(x), rx(x), hx(x), lambdax(x), mx(x)
d = gcd(m, h)[0]
mm = m // d

########################################################## tower field
## Fp2 = Fp[u] / X^2 - alpha, where alpha = -1
## Fp6 = Fp2[v] /X^3 - beta, where beta = alpha + 9
## Fp12 = Fp6[w] / X^2 - gamma, where gamma = v
## full extension field, Fp12 = Fp[w] / X^12 - 18 X^6 + 82
Fp = GF(p)
X = Fp['X'].gen()
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

```

<br />

### verification of final exponentiation

```python=
################################################# verification of final exponentiation
assert(c ** lamb == f * wi)
#################################################
print('verification for final exponentiation done!')
```

<br />

## Verification of Miller Loop

<br />

## Verification of Multiplication of $F_{p^k}$

<br />