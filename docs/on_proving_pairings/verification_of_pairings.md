This note focus on paper [On Proving Pairings](https://eprint.iacr.org/2024/640.pdf), which is a great propsal on recusive snark on pairings (verification pairings within circuit).

<br />

The complete testation code is under [repo](https://github.com/PayneJoe/crypto_research/tree/main/docs/on_proving_pairings). Welcome touch-)

<br />

-------

# Prelimilaries

## BN254 implementation

Public parameters:
- Parameter $x$ of BN-family parameter polynomials
    
    $x = 4965661367192848881$
    
    for
    
    $$
        \begin{aligned}
        p &=  36 \cdot x^4 + 36 \cdot x^3 + 24 \cdot x^2 + 6 \cdot x + 1 \\
        r &= 36 \cdot x^4 + 36 \cdot x^3 + 18 \cdot x^2 + 6 \cdot x + 1 \\
        h &= (p \cdot x^{12} - 1) // r \\
        \lambda &= 6 \cdot x + 2 + p \cdot x - p \cdot x^2 + p \cdot x^3 \\
        \end{aligned}
    $$
- Modulus of base prime field (characteristic) with $254$-bits:

    $$
        p = 21888242871839275222246405745257275088696311157297823662689037894645226208583
    $$
    
- Embedding degree, or the degree of full extension field $F_{p^k}$:
    $$
    k = 12
    $$

<br />

- Elliptic Curve (additive group) over **base prime field** $F_p$:
    $$
        \mathbb{G}_1/E(F_p): y^2 = x^3 + 3
    $$

<br />

- Elliptic Curve (additive group) over extension field $F_{p^2}$ (**D-type twist**):
    $$
        \mathbb{G}_2/E(F_{p^2}): y^2 = x^3 + \frac{3}{\beta}
    $$
    where $\beta = u + 9 \in F_{p^2}$ and $u$ is the generator of $F_{p^2}$.
    
- Largest prime factor of $|E(F_p)|$ with 255-bits:
    $$
        r = 21888242871839275222246405745257275088548364400416034343698204186575808495617
    $$
    
- Target (multiplicative) group with order $r$ defined over $F_{p^k}$:
    $$
        \mathbb{G}_T: F_{p^k}^{\times}[r]
    $$

<br />

We have implemented python versioned BN254 (field/curve/pairing) with python self-contained bigint :

- [field arithemtic](https://github.com/PayneJoe/crypto_research/blob/main/docs/on_proving_pairings/bn254/fields.sage)
- [curve arithemtic](https://github.com/PayneJoe/crypto_research/blob/main/docs/on_proving_pairings/bn254/curves.sage)
- [optimal ate pairing](https://github.com/PayneJoe/crypto_research/blob/main/docs/on_proving_pairings/bn254/optimal_ate.sage)

More details about pairings you can also check another [note](https://hackmd.io/@70xfCGp1QViTYYJh3AMrQg/ryo55eEeC).

<br />

## Equivalent Class

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


# Verification of Final Exponentiation

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

## $f$ is $r$-th residue

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
where $h = \frac{p^k - 1}{r}$, namely $(f_a \cdot f_b)$ must be a $r$-th residue, let $f_a \cdot f_b = f = f_1^r$.

<br />

So we can easily obtain the $r$-th root of $f = f_a \cdot f_b$:
$$
f_1 = (f_a \cdot f_b)^{\frac{1}{r}}
$$
where $\frac{1}{r} =$ inv_mod($r, h$).

<br />

below is my code snippet：
```python=
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
```


<br />

## $f$ is $m'$-th residue

In **optimal ate pairing** the miller loop result is $\lambda$-th residue:
$$
(f^{\lambda})^h = \mu_r
$$
where $\lambda = 6 x + 2 + p - p^2 + p^3 = m \cdot r, h = \frac{p^{12} - 1}{r}$

<br />

Assume $gcd(m, h) = d$ ($d = 3$ especially in BN254), let $d \cdot m' = m$, then $gcd(m', r \cdot h) = 1$, we can easily obtain the $m'$-th root of $f_1$, $f_2 = f_1^{m'}$ with the same method:
$$
f_2 = f_1^{\frac{1}{m'}}
$$
where $\frac{1}{m'}$ = inv_mod($m', r \cdot h$).

<br />

## $f$ is not $d$-th residue

Unfortunately, $f_2$ is not a $d$-th residue, since $gcd(d, h) = d$ ($d = 3$ especially in BN254). Namely we can not obtain the $f_3^3 = f_2$, so that:
$$
f_3^{\lambda} = f_3^{3 \cdot m' \cdot r} = f = f_a \cdot f_b
$$

<br />

Authors of [On Proving Pairings](https://eprint.iacr.org/2024/640.pdf) proposed a creative solution for this.

<br />

## Scaled $f$ is $d$-th residue

We can multiply $f$ with some number say $w_i$, make sure that $f \cdot w_i$ is cubic residue, namely:
$$
f_3^{\lambda} = f_3^{3 \cdot m' \cdot r} = f \cdot w_i
$$

<br />

So that, we can sequentially obtain the $r$-th root of $f \cdot w_i$, $f_1$, and then the $m'$-th root of $f_1$, $f_2$, finally the $d$-th root of $f_2$, $f_3$:
$$
\begin{aligned}
f_1^r &= f \cdot w_i \\
f_2^{m'} &= f_1 \\
f_3^3 &= f_2 \\
\end{aligned}
$$

since $gcd(d, m' \cdot r \cdot h) = 3$, so we can not use the same method as we did above, so the modified **Tonelli Shanks** for cubic root comes:

[Algorith 4](https://eprint.iacr.org/2024/640.pdf) seems not such accurate, below is my implementation.

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
```

<br />

## Witness for Verification of Final Exponentiation

But how can prover obtain these two wintesses $c = f_3$, and $w_i$, so that verifier can easily check:
$$
c^{\lambda} = f \cdot w_i
$$
or 
$$
(c^{-1})^{\lambda} \cdot f \cdot w_i \equiv 1
$$
through 1 miller loop (especially for $(c^{-1})^{\lambda}$) and 2 multiplications on $F_{p^{12}}$. Much simpler than do the **Final Exponentiation** by verifier-self, right?

<br />

According to **[Lemma 2](https://eprint.iacr.org/2024/640.pdf)** Since $f$ is not cubic residue, we only need to find another non-cubic residue. Below is my code snippet:
```python=
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
```

<br />

## Verification of Final Exponentiation

```python=
################################################# verification of final exponentiation
assert(c ** lamb == f * wi)
#################################################
print('verification for final exponentiation done!')
```

<br />

Since $\lambda = 6x + 2 + p - p^2 + p^3$, actually we can easily have $c^{\lambda}$ embeded within a miller loop same as calculating the line evaluations.

<br />

------

# Verification of Miller Loop

## Precomputed Lines

Since pairing applications like **KZG**, the left point $Q$ usualy is fixed, therefore the lines, represented with a tuple $(\alpha, b)$ denotes *slope* and *bias* repectively, can be precomputed for verifier. 

```python=
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
    # bias = y1 - alpha * x1

    return alpha, bias

################################################## cache line parameters for [6x + 2 + p - p^2 + p^3]Q
def line_function(Q, e, lamb):
    ############################################## double-add part, 6x + 2
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

    ############################################# frobenius map part, p - p^2 + p^3
    # Q1 = pi(Q)
    # x = x' * beta^(2 * (p - 1) // 6)
    # y = y' * beta^(3 * (p - 1) // 6))
    pi_1_Q = G2(
        Q.x.conjugate_of().mul(Fp12.beta_pi_1[1]),
        Q.y.conjugate_of().mul(Fp12.beta_pi_1[2]),
        Fp2.ONE())
    assert(pi_1_Q.is_on_curve() == True)
    assert(pi_1_Q == Q.scalar_mul(px(x)))

    # Q2 = pi2(Q)
    # x = x * beta * (2 * (p^2 - 1) // 6)
    # y = y * beta * (3 * (p^2 - 1) // 6) = -y
    pi_2_Q = G2(
        Q.x.mul(Fp12.beta_pi_2[1]),
        Q.y.mul(Fp12.beta_pi_2[2]),
        Fp2.ONE())
    assert(pi_2_Q.is_on_curve() == True)
    assert(pi_2_Q == Q.scalar_mul(px(x) ** 2))

    # Q3 = pi3(Q)
    # x = x' * beta * (2 * (p^3 - 1) // 6)
    # y = y' * beta * (3 * (p^3 - 1) // 6)
    pi_3_Q = G2(
        Q.x.conjugate_of().mul(Fp12.beta_pi_3[1]),
        Q.y.conjugate_of().mul(Fp12.beta_pi_3[2]),
        Fp2.ONE())
    assert(pi_3_Q.is_on_curve() == True)
    assert(pi_3_Q == Q.scalar_mul(px(x) ** 3))

    alpha, bias = line_add(T, pi_1_Q)
    T = T.add(pi_1_Q)
    L.append((alpha, bias))

    assert(T == Q.scalar_mul(e + px(x)))

    alpha, bias = line_add(T, pi_2_Q.negate())
    T = T.add(pi_2_Q.negate())
    L.append((alpha, bias))

    k = e + px(x) - px(x) ** 2
    assert(T == Q.scalar_mul(k if k > 0 else rx(x) - (-k % rx(x))))

    # alpha, bias = line_add(T, pi_3_Q)
    alpha, bias = Fp2.ZERO(), T.x.mul(T.z.inverse().square())
    T = T.add(pi_3_Q)
    L.append((alpha, bias))

    assert(T == Q.scalar_mul(lamb))
    assert(T.is_infinite() == True)
    assert(alpha.is_zero() == True)

    return L
```

<br />

## Verification of Miller Loop with Precomputed Lines

we only need to do two things:
- line evaluation
    
    just few multiplications on $F_{p^2}$
    
- accumulate evaluations

    arithmetic on $F_{p^{12}}$, it can also be delegated to prover through **Randomized Field Arithmetic**
    
<br />

```python=
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
```

<br />

------

# Verification of Multiplication of $F_{p^k}$

Assume there are two elements of $F_{p^k}$, $a, b \in F_{p^k}$, arithmetics (especially multiplication) on them are still expensive for verifier or within circumstances of **recursive snark**. 

<br />

Let $a \cdot b = c$, where coefficients of $a, b, c$ can be treated as 3 vectors. What we want to do is proving operation result of two vectors $a, b$ is another vector $c$?


If you are familier with **Plonk IOP**, you must have the following statement:
$$
\mathcal{R}(a, b) \overset{?}= c
$$

<br />

Recall that:
- Polynomial Code
    
    Code the three vectors $a, b, c$ into a polynomial $c(X)$ through **Lagrange Interpolation**, same as its evaluation polynomial $e(X)$, and the auxilary index polynomial $h(X)$.

- Zero-check Protocol

    if we want to prove $m(X) - n(X) = 0$, it can be transformed into:
    $$
        m(X) - n(X) = q(X) \cdot Z_H(X)
    $$
    where $Z_H(X)$ is zero polynomial which verifier can be easily computed by itself, and $q(X)$ is the quotient polynomial sent from prover.

this part is so-called **Randomized Field Arithmetic**.

<br />

We simply put screenshot from [On Proving Pairings](https://eprint.iacr.org/2024/640.pdf):

![Screen Shot 2024-05-13 at 12.27.48](https://hackmd.io/_uploads/rJFRUG1Q0.png)

<br />

Last but not least, after utilization of **Randomized Field Arithmetic** (IOP) the **line evaluation** part can be simplified further.

<br />

------

# Verification of Pairings

Put it all together (except for verification of multiplication on $F_{p^{12}}$), provided the witness for final exponentation $c$ and $w_i$ by the prover, the verifier only need to finish following four things:

- precomputation of lines

    offchain task, free-cost
    
- line evaluation with precompuated lines

    few multiplications on $F_{p^2}$, cheap enough
    
- accumulation of line evaluation 

    arithmetics on $F_{p^{12}}$ might be a little expensive, but it can be greatly reduced through **Plonk IOP** 
    
- frobenius map
    
    few multiplications on $F_{p^2}$, cheap enough

<br />

show the code of verification part:
```python=
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
```

<br />

The complete code of parings (including prover-side and verifier-side) you can refer [pairing_verify](https://github.com/PayneJoe/crypto_research/blob/main/docs/on_proving_pairings/bn254/pairing_verify.sage).

<br />

------

# References

[1] Original paper: [On Proving Pairings](https://eprint.iacr.org/2024/640.pdf)

[2] Tonelli Shanks for cubic root: [A remark on the computation of cube roots in finite fields](https://eprint.iacr.org/2009/457.pdf)

[3] Pairing implementaiton: [Theory and Practical Implementation of BLS12-381](https://hackmd.io/@70xfCGp1QViTYYJh3AMrQg/ryo55eEeC)

<br />

----

# Touch

- twitter: @pingzhouyuan
- email: joepayne@163.com