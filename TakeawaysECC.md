> This note is mostly adopted from "A graduate course in applied cryptography".

<br />

## Points on Elliptic Curve

There are two methods to find the third point $W$ on curve provided two known points $U$ and $V$:

- **Chord Method** to find $W$, addition of different points with condition
    $$
    \begin{aligned}
        x_1 \ne x_2 &\leftrightarrow U \ne \pm V \\
        \Downarrow \\
        U \boxplus V &= -W
    \end{aligned} 
    $$
- **Tangent Method** to find $W$, addition of same point with condition
    $$
    \begin{aligned}
        \{x_1 = x_2, y_1 \ne 0\} &\leftrightarrow U = V \\
        \Downarrow \\
        U \boxplus V &= -W
    \end{aligned}
    $$

<br />

## Elliptic Curves on Finite Field

First of all, we'll talk about <span style="color:red">Weierstrass form</span> equation of elliptic curve $E$ over finite field $\mathbb{F}_p$, denoted by $E/\mathbb{F}_p$:

$$
y^2 = x^3 + ax + b
$$

with condition $4a^3 + 27b^2 \ne 0$, so that **there would not have two roots for $y = 0$**, only one exists!

<br />

The curve defined above we often call it <span style="color:red">Weierstrass Curve</span>.

<br />

#### base field VS extension of base field
$$
(x_1, y_1) \in \mathbb{F}_{p^e}
$$

- $e = 1$, point $(x_1, y_1)$ is defined over base field $\mathbb{F}_p$
- $e > 1$, point $(x_1, y_1)$ is defined over extension of base field $\mathbb{F}_p$

<br />

#### order of elliptic curve

order of elliptic curve or the number of points on elliptic curve is:

$$
|E(\mathbb{F}_p)| = p^e + 1 - t, |t| \le 2 \sqrt{p^e}
$$

eg. when $p = 11, e = 1$, we get $t = 0$:
$$
E(\mathbb{F}_p) = \{\mathcal{O}, (1, 0), (0, \pm 1), (2, \pm 3), (5, \pm 4), (7, \pm 5), (9, \pm 2)\}
$$

including the **identity** element $\mathcal{O}$.

<br />

:::info

#### How to count?

The statistic method for the number of points is called **SEA algorithm**, whose cost is in time polynomial in $\log(p^e)$, even for a large prime $p$. 

Looks impressive, right?

:::

<br />

#### the addition law

- $x_1 \ne x_2$, two different points, using **chord method**
    - compute the slope of chord line
        $$
            s_c = \frac{y_2 - y_1}{x_2 - x_1} \\
        $$
    - compute new point
        $$
        \begin{aligned}
            x_3 &= s_c^2 - (x_1 + x_2) \\
            y_3 &= s_c * (x_1 - x_3) - y_1 \\
        \end{aligned}
        $$
- $x_1 = x_2$ and $y_1 \ne 0$, same point, using **tangent method**
    - compute the slope of tangent line
        $$
            s_t = \frac{3x_1^2 + a}{2y_1}\\
        $$
    - compute new point
        $$
        \begin{aligned}
            x_3 &= s_t^2 - 2x_1 \\
            y_3 &= s_t * (x_1 - x_3) - y_1 \\
        \end{aligned}
        $$
- $x_1 = x_2$ and $y_1 = -y_2$, negative each other, resulting is the identity element $\mathcal{O}$
    $$
    \begin{aligned}
        x_3 &= 2x_1 \\
        y_3 &= 0 \\
    \end{aligned}
    $$
    
:::info
#### About Identity Element

1. $\mathcal{O}$ occurs when $y = 0$ and only one $x$ exists
2. if $\mathcal{O} \ne P = (x_1, y_1) \in \mathbb{F}_{p^e}$, there must have $-P = (x_1, -y_1)$
:::

<br />

Another two types of curves will be discussed besides **Weierstrass Curve**.

<br />

## Montgomery Curve

$$
Bv^2 = u^3 + A u^2 + u
$$

where $B * (A^2 - 4) \ne 0$.

<br />

#### Pros

Montgomery Curve has computational benefit compared with Weierstrass Curve.

<br />

#### Features

- conversion from **Montgomery Curve** to **Weierstrass Curve**
    $$
    \begin{aligned}
        u &= Bx - A/3 \\
        v &= By \\
    \end{aligned}
    $$
- order of **Montgomery Curve** is always divisable by 4
- not all **Weierstrass Curve** can be converted to **Montgomery Curve**
    
    eg. if the order of Weierstrass Curve is odd, not divisable by 4. P256 Curve is a good example
 
<br />
 
## Edward Curve

$$
x^2 + y^2 = 1 + d x^2 y^2
$$
    
where $d \ne 0, 1$.

<br />

#### Pros
- chord and tangent addition law is easy, no edge cases to consider, simple enough!

    $$
    P \boxplus Q = (\frac{x_1 y_2 + x_2 y_1}{1 + d x_1 x_2 y_1 y_2}, \frac{y_1 y_2 - x_1 x_2}{1 - d x_1 x_2 y_1 y_2})
    $$
- resisting to timing attack

<br />

#### Features

- identity element $\mathcal{O}$ is $(0, 1)$
- if $P = (x_1, y_1)$, then $-P = (-x_1, y_1)$
- converstion between **Weierstrass Curve** and **Edward Curve** is possible
- order of **Edward Curve** is divisible by 4 same as **Montgomery Curve**

<br />

## Discrete Log Problem
On elliptic curve $E/\mathbb{F}_p$, the order of curve $|E(\mathbb{F}_p)| = q$, suppose we have:
$$
\begin{aligned}
Q &= \alpha P \\
\text{or} \\
Q &= P^{\alpha} \\
\end{aligned}
$$

if we want to solve $\alpha$ provided the two points $P$ and $Q$:
$$
\alpha = \log_{P} Q
$$

the time cost of solving this is $O(\sqrt{q})$, we often call it **Discrete Log Problem**.

<br />

#### Insecure Curves

Many curves can solve it in time cost $O(\sqrt{q})$, but exception always comes. There are three kinds of insecure curves:

- all the prime factors of $|E(\mathbb{F}_p)|$ is small
    
    in this condition we must make sure $|E(\mathbb{F}_p)|$ is $q$, $4q$, or $8q$, where $q$ is a **prime**

- $|E(\mathbb{F}_p)| = p$ is unacceptable in cryptos

    the cost of solving DL problem is always in polynomial time

- $|E(\mathbb{F}_p)| = p$ is divisible by $p^{\lambda} - 1$, where $\lambda$ is a small integer

    DL problem can be solved within serval hours, so $\lambda$ must be large enough
 
<br />

## Secp256r1 Curve
Secp256r1 is widely used in internet protocols, the curve is randomly selected, not determined.

#### Definition

- prime modulo $p = 2^{256} - 2^{224} + 2^{192} + 2^{96} - 1$

- Weierstrass form curve $y^2 = x^3 -3x + b$, where $b$ is a 255 bits integer
    $$
        b = 5ac635d8 aa3a93e7 b3ebbd55 769886bc 651d06b0 cc53b0f6 3bce3c3e 27d2604b
    $$
    
- $|E(\mathbb{F}_p)|$ is a prime number, which is also close to $p$

- seed $S$ to generate $b$, so the curve will be determined through this process

<br />

## Secp256k1 Curve

Secp256k1 is widely used in blockchain systems.


#### Definition

- prime modulo $p = 2^{256} - 2^{32} - 2^9 - 2^8 - 2^7 - 2^6 - 2^4 - 1$

- Weierstrass form curve $y^2 = x^3 + 7$, no big number $b$ any more...

- $|E(\mathbb{F}_p)|$ is a prime number, which is also close to $p$

#### Properties

- Koblitz curve

    the magic of **map** functionality is:
    $$
    \begin{aligned}
        \phi(P) &= \lambda P \\
        \Downarrow \\
        \mathcal{O} &= P + \phi(P) + \phi^2(P) \\
        &= (1 + \lambda + \lambda^2)P \\
        \Downarrow \\
        0 &= (1 + \lambda + \lambda^2)
    \end{aligned}
    $$
    
    we can take advantage of **map** to accelerate computation of scalar multiplication:
    $$
    \begin{aligned}
        \alpha P &= (\tau_0 + \tau_1 \lambda + \tau_2 \lambda^2) P \\
        &= \tau_0 P + \tau_1 \phi(P) + \tau_2 \phi^2(P) \\
    \end{aligned}
    $$

    this is so-called **multi-exponentiation**, which is two or three times faster than computing $\alpha P$ directly. Amazing, is it?
    
- Cycle curves accompanied with Secq256k1
    
    order of Secp256k1 is $q = |E(\mathbb{F}_p)|$, and order of Secq256k1 is $p = |E(\mathbb{F}_q)|$
    
<br />
    
#### Security

Since the $|q = E(\mathbb{F}_p)|$ is close to $p$, and both of them are a prime closely to $2^{256}$. So the time cost of solving the DL problem is :
$$
O(\sqrt{p}) \approx O(2^{128})
$$

<span style="color:red">Therefore we can say the security thrength of Secp256k1 curve is 128 bits</span>.

<br />

#### Group Membership Check

Clearly we could see that the addition law of **Weierstrass Curve** is nothing about the constant parameter **b**. 

Another words, if we do not check that the new point $(x_3, y_3)$ satisfies the euqation:
$$
y_3^2 \overset{?}= x_3^2 + ax_3 + b \\
$$

, we can not assure the new added point $(x_3, y_3)$ is actually **on the curve**.

<br />

The well-known **Invalid Curve Attack** is a typical case for that. 


<br />

#### Twist Secure

Curve $E/\mathbb{F}_p$:
$$
y^2 = x^3 + ax + b
$$
has a related curve $\widetilde{E}/\mathbb{F}_p$:
$$
cy^2 = x^3 + ax + b
$$

their order satisfies:
$$
|E(\mathbb{F}_p)| + |\widetilde{E}(\mathbb{F}_p)| = 2p + 2
$$

since $|E(\mathbb{F}_p)| = p + 1 - t$, so $|\widetilde{E}(\mathbb{F}_p)| = p + 1 + t$.


<br />

We can say cuve $E/\mathbb{F}_p$ is **twist secure** if both solving the DL problem for $E/\mathbb{F}_p$ and $\widetilde{E}/\mathbb{F}_p$ is hard. 

Another words, <span style="color:red">both orders of $E/\mathbb{F}_p$ and $\widetilde{E}/\mathbb{F}_p$ must be (big) prime numbers, or small times of big prime number</span>.

<br />

Unfortuntely both Secp256r1 and Secp256k1 curves are not **twist secure**:
- Secp256r1 curve, time cost of solving DL problem for twist Secp256r1 is $\sqrt{34905} = 187$ times easier than Secp256r1
- Secp256k1 curve, time cost of solving DL problem for twist Secp256k1 is $2^{18}$ times easier than Secp256k1 

<br />

## Curve25519

This is a optimized curve on group operation, and it is also twist secure.

<br />

Montgomery curve equation:
$$
By^2 = x^3 + A x^2 + x
$$

where $B = 1, A = 486662$.

<br />

##### Features

-  prime modulo is $p = 2^{255} - 19$, the largest prime number less than $2^{255}$
-  fast scalar multiplication
-  $|E(\mathbb{F}_p)|$ is divisible by 4
    
    more specifically to say it is 8 times of a big prime number, therefore it is **twist secure** too
    
- generator is $(x_1, y_1)$ where $x_1 = 9$

<br />

It has a long story behind the number 486662, two major aspects need to be considered:
1. more smaller A is, faster group operation
2. both order of curve25519 or twist curve25519 are either $4\cross$ of a big prime number of $8\cross$ of that

<br />

## Pairing

three cycle groups $\mathbb{G}_0, \mathbb{G}_1, \mathbb{G}_T$ share the same group order $p$

#### symmetric vs asymmetric

- DDH problem not hold in **symmetric pairing**, $\mathbb{G}_0 = \mathbb{G}_1$
    

- in **asymmetric pairing**, solving discrete log problem must be hard for $\mathbb{G}_T$, or it would be easy to get the solve for $\mathbb{G}_0$ and $\mathbb{G}_1$

<br />

so asymmetric pairing is prefered.

<br />

#### properties of three groups 

- $\mathbb{G}_0$, subgroup $E(\mathbb{F}_p)$ with prime order $q$
- $\mathbb{G}_1$, subgroup $E(\mathbb{F}_{p^d})$ with prime order $q$,  where $d > 0$ and $\mathbb{G}_0 \cap \mathbb{G}_1 = \mathcal{O}$
- $\mathbb{G}_T$, **multiplicative subgroup** over $\mathbb{F}_{p^d}$ with prime order $q$

<br />

#### about embedding degree $d$

- paring friendly, $d < 16$, need to be small for computing efficiency
- since $p^d < p$, the bit size of $\mathbb{F}_p$ is smaller than $\mathbb{F}_{p^d}$, so computing on $\mathbb{G}_0$ is much cheaper than $\mathbb{G}_1$

<br />

#### bn256 vs bls381

- both have embedding degree $d = 12$
- both have group order $q$ with bit size 256
- bit size of prime field $\mathbb{F}_p$ in bn254 is 254, and bls381 is 381
- bls381 is more secure than bn256
    
    solving the DL problem for $\mathbb{G}_T$ in bls381 is harder, while computing pairing is more expensive

<br />

#### pairing and exponentiation