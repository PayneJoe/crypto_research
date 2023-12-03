## Definition

single exponentiation:
$$
g^a
$$

where:
$$
g \in E(\mathbb{F}_p), a \in \mathbb{Z}_q
$$

> we can also say that $a$ is over scalar field $\mathbb{F}_q$, point $g$ is over base field $\mathbb{F}_p$

<br />

if $g$ and $a$ both vary, we call it **variable based exponentiation**; if only $a$ varies, we call it **fixed based exponentiation**.

<br />

## Group Law on Weierstrass Model

#### point addition

condition that:
$$
x_1 \ne x_2 \leftrightarrow P \ne Q
$$

**chord method** applies:
$$
\begin{aligned}
\lambda &= \frac{y_2 - y_1}{x_2 - x_1} \\
\downarrow \\
x_3 &= \lambda^2 - x_1 - x_2 \\
y_3 &= \lambda (x_1 - x_3) - y_1 \\
\downarrow \\
R &= (x_3, y_3)
\end{aligned}
$$

<br />

field operation cost would be:
$$
6A + 2M + 1I
$$

<br />

#### point doubling

condition that:
$$
x_1 = x_2, y_1 = y_2
$$

**tangent method** applies:
$$
\begin{aligned}
\lambda &= \frac{3x_1^2 + a}{2y_1} \\
\downarrow \\
x_3 &= \lambda^2 - 2x_1 \\
y_3 &= \lambda (x_1 - x_3) - y_1 \\
\downarrow \\
R &= (x_3, y_3) \\
\end{aligned}
$$

<br />

field operation cost:
$$
4A + 2M + 1I
$$

<br />

## Variable Based Single Exponentiation

#### vallina (left-to-right) method 
the input, one is binary representation of a scalar field $a$, the other is point $g$ on elliptic curve:
$$
a = (a_{l - 1},...,a_0)_2 \in \mathbb{F}_q, g \in E(\mathbb{F}_p)
$$

expected output is:
$$
[a]g
$$

<br />

the **left-to-right** algorithm would be like:
```python=
def single_exponentiation(a, g):
    acc = AffinePoint(0, 0)
    for i in range(l - 1, -1, -1):
        acc = 2 * acc
        if a[i] == 1:
            acc = acc + g
    return acc
```

<br />

**group operation cost** would be:
$$
\frac{l}{2} * A + l * D
$$

<br />

apply the addition law of Weierstrass Model, the **field operation cost** would be:
$$
\frac{l}{2} * (6A + 2M + 1I) + l * (4A + 2M + 1I)
$$

<br />

:::info
Notice:
You may have already found that the *doubling* and *addition* operations were merely equal, as the number of *addtion* operation is the non-zeros in $a_i$.

So, do we have a more *fair* solution and make them even?
:::

<br />

#### NAF method

NAF, so-called non-adjacent form, is a another encoding/representation method for scalar $a$, beyound binary representation.

<br />

The expectation of NAF is simple:
$$
\begin{align}
a_i \in \{-1, 0, +1\} \\
a_i * a_{i + 1} = 0 \\
a_{l - 1} \ne 0 \\
\end{align}
$$

<br />

The property of NAF is:
- $NAF(k)$ is unique regarding a scalar $k$

- has the fewest nonzero digits, including $-1$ and $+1$

- at most **one more than** the length of binary representation

- $k$ falls into a fixed range $\frac{2^l}{3} < k < \frac{2^{l + 1}}{3}$

- the average number of non-zeros is $\frac{l}{3}$, less than $\frac{l}{2}$ in *binary representation*
    
    this is why NAF representation is superier than *binary representation*

<br />

So, how to encode scalar $a$, here is how it goes:
```python=
def NAF(a):
    i = 0
    naf_a = []
    while a >= 1:
        if a % 2 == 1:
            a_i = 2 - k % 4
        else:
            a_i = 0
        a = (a - a_i) / 2
        i += 1
        naf_a.append(a_i)
    return naf_a
        
```

<br />

Next, we will show how exponentiation goes with NAF represented scalar $a$:
```python=
def single_exponentiation(naf_a, g):
    acc = AffinePoint(0, 0)
    for i in range(l - 1, -1, -1):
        acc = 2 * acc
        if naf_a[i] == 1:
            acc = acc + g
        elif naf_a[i] == -1:
            acc = acc - g  
    return acc
```
<br />

**group operation cost** would be:
$$
\frac{l}{3} * A + l * D
$$

:::info
Comparation:
- *doubling operation* is same as that of *vallina method*
- *addition operation* is slightly less than that of *vallinda mehod*, $\frac{l}{3}$ verus $\frac{l}{2}$
:::

<br />

Still one more question exists, why the average number of non-zeros is $\frac{l}{3}$ or the probability of that is 1/3? Let's prove that...

Assume the number zeros between two non-zero digits is $t$:

- $l = 3$
    $$
    \begin{aligned}
        &a_2 = \pm 1, a_1 = 0 \rightarrow P(a_0 = 1) = \frac{1}{2} \\
        &\Downarrow \\
        &P(t = 1) = \frac{1}{2}
    \end{aligned}
    $$

- $l = 4$
    $$
    \begin{aligned}
        &a_3 = 1, a_2 = 0 \rightarrow P(a_1 = 0, a_0 = 1) = \frac{1}{4} \\
        &\Downarrow \\
        &P(t = 2) = \frac{1}{4}
    \end{aligned}
    $$

- $l = 5$
    $$
    \begin{aligned}
        &a_4 = 1, a_3 = 0 \rightarrow P(a_2 = 0, a_1 = 0, a_0 = 1) = \frac{1}{8} \\
        &\Downarrow \\
        &P(t = 3) = \frac{1}{8}
    \end{aligned}
    $$

...

So, we can summarize that the expected zeros between two non-zero digits is:
$$
\begin{aligned}
E(t) &= 1 * \frac{1}{2} + 2 * \frac{1}{4} + 3 * \frac{1}{8} + ... = \sum_{i = 1}^\inf \frac{i}{2^i} \\
\\
2E(t) - E(t) &= E(t) = 1 + \sum_{i = 1}^\inf \frac{1}{2^i} = 2 \\
\end{aligned}
$$

therefore the **expected** NAF sequence digits would be like:
$$
1 00 1 00 1 00 ...
$$

Last of all, we can conclude that the probability of non-zero digits is $\frac{1}{3}$, or the average number of non-zeros is $\frac{l}{3}$.

<br />

#### Windowed NAF method

If there is a representation promising that at most one non-zero digit exsits within $w$ length of digit sequence, the *group addition* cost would be cheapter:
$$
\frac{l}{w + 1} \le \frac{l}{2 + 1}
$$

<br />

Let's see what the window-$w$ NAF looks like:
```python=
def window_naf(w, a):
    w_naf_a = []
    i = 0
    while a >= 1:
        if a % 2 == 1:
            a_i = a % pow(2, w)
        else:
            a_i = 0
        a = (a - a_i)/2
        w_naf_a.append(a_i)
        i += 1
    return w_naf_a
```

<br />

the NAF elements $a_i$ will be odd number, we need to do some precomputation to accomplish exponentiation:

```python=
def single_exponentiation(w, w_naf_a, g):
    p = [i * g for i in range(1, pow(2,w - 1), 2)]
    acc = AffinePoint(0, 0)
    for i in range(l - 1, -1, -1):
        acc = 2 * acc
        if w_naf_a[i] > 0:
            acc = acc + p[w_naf_a[i]]
        elif w_naf_a[i] < 0:
            acc = acc - p[w_naf_a[i]]
    return acc
```

<br />

**group operation** cost would be:
$$
((2^{w - 2} - 1) * A + 1 * D) + (\frac{l}{w + 1} * A + l * D)
$$

<br />


## Fixed Based Single Exponentiation
