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

chord method applies:
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

tangent method applies:
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
    acc = (0, 0)
    for i in range(l - 1, -1, -1):
        acc = 2 * acc
        if a[i] == 1:
            acc = acc + g
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
$$
$$

<br />

## Fixed Based Single Exponentiation
