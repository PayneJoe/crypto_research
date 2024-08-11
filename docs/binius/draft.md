
- addition is **XOR** operation on bits

- double is **ZERO**, so

    $$
    (x + y)^2 = x^2 + y^2
    $$

- in hash-based scenarios, such as **SHA256/KECCAK256**, binary field is more favorable (efficient) than prime field

<br />

## Multivariate Polynomial and Lagrange 

![multivariate_poly.drawio](https://hackmd.io/_uploads/HJe1MDe9R.png)

<br />

Any $f: \{0, 1\}^n \rightarrow F_p$ can be encoded into $\hat{f}$, which is so-called **MLE**:
$$
\begin{aligned}
\hat{f}(x_1, x_2, ..., x_v) &= \sum_{w \in \{0, 1\}^v} f(w) \cdot \chi_w(x_1, x_2, ..., x_v) \\
\chi_w(x_1, x_2, ..., x_v) &= \prod_{i = 1}^v (x_i \cdot w_i + (1 - x_i) \cdot (1 - w_i)) \\
\end{aligned}
$$
where $\chi_{w}(x_1, x_2, ..., x_v)$ is so-called **identity** polynomial through **Lagrange Interpolation**.

<br />

Evaluations with vector length-$2^n$
$$
e = [f_0, f_1, f_2, ..., f_{2^n - 1}]
$$
can be encoded into $\hat{f}$ defined over entire finte field $F_p$.

<br />

----

Equality **identity** polynomial:
$$
\widehat{eq}(X_0, X_1, ..., X_{l - 1}, Y_0, Y_1, ..., Y_{l - 1}) = \prod_{i = 0}^{l - 1} X_i \cdot Y_i + (1 - X_i) \cdot (1 - Y_i)
$$

<br />

----

There is a point defined over finite field $r = (r_0, r_1, ..., r_{l - 1}) \in K^l$, we can pre-compute its **identity evaluation** table within time $O(2^l)$ through dynamic programming method: 
$$
A^{j}_{(w_1, w_2, ..., w_j)}(x) = A^{j - 1}_{(w_1, w_2, ..., w_{j - 1})}(x) \cdot [w_j \cdot x_j + (1 - w_j) \cdot (1 - x_j)]
$$
where $A^{j}_{(w_1, w_2, ..., w_j)}(r)$ is evaluating **identity polynomial** $\hat{e}_w(x)$ on point $r$ at basis $(w_1, w_2, ..., w_j) \in \mathcal{B}^j$.

<br />

For example, $l = 2$ and $x = r$:

- initialization
    
    $$
    \def\arraystretch{1.5}
       \begin{array}{c:c}
       w_0 & A^0_{(w_0)}(r) \\ \hline
       0 & 1 \\ \hdashline
       1 & 1 \\ \hdashline
    \end{array}
    $$
    
- first iteration for $j = 1$

    $$
    \def\arraystretch{1.5}
       \begin{array}{c:c:c}
       w_0 & w_1 & A^{0}_{(w_0)}(r) & w_1 \cdot x_1 + (1 - w_1) \cdot (1 - x_1) & A^{1}_{(w_0, w_1)}(r) \\ \hline
       0 & 0 & 1 & (1 - r_1) & (1 - r_1) \\ \hdashline
       1 & 1 & 1 & r_1 & r_1 \\ \hdashline
    \end{array}
    $$

- second iteration for $j = 2$

    $$
    \def\arraystretch{1.5}
       \begin{array}{c:c:c}
       w_1 & w_2 & A^1_{(w_1)}(r) & w_2 \cdot x_2 + (1 - w_2) \cdot (1 - x_2) & A^{2}_{(w_1, w_2)}(r) \\ \hline
       0 & 0 & (1 - r_1) & (1 - r_2) & \color{red}{(1 - r_1) \cdot (1 - r_2)} \\ \hdashline
       1 & 0 & r_1 & (1 - r_2) & \color{red}{x_1 \cdot (1 - r_2)} \\ \hdashline
       0 & 1 & (1 - r_1) & x_2 & \color{red}{(1 - r_1) \cdot r_2} \\ \hdashline
       1 & 1 & r_1 & r_2 & \color{red}{r_1 \cdot r_2} \\ \hdashline
    \end{array}
    $$

<br />

So we finally get the point $r = (r_1, r_2, ..., r_l)$ evaluated on **identity polynomial** $\hat{e}_w(x)$, and this evaluation table is what we called **tensor product expansion**, denotes as:
$$
\bigotimes_{i = 0}^{l - 1}(1 - r_i, r_i)
$$

<br />

## Univariate Polynomial and FFT

<br />

## Reed-Solomon Code

$S = \{s_0, s_1, ..., s_{n - 1} \}$ with length $n$, which is a subset of alphabet(finite field) $K$. There's a message with length $k \le n$, then the **Reed-Solomon Code** $RS_{K, S}[n, k]$ can be denoted as:
$$
C = \{ p(s_0), p(s_1), ..., p(s_{n - 1}) | p(X) \in K[X]^{<k} \}
$$
where $p(X)$ is a polynomial defined over finite field $K$ with degree less than $k$, and $p(s_i)$ is $p(X)$ evaluated on point $s_i \in S \subset K$.

<br />

#### Schwartz Lemma 

Two polynomials $p_a(x)$ and $p_b(x)$ with maximal degree-$n$, defined over domain $S$ with size-$p$. Then equal probability of these two polynomials is $\frac{n}{p}$. Using **Schwartz Lemma** we can easily prove that:

<br />

prove :

$$
p_a(x) = p_b(x), x \in S
$$
equally, get roots of polynomial:
$$
p(x) = p_a(x) - p_b(x), x \in S
$$
there will be $n$ roots here, so the probability of $p(x) = 0$ is $\frac{n}{p}$.

<br />

#### Distance of Reed-Solomon Code

- the original message $a = [a_0, a_1, ..., a_{k - 1}]$ will be encoded into a polynomial $t_a(x)$ through **FFT** interpolation, defined over a domain $D_t$ and $|D_t| = k$

- evaluate $t_a(x)$ on a larger domain $D$, where $|D| = r \cdot k = n$, obtain $n = r \cdot k$ evaluations, which is so-called **codeword**
    $$
    a' = [\color{red}{a_0}, a_{0,0}, a_{0,1}, ..., \color{red}{a_1}, a_{1,0}, a_{1,1}, ..., \color{red}{a_{k - 1}}, a_{k - 1, 0}, a_{k - 1, 1}, ...a_{k - 1, r - 1}]
    $$

- **codeword** $a'$ can also be encoded into a polynomial $t'(x)$ through **FFT** interpolation

<br />

There's another message $b = [b_0, b_1, ..., b_{k - 1}]$ with equal length $k$, and its corresponding **codeword** $b' = [\color{red}{b_0}, b_{0,0}, b_{0,1}, ..., \color{red}{b_1}, b_{1,0}, b_{1,1}, ..., \color{red}{b_{k - 1}}, b_{k - 1, 0}, b_{k - 1, 1}, ...b_{k - 1, r - 1}]$.

<br />

If these two messages are equal $a = b$, namely:
$$
[a_0, a_1, ..., a_{k - 1}] = [b_0, b_1, ..., b_{k - 1}] 
$$
then we must have $t_a(x) = t_b(x)$ where $x \in D_t$, namely:
$$
t(x) = t_a(x) - t_b(x) = 0, x \in D_t
$$
there're $k - 1$ roots for polynomial $t(x)$ on domain $D_t$. And consequently their **codeword** are also equal anyway $a' = b'$.

<br />

But if they are different $a \ne b$, then there will be at most $k - 1$ roots for $t(x)$ on domain $D_t$, namely $t_a(\alpha) = t_b(\alpha), \alpha \in D'_t$, where $|D'_t| = k - 1$. For example:
$$
\begin{aligned}
t_a(x) &= 2 + x + x^2 \\
t_b(x) &= x \\
\end{aligned}
$$
then $t(x) = 2 + x^2$, there're at most $k - 1$ roots for $t(x) = 0$, consequently there will be $k - 1$ equal symbols (elements) between the two **codeword** $a'$ and $b'$. 

<br />

Finally the Hanming Distance of **Reed-Solomon Code** will be $n - (k - 1)$.

<br />

## Binary Towers

Binary field $T_2$ is extended from binary field $T_1$:
$$
T_2 = T_1[X_2] / X_2^2 + X_2 \cdot X_1 + 1
$$
where $X_1$ is the generator of binary extension field $T_1$, and $X_2$ is the generator of binary extension field $T_2$.

<br />

snippet code for binary towers:
```pyhton3
## T0
T0 = GF(2)
print('T0: {}'.format(T0.list()))

## T1, T2, T3, T4
X1 = T0.gen()
X2 = T0['X1'].gen()
T_extends = [T0]
for i in range(1, 5):
    var = 'X{}'.format(i)
    pol = X2^2 + X1 * X2 + 1
    assert(pol.is_irreducible() == True)
    Ti = T_extends[-1].extension(pol, var)
    T_extends.append(Ti)
    X1 = T_extends[-1].gen()
    X2 = T_extends[-1][var].gen()
```

<br />

## Arithmetics of Binary Field

#### Multiplication

Assume binary field $F_2^{(n)} = F_2^{(n - 1)}[X_n] / X_n^2 + X_n \cdot X_{n - 1} + 1$, two binary field element $a$ and $b$ denoted as:
$$
\begin{aligned}
a &= a_0 + a_1 \cdot X_n \\
b &= b_0 + b_1 \cdot X_n \\
\end{aligned}
$$
where $a_0$ and $a_1$ are the lower bits and higher bits of $a$ repectively. For example, $a = 14 = 2 + 3 * 2^2$, where $a_0 = 2, a_1 = 3$.

<br />

Let's take a look at how multiplication works:
$$
\begin{aligned}
a \cdot b &= (a_0 + a_1 \cdot X_n) \cdot (b_0 + b_1 \cdot X_n) \\
&= a_0 \cdot b_0 + a_1 \cdot b_1 \cdot X_n^2 + (a_0 \cdot b_1 + a_1 \cdot b_0) \cdot X_n \\ 
\end{aligned}
$$
since $X_n^2 = X_{n - 1} \cdot X_n + 1$, we have:
$$
a \cdot b = (\boxed{a_0 \cdot b_0 + a_1 \cdot b_1}) + (\boxed{a_1 \cdot b_1 \cdot X_{n - 1} + a_0 \cdot b_1 + a_1 \cdot b_0}) \cdot X_n
$$
where the result consists two parts, both of them are elements of $F_{n - 1}$. 

<br />

According to **Karatsuba** method:
$$
(a_0 \cdot b_1 + a_1 \cdot b_0) = (a_0 + a_1) \cdot (b_0 + b_1) - a_0 \cdot b_0 - a_1 \cdot b_1 \\
$$
since binary field $-x = x$, so we have: 
$$
(a_0 \cdot b_1 + a_1 \cdot b_0) = (a_0 + a_1) \cdot (b_0 + b_1) + a_0 \cdot b_0 + a_1 \cdot b_1
$$

<br />

If we ignore the cost of addition, the total cost of $F_n$ multiplication equals **FOUR** times of cost $F_{n - 1}$ multiplication.
$$
\def\arraystretch{1.5}
   \begin{array}{c:c:c}
   Symbol & Expr \\ \hline
   T0 & a_0 \cdot b_0 \\
   T1 & a_1 \cdot b_1 \\
   T2 & T1 \cdot X_{n - 1} \\
   T3 & (a_0 + a_1) \cdot (b_0 + b_1) \\
\end{array}
$$
the result will be:
$$
a \cdot b = (T0 + T1) + (T2 + T3 + T0 + T1) \cdot X_n
$$

<br />

## References

[1] [Succinct Arguments over Towers of Binary Fields](https://eprint.iacr.org/2023/1784.pdf)