
#### Even Parity Check

- can only be detect sigle error
- can not fix error

<br />

#### MLE Decoding

- detect errors
- correct errors

with probability distribution

<br />

#### Block Codes

Lemma :
$$
d(\bold{x}, \bold{y}) \le d(\bold{x}, \bold{z}) + d(\bold{z}, \bold{y})
$$

<br />

- the minimus distance of codewords must be greater than tolerance of error correction 

- a code scheme with $d_{min} = 2 n + 1$, then it can correct at most $n$ errors, and detect at most $2n$ errors

<br />

#### Linear Code

Group code is a set of codewords consisting a group, which including the **zero** or **identity** element.

- distance of code is defined to be the minimum distance between any two codewords, $d_{min}$, it determines the capability of detecting and correcting errors

- weight of code is defined to be the minimus number of $1$ within all codewords except *zero-codeword*, $w$

- compute the distance of codewords is costly, say $\binom{n}{2} = \frac{n (n - 1)}{2}$, fortunately $d_{min}$ is dependent with weights

    $$
    d(\bold{x}, \bold{y}) = w(\bold{x} + \bold{y})
    $$

<br />

Theorem:
$$
d_{min} = min\{w(\bold{x}): \bold{x} \ne \bold{0} \}
$$
the distance of **group code** is the minimum weight of non-zero codewords, therefore the cost of computing distance becomes $n$, greatly reduced from $\frac{n (n - 1)}{2}$.

<br />

How to generate a good group code?

<br />

Theorem:

**Let $H$ be a matrix $\mathbb{M}_{m \times n}(\mathbb{Z}_2)$, then its null space is a group code**:
$$
Null(H) = \{\bold{x}: \langle H \cdot \bold{x} \rangle = \bold{0} \}
$$

For example:
$$
H = 
\begin{bmatrix} 
0 & 1 & 0 & 1 & 0 \\
1 & 1 & 1 & 1 & 0 \\
0 & 0 & 1 & 1 & 1 \\
\end{bmatrix}
$$

then the group code would be $\bold{x} = \{ 00000, 11110, 10101, 01011 \}$, with distance $3$.

<br />

But how to detect errors?

If the received word is $10001$, then we can easily check that it's not a codeword, since $H \cdot \begin{bmatrix} 1 \\ 0 \\ 0 \\ 0 \\ 1 \end{bmatrix} = \begin{bmatrix} 0 \\ 1 \\ 1 \end{bmatrix} \ne \bold{0}$.

<br />

Is there a efficient method to generate group code?

<br />

#### Generator Matrix

Original words are represented as column vector $\bold{x}_{n - m, 1}$.

<br />

Split $H_{m \times n}$, $n > m$, into two parts:
$$
H_{m \times n} = (A | I_m)
$$
where $A$ is $m \times (n - m)$ partial matrix, and $I_m$ is a $m \times m$ identity matrix.

<br />

Let generator matrix:
$$
G_{n \times (n - m)} = (\frac{I_{n - m}}{A})
$$

<br />

Then we must have:
$$
H \cdot \bold{y} = \bold{0}
$$
where $\bold{y} = \langle G, \bold{x} \rangle$ is what we want, group code. Why it holds for that? Let's prove it:
$$
\begin{aligned}
H \cdot \bold{y} &= (A | I_m) \cdot (\frac{I_{n - m}}{A}) \cdot \bold{x} \\
&= (A | I_m) \cdot (\frac{\bold{x}}{A \cdot \bold{x}}) \\
&= \langle A, \bold{x} \rangle + \langle I_m, A \cdot \bold{x} \rangle \\
&= \langle A, \bold{x} \rangle + \langle A, \bold{x} \rangle  \\
&= \bold{0}
\end{aligned}
$$

<br />

For example, we have $3$-bits message for transmission, which $n - m = 3$, and assume length of codeword is $n = 2 * 3 = 6$, and partial matrix: 
$$
A = 
\begin{bmatrix}
0 & 1 & 1 \\
1 & 1 & 0 \\
1 & 0 & 1 \\
\end{bmatrix}
$$

then we can obtain the entire group code $\bold{y}$ with above method:
$$
\def\arraystretch{1.5}
   \begin{array}{c:c}
   \bold{x} & G \cdot \bold{x} \\ \hline
   000 & 000\color{red}{000} \\ \hdashline
   001 & 001\color{red}{101} \\ \hdashline
   010 & 010\color{red}{110} \\ \hdashline
   011 & 011\color{red}{100} \\ \hdashline
   100 & 100\color{red}{011} \\ \hdashline
   101 & 101\color{red}{001} \\ \hdashline
   110 & 110\color{red}{010} \\ \hdashline
   111 & 111\color{red}{000} \\ \hdashline
\end{array}
$$
the time complexity of this is $O(2^m)$.

<br />

Let's take a closer look at:
$$
\begin{aligned}
\bold{y} &= G \cdot \bold{x} = (\frac{I_{n - m}}{A}) \cdot \bold{x}  \\
&= (\frac{\bold{x}}{A \cdot \bold{x}}) \\
\end{aligned}
$$
the above part of $\bold{y}$ is orginal message $\bold{x}$ **information bits**, the below part of $\bold{y}$ is **check bits**.