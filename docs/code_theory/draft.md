
> This is note for **Code Theory**-)

<br />

## Even Parity Check

- can only be detect sigle error
- can not fix error

<br />

## MLE Decoding

- detect errors
- correct errors

with probability distribution

<br />

## Block Codes

**Lemma**:
$$
d(\mathbf{x}, \mathbf{y}) \le d(\mathbf{x}, \mathbf{z}) + d(\mathbf{z}, \mathbf{y})
$$

<br />

Two properties:
- the minimus distance of codewords must be greater than tolerance of error correction 

- a code scheme with $d_{min} = 2 n + 1$, then it can correct at most $n$ errors, and detect at most $2n$ errors

<br />

## Linear Code

Group code is a set of codewords consisting a group, which including the **zero** or **identity** element.

- distance of code is defined to be the minimum distance between any two codewords, $d_{min}$, it determines the capability of detecting and correcting errors

- weight of code is defined to be the minimus number of $1$ within all codewords except *zero-codeword*, $w$

- compute the distance of codewords is costly, say $\binom{n}{2} = \frac{n (n - 1)}{2}$, fortunately $d_{min}$ is dependent with weights

    $$
    d(\mathbf{x}, \mathbf{y}) = w(\mathbf{x} + \mathbf{y})
    $$

<br />

**Theorem**:
$$
d_{min} = min\{w(\mathbf{x}): \mathbf{x} \ne \mathbf{0} \}
$$
the distance of **group code** is the minimum weight of non-zero codewords, therefore the cost of computing distance becomes $n$, greatly reduced from $\frac{n (n - 1)}{2}$.

<br />

How to generate a good group code?

<br />

**Theorem**:

**Let $H$ be a matrix $\mathbb{M}_{m \times n}(\mathbb{Z}_2)$, then its null space is a group code**:
$$
Null(H) = \{\mathbf{x}: \langle H \cdot \mathbf{x} \rangle = \mathbf{0} \}
$$

<br />

For example:
$$
H = 
\begin{bmatrix} 
0 & 1 & 0 & 1 & 0 \\
1 & 1 & 1 & 1 & 0 \\
0 & 0 & 1 & 1 & 1 \\
\end{bmatrix}
$$

then the group code would be $\mathbf{x} = \{ 00000, 11110, 10101, 01011 \}$, with distance $3$.

<br />

But how to detect errors?

If the received word is $10001$, then we can easily check that it's not a codeword, since $H \cdot \begin{bmatrix} 1 \\ 0 \\ 0 \\ 0 \\ 1 \end{bmatrix} = \begin{bmatrix} 0 \\ 1 \\ 1 \end{bmatrix} \ne \mathbf{0}$.

<br />

Is there a efficient method to generate group code?

<br />

## Generator Matrix

Original words are represented as column vector $\mathbf{x}_{n - m, 1}$.

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
H \cdot \mathbf{y} = \mathbf{0}
$$
where $\mathbf{y} = \langle G, \mathbf{x} \rangle$ is what we want, group code. Why it holds for that? Let's prove it:
$$
\begin{aligned}
H \cdot \mathbf{y} &= (A | I_m) \cdot (\frac{I_{n - m}}{A}) \cdot \mathbf{x} \\
&= (A | I_m) \cdot (\frac{\mathbf{x}}{A \cdot \mathbf{x}}) \\
&= \langle A, \mathbf{x} \rangle + \langle I_m, A \cdot \mathbf{x} \rangle \\
&= \langle A, \mathbf{x} \rangle + \langle A, \mathbf{x} \rangle  \\
&= \mathbf{0}
\end{aligned}
$$
Besides, we also have $H \cdot G = \mathbf{0}$:
$$
H \cdot G = (A | I_m) \cdot (\frac{I_{n - m}}{A}) = \langle A, I_{n - m}\rangle + \langle I_m, A\rangle = A + A = \textbf{0}
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

then we can obtain the entire group code $\mathbf{y}$ with above method:
$$
\def\arraystretch{1.5}
   \begin{array}{c:c}
   \mathbf{x} & G \cdot \mathbf{x} \\ \hline
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
\mathbf{y} &= G \cdot \mathbf{x} \\ 
&= (\frac{I_{n - m}}{A}) \cdot \mathbf{x}  \\
&= (\frac{\mathbf{x}}{A \cdot \mathbf{x}}) \\
\end{aligned}
$$
the above part of $\mathbf{y}$ is orginal message $\mathbf{x}$ **information bits**, the below part of $\mathbf{y}$ is **check bits**.

<br />

## General Decode

There's an important concept in decoding a receiving word, **syndrome**. For example, there's a binary matrix:
$$
H = 
\begin{pmatrix}
1 & 1 & 1 & 0 & 0 \\
0 & 1 & 0 & 1 & 0 \\
1 & 0 & 0 & 0 & 1 \\
\end{pmatrix}
$$

And two receiving words, $\mathbf{x_1} = \begin{pmatrix} 1 \\ 1 \\ 0 \\ 1 \\ 1 \end{pmatrix}$ and $\mathbf{x_2} = \begin{pmatrix} \color{red}{0} \\ 1 \\ 0 \\ 1 \\ 1 \end{pmatrix}$. Then we can prove that $\mathbf{x_1}$ is correct codeword, while the $\mathbf{x_2}$ is not. Since we have:
$$
H \cdot \mathbf{x_1} = \mathbf{0} \\ 
H \cdot \mathbf{x_2} = \begin{pmatrix} 1 \\ 0 \\ 1 \end{pmatrix} \ne \mathbf{0}
$$
and $H \cdot \mathbf{x_2}$ is just the first column of $H$, we can conclude that the error term is first element of $\mathbf{x_2}$, if we flip the $0$ to $1$, then it will the very $\mathbf{x_1}$.

<br />

So, we call $H \cdot \mathbf{x}$ is the **syndrome** of receiving word $\mathbf{x}$.

<br />

In general, 
$$
\mathbf{x} = \mathbf{c} + \mathbf{e}
$$ 
where $\mathbf{x}$ is receiving word, $\mathbf{c}$ is sending codeword, $\mathbf{e}$ is error term. Since $H \cdot \mathbf{c} = \mathbf{0}$, so $H \cdot \mathbf{x} = H \cdot (\mathbf{c} + \mathbf{e}) = H \cdot \mathbf{e}$. In another words, **syndrome** is independent of sending codeword, it's about the error term.

<br />

## Coset Decode

Assuming we have the same binary matrix:
$$
H = 
\begin{pmatrix}
1 & 1 & 1 & 0 & 0 \\
0 & 1 & 0 & 1 & 0 \\
1 & 0 & 0 & 0 & 1 \\
\end{pmatrix}
$$
we can easily obtain the coresponding codewords $C = G \cdot \mathbf{x}$ with the help of **generator matrix** $G$:
$$
\def\arraystretch{1.5}
   \begin{array}{c:c}
   \mathbf{x} & G \cdot \mathbf{x} \\ \hline
   00 & 00000 \\ \hdashline
   01 & 01101 \\ \hdashline
   10 & 10011 \\ \hdashline
   11 & 11110 \\ \hdashline
\end{array}
$$

<br />

Since $C$ is *group code* with order $2^2$, which is also a subgroup of $\mathbb{Z}_2^5$. So there're $2^{5 - 2} = 2^3$ cosets. How does the coset come from? that's 
$$
\mathbb{Z}_2^5 / C = \{\mathbf{e} + C: \mathbf{e} \in \mathbb{Z}_2^5\}
$$
where $\mathbf{e}$ is elements from $\mathbb{Z}_2^5$ distinct from $C$.
According the definition of **coset**, we can also say that $\mathbb{Z}_2^5$ is divided into $2^3$ cosets by group $C$, each coset with order $2^2$. And $\mathbf{e}$ is what we called the **representative** of coset.

<br />

Coset is not we want to emphasize, the relation between coset representative $\mathbf{e}$ and **syndrome** $\mathbf{s}$ is, so what is the relationship beween them? 

<br />

Assume $\mathbf{x} = \mathbf{c} + \mathbf{e}$, since $\mathbf{c} \in C$, we have the **syndrome** $\mathbf{s} = H \cdot \mathbf{x} = H \cdot \mathbf{e}$, the decode table looks like:
$$
\def\arraystretch{1.5}
   \begin{array}{c:c}
   \mathbf{e} & \mathbf{s} \\ \hline
   00000 & 000 \\ \hdashline
   00001 & 001 \\ \hdashline
   00010 & 010 \\ \hdashline
   10000 & 011 \\ \hdashline
   00100 & 100 \\ \hdashline
   01000 & 101 \\ \hdashline
   00110 & 110 \\ \hdashline
   10100 & 111 \\ \hdashline 
\end{array}
$$
<br />

We can obtain the error term or representative of coset, $\mathbf{e}$,after looking up the above decode table, treating syndrome $s$ as query.

<br />

Finally we can get the orginal codeword $\mathbf{c}$:
$$
\begin{aligned}
\mathbf{c} = \mathbf{x} - \mathbf{e} \\
= \mathbf{x} + \mathbf{e}
\end{aligned}
$$

<br />

#### References

[1] [Abstract Algebra: Theory and Applications](http://abstract.ups.edu/)