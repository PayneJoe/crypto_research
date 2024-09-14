## Core Concepts

#### Linear Code

If a code $\mathcal{C}$ is also a **group**, then we call it *linear code*.


A great property of *linear code* is that:
- the minimum distance of code equals the minimum weight 
    $$
        d(\mathbf{x}, \mathbf{y}) = w(\mathbf{x} - \mathbf{y})
    $$
    since $\mathbf{x}, \mathbf{y} \in \mathcal{C}$ and <span style="color: red">$\mathbf{x} - \mathbf{y} \in \mathcal{C}$</span> ($\mathcal{C}$ is a group), then we must have:
    $$
        d_{min}(\mathcal{C}) = w_{min}(\mathcal{C})
    $$

<br />

#### Dimension of Kernel

there is a map $\varphi: \mathbb{R}^3 \rightarrow \mathbb{R}^2$, for example, given a matrices:
$$
A = 
\begin{pmatrix}
1 & 2 & 3 \\ 
4 & 5 & 6 \\
\end{pmatrix}
$$

and a vector $\mathbf{x} \in \mathbb{R}^3$, we have the kernel representation:
$$
K(\varphi) = \{\mathbf{x} | A \cdot \mathbf{x} = \mathbf{0}\}
$$

equivalently:
$$
\begin{aligned}
x_1 + 2 x_2 + 3 x_3 = 0 \\
4 x_1 + 5 x_2 + 6 x_3 = 0 \\
\end{aligned}
$$
then we have solution $x_1 = -x_3, x_2 = 2 x_3$. Finally we have the representation of $\mathbf{x}$:
$$
\mathbf{x} = 
\begin{pmatrix}
-x \\
2x \\
x \\
\end{pmatrix}
$$

<br />

<span style="color: red">So, kernel $K(\varphi)$ is spaned by a vector $(-1, 2, 1)$, the dimension of kernel is 1. This is a 1-dimensional subspace of $\mathbb{R}^3$.</span>

<br />

#### Rank-Nullity Theorem

$$
Dim(D) = rank(\varphi) + nullity(\varphi)
$$

where:

- <span style="color: red"> nullity($\varphi$) </span>: dimension of kernel $K(\varphi)$, say $dim(\mathbf{x})$, it's $1$ for above example

- <span style="color: red"> rank($\varphi$) </span>: dimension of image $Img(\varphi)$, say $dim(A \cdot \mathbf{x})$, it's $2$ for above example

- <span style="color: red"> dim($D$) </span>: dimension of domain, say $dim(\mathbb{R}^3)$,it's $3$ for above example

<br />

#### Dual Codes

> all defined over binary field $F_2$.

<br />

Generator matrices $G$ of an $[n, k]$ code $\mathcal{C}$:
$$
G = [I_k | A]
= \begin{bmatrix}
1 & 0 & 0 & 0 \mid 0 & 1 & 1 \\
0 & 1 & 0 & 0 \mid 1 & 0 & 1 \\
0 & 0 & 1 & 0 \mid 1 & 1 & 0 \\
0 & 0 & 0 & 1 \mid 1 & 1 & 1 \\
\end{bmatrix}
$$

<br />

parity check matrces $H$ of code $\mathcal{C}$:
$$
H = [-A^T | I_{n - k}] = [A^T | I_{n - k}]
= \begin{bmatrix}
0 & 1 & 1 & 1 \mid 1 & 0 & 0 \\
1 & 0 & 1 & 1 \mid 0 & 1 & 0 \\
1 & 1 & 0 & 1 \mid 0 & 0 & 1 \\
\end{bmatrix}
$$

<br />

message $\mathbf{m} = \begin{pmatrix}x_0 \\ x_1 \\ x_2 \\ x_3 \end{pmatrix}$, codeword $\mathbf{x}$ of code $\mathcal{C}$ looks like:
$$
\mathbf{x} = G^T \cdot \mathbf{m} =
\begin{pmatrix}
x_0 \\
x_1 \\
x_2 \\
x_3 \\
x_1 + x_2 + x_3 \\
x_0 + x_2 + x_3 \\
x_0 + x_1 + x_3 \\
\end{pmatrix}
$$

<br />

$H$ is also the generator matrix of some code $\mathcal{C}^{\perp}$, message $\mathbf{m'} = \begin{pmatrix}x_0' \\ x_1' \\ x_2' \end{pmatrix}$, then codeword $\mathbf{x'}$ of code $\mathcal{C}^{\perp}$ looks like:
$$
\mathbf{x'} = 
\begin{pmatrix}
x_1' + x_2' \\
x_0' + x_2' \\
x_0' + x_1' \\
x_0' + x_1' + x_2' \\
x_0' \\
x_1' \\
x_2' \\
\end{pmatrix}
$$

Two conclusions:
- <span style="color: red">$G$ and $H$ are generator and parity check matrices of code $\mathcal{C}$, at the same time, $H$ and $G$ are also generator and parity check matrices of code $\mathcal{C}^{\perp}$</span>
    $$
    \begin{aligned}
    H \cdot (G^T \cdot \mathbf{m}) = H \cdot \mathbf{x} = \mathbf{0} \\
    G \cdot (H^T \cdot \mathbf{m'}) = G \cdot \mathbf{x'}= \mathbf{0} \\
    \end{aligned}
    $$

- <span style="color: red">codeword $\mathbf{x}$ is conjugative with codeword $\mathbf{x'}$, say $\langle \mathbf{x} \cdot \mathbf{x'} \rangle = 0$</span>
    $$
    \mathcal{C}^{\perp} = \{\mathbf{x} \in \mathbb{F}_q^n | \mathbf{x} \cdot \mathbf{c} = 0, \mathbf{c} \in \mathcal{C} \}
    $$

- <span style="color: red">If $\mathcal{C} = \mathcal{C}^{\perp}$, then we can say $\mathcal{C}$ is *self-dual* code</span>

<br />

What is *self-dual* anyway?

<br />

Let's take the $[8, 4]$ code $\mathcal{C}$ as example for illustration:
$$
G = 
\begin{bmatrix}
1 & 0 & 0 & 0 \mid 0 & 1 & 1 & 1 \\
0 & 1 & 0 & 0 \mid 1 & 0 & 1 & 1 \\
0 & 0 & 1 & 0 \mid 1 & 1 & 0 & 1 \\
0 & 0 & 0 & 1 \mid 1 & 1 & 1 & 0 \\
\end{bmatrix},
H = 
\begin{bmatrix}
0 & 1 & 1 & 1 \mid 1 & 0 & 0 & 0\\
1 & 0 & 1 & 1 \mid 0 & 1 & 0 & 0\\
1 & 1 & 0 & 1 \mid 0 & 0 & 1 & 0\\
1 & 1 & 1 & 0 \mid 0 & 0 & 0 & 1\\
\end{bmatrix}
$$

<br />

Then we have codeword $\mathbf{x}$ of code $\mathcal{C}$ and codeword $\mathbf{x'}$ of code $\mathcal{C}^{\perp}$:
$$
\mathbf{x} = \begin{pmatrix}
x_0 \\
x_1 \\
x_2 \\
x_3 \\
x_1 + x_2 + x_3 \\
x_0 + x_2 + x_3 \\
x_0 + x_1 + x_3 \\
x_0 + x_1 + x_2 \\
\end{pmatrix},
\mathbf{x'} = \begin{pmatrix}
x_1 + x_2 + x_3 \\
x_0 + x_2 + x_3 \\
x_0 + x_1 + x_3 \\
x_0 + x_1 + x_2 \\
x_0 \\
x_1 \\
x_2 \\
x_3 \\
\end{pmatrix}
$$

<br />

then we can obtain the table of code $\mathcal{C}$:
$$
\def\arraystretch{1.5}
   \begin{array}{c:c:c}
   \mathbf{m} & \mathbf{x} \\ \hline
   0000 & 00000000 \\ \hdashline
   1000 & 10000111 \\ \hdashline
   0100 & 01001011 \\ \hdashline
   1100 & 11001100 \\ \hdashline
   0010 & 00101101 \\ \hdashline
   1010 & 10101010 \\ \hdashline
   0110 & 01100110 \\ \hdashline
   1110 & 11100001 \\ \hdashline
   0001 & 00011110 \\ \hdashline
   1001 & 10011001 \\ \hdashline
   0101 & 01010101 \\ \hdashline
   1101 & 11010010 \\ \hdashline
   0011 & 00110011 \\ \hdashline
   1011 & 10110100 \\ \hdashline
   0111 & 01111000 \\ \hdashline
   1111 & 11111111 \\ \hdashline
   \end{array}
$$

and the table of code $\mathcal{C}^{\perp}$:
$$
\def\arraystretch{1.5}
   \begin{array}{c:c:c}
   \mathbf{m} & \mathbf{x'} \\ \hline
   0000 & 00000000 \\ \hdashline
   1000 & 01111000 \\ \hdashline
   0100 & 10110100 \\ \hdashline
   1100 & 11001100 \\ \hdashline
   0010 & 11010010 \\ \hdashline
   1010 & 10101010 \\ \hdashline
   0110 & 01100110 \\ \hdashline
   1110 & 00011110 \\ \hdashline
   0001 & 11100001 \\ \hdashline
   1001 & 10011001 \\ \hdashline
   0101 & 01010101 \\ \hdashline
   1101 & 00101101 \\ \hdashline
   0011 & 00110011 \\ \hdashline
   1011 & 01001011 \\ \hdashline
   0111 & 10000111 \\ \hdashline
   1111 & 11111111 \\ \hdashline
   \end{array}
$$

<br />

<span style="color: red">We can say code $\mathcal{C}$ is *self-dual* since $\mathcal{C} = \mathcal{C}^{\perp}$</span>.

<br />

Let's take a closer look at the $[8, 4]$ code $\mathcal{C}$, what properties can we find:
1. length of code is a even number, $|\mathcal{C}| = n = 2^4 = 16$
2. dimension of code $dim(\mathcal{C}) = n / 2= 8$
3. minimum distance $d_{min}(\mathcal{C}) = 4$
4. inner product of any two codewords is $0$, $\langle \mathbf{x}, \mathbf{y} \rangle = 0, \mathbf{x}, \mathbf{y} \in \mathcal{C}$

    $$
    \langle \mathbf{x}, \mathbf{y} \rangle = \sum_i x_i \cdot y_i \mod 2 = 0
    $$

From property of $\boxed{4}$, <span style="color: red">we can also conclude that code $\mathcal{C}$ is also *self-orthogonal*</span>.

<br />

#### Self-Orthogonal Binary Code

Let's dig into self-orthogonal a little more, and see what properties does it have.

<br />

Another example of $[7, 4]$ code $\mathcal{C}$, with parity check matrix:
$$
H = 
\begin{bmatrix}
1 & 1 & 0 & 1 \mid 1 & 0 & 0 \\ 
0 & 1 & 1 & 1 \mid 0 & 1 & 0 \\
1 & 1 & 1 & 0 \mid 0 & 0 & 1 \\
\end{bmatrix}
$$
$H$ is also generator matrix of $[7, 3]$ code $\mathcal{C}^{\perp}$, then the code $\mathcal{C}^{\perp}$ would be like:
$$
\def\arraystretch{1.5}
   \begin{array}{c:c:c}
   \mathbf{m} & \mathbf{x'} \\ \hline
   000 & 0000000 \\ \hdashline
   100 & 1101100 \\ \hdashline
   010 & 0111010 \\ \hdashline
   110 & 1010110 \\ \hdashline
   001 & 1110001 \\ \hdashline
   101 & 0011101 \\ \hdashline
   011 & 1001011 \\ \hdashline
   111 & 0100111 \\ \hdashline
   \end{array}
$$

<br />

Generator matrix of $[7, 4]$ code $\mathcal{C}$ is:
$$
G = 
\begin{bmatrix}
1 & 0 & 0 & 0 \mid 1 & 0 & 1 \\
0 & 1 & 0 & 0 \mid 1 & 1 & 1 \\
0 & 0 & 1 & 0 \mid 0 & 1 & 1 \\
0 & 0 & 0 & 1 \mid 1 & 1 & 0 \\
\end{bmatrix}
$$
then the code $\mathcal{C}$ would be like:
$$
\def\arraystretch{1.5}
   \begin{array}{c:c:c}
   \mathbf{m} & \mathbf{x} \\ \hline
   0000 & 0000000 \\ \hdashline
   1000 & 1000101 \\ \hdashline
   0100 & 0100111 \\ \hdashline
   1100 & 1100010 \\ \hdashline
   0010 & 0010011 \\ \hdashline
   1010 & 1010110 \\ \hdashline
   0110 & 0110100 \\ \hdashline
   1110 & 1110001 \\ \hdashline
   0001 & 0001110 \\ \hdashline
   1001 & 1001011 \\ \hdashline
   0101 & 0101001 \\ \hdashline
   1101 & 1101100 \\ \hdashline
   0011 & 0011101 \\ \hdashline
   1011 & 1011000 \\ \hdashline
   0111 & 0111010 \\ \hdashline
   1111 & 1111111 \\ \hdashline
   \end{array}
$$

We call code $\mathcal{C}^{\perp}$ is *self-orthogonal* since $\mathcal{C}^{\perp} \subseteq \mathcal{C}$.

<br />

## Theorems

#### Minimum Distance VS Partity Check Matrix

<span style="color: red">The minimum distance of a linear code $\mathcal{C}$ equals the least number of independent columns of parity check matrix</span>.

<br />

For example, a parity check matrix:
$$
H = 
\begin{bmatrix}
0 & 0 & 0 & 1 & 1 & 1 & 1 \\
0 & 1 & 1 & 0 & 0 & 1 & 1 \\
1 & 0 & 1 & 0 & 1 & 0 & 1 \\
\end{bmatrix}
$$

The minimum independent columns is $3$, such as column $1, 2, 3$, sum of these three columns $[0, 0, 1] + [0, 1, 0] + [0, 1, 1] = [0, 0, 0]$. So the minimum distance of the cooresponding code to $H$ is $3$.

<br />

## References

[1] [Fundamentals of Error Correcting Code](https://theswissbay.ch/pdf/Gentoomen%20Library/Information%20Theory/Coding%20Theory/Fundamentals%20of%20Error-Correcting%20Codes%20-%20W.%20Cary%20Huffman.pdf)