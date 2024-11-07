# Permute Basis

Based on [https://doi.org/10.1063/1.3518900]

|            Basis             | TrAnslation | Parity | Spin-Flip |      |
| :--------------------------: | :---------: | :----: | :-------: | :--: |
|       `ProjectedBasis`       |      ✗      |   ✗    |     ✗     |      |
|     `TranslationalBasis`     |      ✓      |   ✗    |     ✗     |      |
|        `ParityBasis`         |      ✗      |   ✓    |     ✗     |      |
|         `FlipBasis`          |      ✗      |   ✗    |     ✓     |      |
|   `TranslationParityBasis`   |      ✓      |   ✓    |     ✗     |      |
|    `TranslationFlipBasis`    |      ✓      |   ✗    |     ✓     |      |
|      `ParityFlipBasis`       |      ✗      |   ✓    |     ✓     |      |
| `TranslationParityFlipBasis` |      ✓      |   ✓    |     ✓     |      |



## Translational Symmetry

For a periodic chain, we define the translation operator as moving the spins one step cyclically to the "right";
$$
T\left|S_1^z, \ldots, S_{L}^z\right\rangle=\left|S_{L}^z, S_1^z, \ldots, S_{L-1}^z\right\rangle
$$
We can construct momentum states $|\Psi(k)\rangle$, which by definition are eigenstates of the translation operator,
$$
T|\Psi(k)\rangle=\mathrm{e}^{i k}|\Psi(k)\rangle
$$
Here the allowed momenta are $k=2 n \pi / L$, with $n=0, \ldots, L-1$, following from the fact that $T^L=1$. States with different $k$ form their own individually diagonalizable blocks of the Hamiltonian. 

### Normalization

A momentum state can be constructed using a reference state $|a\rangle$ (a single state in the $z$-component basis) and all its translations:
$$
|a(k)\rangle=\frac{1}{R_a} \sum_{r=0}^{L-1} \mathrm{e}^{-i k r} T^r|a\rangle .
$$
It can easily be verified that  $T|a(k)\rangle=\mathrm{e}^{i k}|a(k)\rangle$, which is the definition of a momentum state. If all the translated states $T^r|a\rangle$ are distinct, the normalization constant is just $\sqrt{L}$. Some reference states have periodicities less than $L$, and this affects normalization. The periodicity of a state is defined as the smallest integer $N_a$ for which
$$
T^{N_a}|a\rangle=|a\rangle, \quad N_a \in\{1, \ldots, L\} .
$$
If $N_a < L$, then there are multiple copies of the same state in the sum, and the normalization constant must be modified accordingly. 

> **Claim:** For given $|a\rangle$, the allowed momenta are those for which $k N_a = 2\pi m$. For the allowed momenta, the normalization constant is $R_a = L/\sqrt{N_a}$.
>
> **Proof.** The periodicity of the representative has to be compatible with the momentum in order to be a viable state. The compatibility is related to normalizability. The sum of phase factors associated with the representative state $|a\rangle$ is
> $$
> F\left(k, N_a\right)=\sum_{n=0}^{L/N_a-1} \mathrm{e}^{-i k n N_a}
> = \begin{cases}
> L / N_a, & k N_a = 2 \pi m \\ 
> 0, & \text { otherwise }
> \end{cases}.
> $$
> The normalization constant is then
> $$
> R_a^2=\langle a(k) | a(k)\rangle=N_a\left|F\left(k, N_a\right)\right|^2.
> $$
> Therefore, if $F\left(k, N_a\right)=0$, no state with momentum $k$ can be defined using the reference state $|a\rangle$. Thus,  
> $$
> k=\frac{2 \pi}{R_a} m, \quad
> R_a = \frac{L}{\sqrt{N_a}}.
> $$

We remark that the summation in the definition of $|a(k)\rangle$ can also be written as
$$
|a(k)\rangle = \frac{1}{\sqrt{N_a}}\sum_{r=0}^{N_a-1} e^{-ikr}T^r|a\rangle.
$$
The coefficient is related to normalization $N_a^{-1/2} = R_a/L$.

### Matrix elements

Consider the translational invariant Hamiltonian $H=\sum_{j=1}^L h_j$. We need to find the state resulting when $H$ acts on the momentum states. Since $[H, T]=0$ we can write
$$
H|a(k)\rangle = \frac{1}{R_a} \sum_{r=0}^{L-1} \mathrm{e}^{-i k r} T^r H|a\rangle
= \frac{1}{R_a} \sum_{j=1}^L \sum_{r=0}^{L-1} \mathrm{e}^{-i k r} T^r h_j|a\rangle.
$$
We need to operate with the Hamiltonian operators $h_j$ only on the reference state. We can write $h_j|a\rangle= \sum_{b'}h_j(b',a)\left|b^{\prime}\right\rangle$. The prime in $\left|b^{\prime}\right\rangle$ is there to indicate that this new state is not necessarily one of the reference states used to define the basis and, therefore, a momentum state should not be written directly based on it. Provided that $\left|b^{\prime}\right\rangle$ is compatible with the momentum, there must be a reference state $\left|b\right\rangle$ which is related to it by 
$$
\left|b\right\rangle=T^{l(b')}\left|b^{\prime}\right\rangle.
$$
Using this relation we have
$$
H|a(k)\rangle=\sum_{j,b'} \frac{h_j(b',a)}{R_a} \sum_{r=0}^{L-1} \mathrm{e}^{-i k r} T^{r-l(b')}\left|b_j\right\rangle
=\sum_{j,b'} h_j(b',a) \mathrm{e}^{-i k l(b')} \frac{R_{b}}{R_a}\left|b(k)\right\rangle.
$$
We thus obtain the matrix element
$$
\langle a(k)|H|b(k)\rangle = \sum_{j=1}^L\sum_{b'} h_j(b',a) \mathrm{e}^{-i k l(b')} \frac{R_{b}}{R_a}.
$$

## Reflection and Spin-Flip Symmetry

The spatial reflection operator is defined as
$$
P\left|S_1^z, \ldots, S_{L}^z\right\rangle=\left|S_{L}^z, \ldots, S_1^z\right\rangle
$$
For an eigenstate of $P|\Psi(p)\rangle=p|\Psi(p)\rangle$, where $p= \pm 1$ since $P^2=1$. We will use $T$ and $P$ for block-diagonalization, although they **cannot always be used simultaneously** because $[T, P]=0$ only in a sub-space of the Hilbert space. For a system with open boundaries, $T$ is not defined, but $P$ can be used.

For the special (and most important) case $m_z=0$ (for even $N$ ), we can block-diagonalize using a discrete subset of all the possible rotations in spin-space; the spin-flip symmetry, i.e., invariance with respect to flipping all the spins. This is defined formally by an operator we call $Z$;
$$
Z\left|S_1^z, \ldots, S_{L}^z\right\rangle=\left|-S_0^z,-S_1^z, \ldots,-S_{L}^z\right\rangle
$$
For this operator we again have $Z^2=1$ and the eigenvalues $z= \pm 1$. Since $Z$ commutes with both $P$ and $T$, it can be used together with these operators to further block-diagonalize $H$.

### States with Parity

Consider the following extension of the momentum state:
$$
|a(k, p)\rangle=\frac{1}{R_a} \sum_{r=0}^{N-1}\sum_{s=0}^1 \mathrm{e}^{-i k r}p^s T^r P^s|a\rangle,
$$
where $p= \pm 1$. Clearly, this is a state with momentum $k$, but is it also an eigenstate of $P$ with parity $p$. We can check this by explicit operation with $P$, using $P^2=1, p^2=1$, and the relationship $P T=T^{-1} P$ :
$$
P|a(k, p)\rangle =\frac{1}{R_a} \sum_{r=0}^{N-1} \mathrm{e}^{-i k r} T^{-r}(P+p)|a\rangle
 = \frac{p}{R_a} \sum_{r=0}^{N-1} \mathrm{e}^{i k r} T^r(1+p P)|a\rangle .
$$
This is not exactly of the original form unless $k=0$ or $\pi$, for which $\mathrm{e}^{i k r}=\mathrm{e}^{-i k r}$. Thus, in these two special cases, parity and translational invariance can be used simultaneously for block-diagonalization and $|a(k, p)\rangle$ is indeed a momentum state with parity $p$ (or, in other words, $[T, P]=0$ in the sub-spaces with momenta $k=0$ and $\pi$ ).

### General Abelian Symmetry

The Abelian symmetry $G=Z_L\times Z_2$

