# crypto_research

crypto research stuffs from numeric arithmetics to ZK-applied protocols all written with dirty hands.

<br />

## Code Structure

- docs
    - bls12_381
        note and python versioned bls12-381 implementation [$`\color{green}\checkmark`$]

    - on_proving_pairings 
        note and python versioned bn254 implementation [$`\color{green}\checkmark`$]

- ecc 
    - integer_arithmetic 

        - basic arithmetics on big integers

          - add/substruction/multiplication/division [$`\color{green}\checkmark`$]

          - euclid extended gcd/lehmer extended gcd [$`\color{green}\checkmark`$]

    - finite_field_arithmetic 

        - basic arithmetics on base field $F_q$ and its instantiation

          - add/substruction/multiplication/division/inversion/modulo/exponentiation/sqrt/square [$`\color{green}\checkmark`$]
          - field implementation for pallas/vasta curves [$`\color{green}\checkmark`$] 

        - basic arithmetics on extension field $F_{q^k}$ of $F_q$

            - quadratic extension $F_{q^2}/F_{q}$ [$`\color{green}\checkmark`$]

            - cubic extension $F_{q^3}/F_{q}$ [$`\color{green}\checkmark`$]

            - cyclotomic [$`\color{green}\checkmark`$]

            - twist/untwist $\Phi: F_{q^k} \mapsto F_{q^{k / d}}$ [$`\color{green}\checkmark`$]

            - frobenius map [$`\color{green}\checkmark`$]
            - field implementation for BLS12/MNT/BN pairing-family curves [$`\color{green}\checkmark`$]

    - elliptic_curve_arithmetic 

        - neccessary arithmetics on elliptic curves over base field $F_q$ 

          - add/doubling/scalar_mul/... [$`\color{green} \checkmark`$]

        - neccessary arithmetics on pairing-friendly elliptic curves over extension field $F_{q^k}$ and its instantiation

          - add/doubling/scalar_mul/... [$`\color{green} \checkmark`$]

          - BLS12/MNT/BN pairing-friendly curves [$`\color{green} \checkmark`$]

    - hyperelliptic_curve_arithmetic [TODO] 

    - special_curve_arithmetic [TODO]

    - pairings 
        -  Weil Pairing [$`\color{green}\checkmark`$]

        - Tate Pairing [$`\color{green}\checkmark`$]

        - Ate Pairing [Ongoing]

        - Optimal Ate Pairings [Ongoing]

    - ...

- hash
    - shake128(variable output length) [$`\color{green}\checkmark`$]
    - poseidon [TODO]

- pcs
    - sparse_polynomial [$`\color{green}\checkmark`$]
    - IPA [$`\color{green} \checkmark`$]
    - KZG [Ongoing]

- recursive snark
    ...

<br />

## Notes 

#### Small Fields

- [Circle Stark from A New Starker
](https://hackmd.io/@70xfCGp1QViTYYJh3AMrQg/SJrem_VKC)
- [Binius](https://hackmd.io/@70xfCGp1QViTYYJh3AMrQg/r1iOQbaFR)

#### Pairings

- [Theory and Practical Implementation of BLS12-381](https://hackmd.io/@70xfCGp1QViTYYJh3AMrQg/ryo55eEeC)
- [On Proving Pairings](https://hackmd.io/@70xfCGp1QViTYYJh3AMrQg/S1cU7YJGC)

#### Proving System
- [Plonk From Scratch](https://hackmd.io/@70xfCGp1QViTYYJh3AMrQg/HJzwPUU7a)
- [Relaxed Plonk Step by Step
](https://hackmd.io/@70xfCGp1QViTYYJh3AMrQg/BkT0ayKmT)
- [Cross Term in Plonk](https://hackmd.io/@70xfCGp1QViTYYJh3AMrQg/SkDf2nIzp)

<br />

## Credits

[1] [Handbook of Elliptic and Hyperelliptic Curve Cryptography](https://blkcipher.pl/assets/pdfs/Handbook_of_Elliptic_and_Hyperelliptic_Curve_Cryptography.pdf)

[2] [Guide to Elliptic Curve Cryptography](http://tomlr.free.fr/Math%E9matiques/Math%20Complete/Cryptography/Guide%20to%20Elliptic%20Curve%20Cryptography%20-%20D.%20Hankerson,%20A.%20Menezes,%20S.%20Vanstone.pdf)

[3] [Pairings For Beginners](https://static1.squarespace.com/static/5fdbb09f31d71c1227082339/t/5ff394720493bd28278889c6/1609798774687/PairingsForBeginners.pdf)

[4] [Algorithms for Modern Hardware](https://en.algorithmica.org/hpc/)

[5] [IPA PCS](https://hackmd.io/@arijitdutta67/r1ZFKoHy2#Accumulation-of-IPA-PCS-and-Recursive-Process-in-Aztec-3)

[6] [Guide to Pairing-based Cryptography](https://www.math.u-bordeaux.fr/~damienrobert/csi/book/book.pdf)

[7] [A note on twists for pairing friendly curves](http://indigo.ie/%7Emscott/twists.pdf)

[8] [Pairing-Friendly Elliptic Curves of Prime Order](https://eprint.iacr.org/2005/133.pdf)

[9] [Optimal Pairings](https://eprint.iacr.org/2008/096.pdf)

[10] [On the Computation of the Optimal Ate Pairing at the 192-bit Security Level](https://eprint.iacr.org/2016/130.pdf)

[11] [Efficient Final Exponentiation via Cyclotomic Structure for Pairings over Families of Elliptic
Curves](https://eprint.iacr.org/2020/875.pdf)

[12] [Faster Squaring in the Cyclotomic Subgroup of Sixth Degree Extensions](https://eprint.iacr.org/2009/565.pdf)

[13] [A Guide to Plane Algebraic Curves]()

[14] [Circle Stark](https://eprint.iacr.org/2024/278.pdf)

[15] [Binius](https://eprint.iacr.org/2023/1784.pdf)