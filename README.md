# crypto_research

crypto research stuffs from numerics to ZKP related protocols all written with dirty hands.

<br />

## Code Structure
- abstract_algebraic [Ongoing]
    - group

- ecc 
    - integer_arithmetic [Done]

        specificially for basic arithmetics on big integers

        - add/substruction/multiplication/division

        - euclid extended gcd/lehmer extended gcd

    - finite_field_arithmetic [Done]

        specificailly for basic arithmetics on finite field

        - add/substruction/multiplication/division/inversion/modulo/exponentiation/sqrt/square

    - elliptic_curve_arithmetic [Done]

        specially for neccessary arithmetics on elliptic curves

        - add/doubling/scalar_mul/...

    - hyperelliptic_curve_arithmetic [TODO] 
    - special_curve_arithmetic [TODO]
    - pairings [TODO]
    - ...

- hash
    - shake128(variable output length) [Done]
    - poseidon [TODO]

- pcs
    - sparse_polynomial [Done]
    - ipa [Done]
    - pb [Ongoing]

- recursive snark
    ...

<br />

## Notes 

- [Integer Arithmetic](https://hackmd.io/@70xfCGp1QViTYYJh3AMrQg/rkF-5hHwT)
- [Foundations of Elliptic Curves
](https://hackmd.io/@70xfCGp1QViTYYJh3AMrQg/HJ7rcsY4a)

<br />

## Credits

[1] Handbook of Elliptic and Hyperelliptic Curve Cryptography

[2] Guide to Elliptic Curve Cryptography

[3] Pairings For Beginners

[4] [Algorithms for Modern Hardware](https://en.algorithmica.org/hpc/)

[5] [IPA PCS](https://hackmd.io/@arijitdutta67/r1ZFKoHy2#Accumulation-of-IPA-PCS-and-Recursive-Process-in-Aztec-3)