# crypto_research

crypto research stuffs from numeric arithmetics to ZK-applied protocols all written with dirty hands.

<br />

## Code Structure
- ecc 
    - integer_arithmetic 

        - basic arithmetics on big integers

          - add/substruction/multiplication/division [\checkmark]

          - euclid extended gcd/lehmer extended gcd [\checkmark]

    - finite_field_arithmetic 

        - basic arithmetics on base field $F_q$

          - add/substruction/multiplication/division/inversion/modulo/exponentiation/sqrt/square [\checkmark]

        - basic arithmetics on extension field $F_{q^k}$ of $F_q$

            - quadratic extension [Ongoing]

            - cubic extension [Ongoing]

            - sextic extension [Ongoing]

            - cyclotomic [Ongoing]

            - twist/untwist [Ongoing]

    - elliptic_curve_arithmetic 

        - neccessary arithmetics on elliptic curves over base field $F_q$ 

          - add/doubling/scalar_mul/... [\checkmark]

        - neccessary arithmetics on pairing-friendly elliptic curves over extension field $F_{q^k}$

          - add/doubling/scalar_mul/... [Ongoing]

    - hyperelliptic_curve_arithmetic [TODO] 

    - special_curve_arithmetic [TODO]

    - pairings 

        - Tate Pairing [Ongoing]

        - Ate Pairing [Ongoing]

        - Optimal Ate Pairings [Ongoing]

    - ...

- hash
    - shake128(variable output length) [\checkmark]
    - poseidon [TODO]

- pcs
    - sparse_polynomial [\checkmark]
    - IPA [\checkmark]
    - KZG [Ongoing]

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

[3] Pairings For Beginners [https://static1.squarespace.com/static/5fdbb09f31d71c1227082339/t/5ff394720493bd28278889c6/1609798774687/PairingsForBeginners.pdf]

[4] [Algorithms for Modern Hardware](https://en.algorithmica.org/hpc/)

[5] [IPA PCS](https://hackmd.io/@arijitdutta67/r1ZFKoHy2#Accumulation-of-IPA-PCS-and-Recursive-Process-in-Aztec-3)