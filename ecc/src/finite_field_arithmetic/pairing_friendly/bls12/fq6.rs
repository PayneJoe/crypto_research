use crate::finite_field_arithmetic::bigint::BigInt;
use crate::finite_field_arithmetic::pairing_friendly::bls12::{fq::Fq, fq2::Fq2};
use crate::finite_field_arithmetic::pairing_friendly::cubic_extension::{
    CubicExtension, CubicExtensionConfig,
};

const NUM_LIMBS: usize = 6;

#[derive(Copy, Clone)]
pub struct Fq6Config;

impl CubicExtensionConfig<NUM_LIMBS> for Fq6Config {
    type BasePrimeField = Fq<NUM_LIMBS>;
    type BaseField = Fq2;
    type FrobCoeff = Fq2;

    // absolute degree of current (extension) field
    const DEGREE_OVER_BASE_PRIME_FIELD: usize = 6;

    // Fq12 = Fq6[w] / X^2 - gamma, where gamma = v
    // Fq6 = Fq2[v] / X^3 - beta, where beta = u + 1
    // Fq2 = Fq[u] / X^2 - alpha, where alpha = -1
    // beta = Fq2([1, 1])
    const NON_CUBIC_RESIDUAL: Self::BaseField = Fq2 {
        c0: Fq(BigInt([
            8505329371266088957,
            17002214543764226050,
            6865905132761471162,
            8632934651105793861,
            6631298214892334189,
            1582556514881692819,
        ])),
        c1: Fq(BigInt([
            8505329371266088957,
            17002214543764226050,
            6865905132761471162,
            8632934651105793861,
            6631298214892334189,
            1582556514881692819,
        ])),
    };

    // beta^{(p^d - 1)/3} for c1
    const FROBENIUS_COEFF_C1: &'static [Self::FrobCoeff] = &[
        Fq2 {
            c0: Fq(BigInt([
                8505329371266088957,
                17002214543764226050,
                6865905132761471162,
                8632934651105793861,
                6631298214892334189,
                1582556514881692819,
            ])),
            c1: Fq(BigInt([0, 0, 0, 0, 0, 0])),
        },
        Fq2 {
            c0: Fq(BigInt([0, 0, 0, 0, 0, 0])),
            c1: Fq(BigInt([
                14772873186050699377,
                6749526151121446354,
                6372666795664677781,
                10283423008382700446,
                286397964926079186,
                1796971870900422465,
            ])),
        },
        Fq2 {
            c0: Fq(BigInt([
                3526659474838938856,
                17562030475567847978,
                1632777218702014455,
                14009062335050482331,
                3906511377122991214,
                368068849512964448,
            ])),
            c1: Fq(BigInt([0, 0, 0, 0, 0, 0])),
        },
        Fq2 {
            c0: Fq(BigInt([0, 0, 0, 0, 0, 0])),
            c1: Fq(BigInt([
                8505329371266088957,
                17002214543764226050,
                6865905132761471162,
                8632934651105793861,
                6631298214892334189,
                1582556514881692819,
            ])),
        },
        Fq2 {
            c0: Fq(BigInt([
                14772873186050699377,
                6749526151121446354,
                6372666795664677781,
                10283423008382700446,
                286397964926079186,
                1796971870900422465,
            ])),
            c1: Fq(BigInt([0, 0, 0, 0, 0, 0])),
        },
        Fq2 {
            c0: Fq(BigInt([0, 0, 0, 0, 0, 0])),
            c1: Fq(BigInt([
                3526659474838938856,
                17562030475567847978,
                1632777218702014455,
                14009062335050482331,
                3906511377122991214,
                368068849512964448,
            ])),
        },
    ];
    // beta^{2 * (p^d - 1)/3} for c2
    const FROBENIUS_COEFF_C2: &'static [Self::FrobCoeff] = &[
        Fq2 {
            c0: Fq(BigInt([
                8505329371266088957,
                17002214543764226050,
                6865905132761471162,
                8632934651105793861,
                6631298214892334189,
                1582556514881692819,
            ])),
            c1: Fq(BigInt([0, 0, 0, 0, 0, 0])),
        },
        Fq2 {
            c0: Fq(BigInt([
                9875771541238924739,
                3094855109658912213,
                5802897354862067244,
                11677019699073781796,
                1505592401347711080,
                1505729768134575418,
            ])),
            c1: Fq(BigInt([0, 0, 0, 0, 0, 0])),
        },
        Fq2 {
            c0: Fq(BigInt([
                14772873186050699377,
                6749526151121446354,
                6372666795664677781,
                10283423008382700446,
                286397964926079186,
                1796971870900422465,
            ])),
            c1: Fq(BigInt([0, 0, 0, 0, 0, 0])),
        },
        Fq2 {
            c0: Fq(BigInt([
                4897101644811774638,
                3654671041462534141,
                569769440802610537,
                17053147383018470266,
                17227549637287919721,
                291242102765847046,
            ])),
            c1: Fq(BigInt([0, 0, 0, 0, 0, 0])),
        },
        Fq2 {
            c0: Fq(BigInt([
                3526659474838938856,
                17562030475567847978,
                1632777218702014455,
                14009062335050482331,
                3906511377122991214,
                368068849512964448,
            ])),
            c1: Fq(BigInt([0, 0, 0, 0, 0, 0])),
        },
        Fq2 {
            c0: Fq(BigInt([
                17076301903736715834,
                13907359434105313836,
                1063007777899403918,
                15402659025741563681,
                5125705813544623108,
                76826746747117401,
            ])),
            c1: Fq(BigInt([0, 0, 0, 0, 0, 0])),
        },
    ];

    fn multiply_frobenius_coeff(c1: &mut Self::BaseField, c2: &mut Self::BaseField, power: usize) {
        *c1 = *c1 * Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
        *c2 = *c2 * Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
    }
}

pub type Fq6 = CubicExtension<NUM_LIMBS, Fq6Config>;
impl Fq6 {}

mod tests {
    use super::*;
    #[test]
    fn test_from() {}

    #[test]
    fn test_addition() {}

    #[test]
    fn test_substraction() {}

    #[test]
    fn test_multiplication() {}

    #[test]
    fn test_square() {}

    #[test]
    fn test_square_root() {}

    #[test]
    fn test_inverse() {}
}
