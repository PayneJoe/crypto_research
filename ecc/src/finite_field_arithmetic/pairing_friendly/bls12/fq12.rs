use crate::finite_field_arithmetic::bigint::BigInt;
use crate::finite_field_arithmetic::pairing_friendly::bls12::{
    fq::Fq, fq2::Fq2, fq6::Fq6, fq6::Fq6Config,
};
use crate::finite_field_arithmetic::pairing_friendly::cubic_extension::CubicExtensionConfig;
use crate::finite_field_arithmetic::pairing_friendly::cyclotomic::CyclotomicGroup;
use crate::finite_field_arithmetic::pairing_friendly::field::Field;
use crate::finite_field_arithmetic::pairing_friendly::quadratic_extension::{
    QuadraticExtension, QuadraticExtensionConfig,
};

const NUM_LIMBS: usize = 6;
const MAX_SCALAR_LIMBS: usize = 8;

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct Fq12Config;

impl QuadraticExtensionConfig<NUM_LIMBS> for Fq12Config {
    type BasePrimeField = Fq<NUM_LIMBS>;
    type BaseField = Fq6;
    type FrobCoeff = Fq2;

    // absolute degree of current (extension) field
    const DEGREE_OVER_BASE_PRIME_FIELD: usize = 12;

    // Fq12 = Fq6[w] / X^2 - gamma, where gamma = v
    // Fq6 = Fq2[v] / X^3 - beta, where beta = u + 1
    // Fq2 = Fq[u] / X^2 - alpha, where alpha = -1
    // gamma = Fq6([0, 1, 0])
    const NON_QUADRATIC_RESIDUAL: Self::BaseField = Fq6 {
        c0: Fq2 {
            c0: Fq(BigInt([0, 0, 0, 0, 0, 0])),
            c1: Fq(BigInt([0, 0, 0, 0, 0, 0])),
        },
        c1: Fq2 {
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
        c2: Fq2 {
            c0: Fq(BigInt([0, 0, 0, 0, 0, 0])),
            c1: Fq(BigInt([0, 0, 0, 0, 0, 0])),
        },
    };

    // precomputed constant coefficients of frobenius map on quadratic extension
    // namely alpha^{(p^d - 1)/2}
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
            c0: Fq(BigInt([
                506819140503852133,
                14297063575771579155,
                10946065744702939791,
                11771194236670323182,
                2081670087578406477,
                644615147456521963,
            ])),
            c1: Fq(BigInt([
                12895611875574011462,
                6359822009455181036,
                14936352902570693524,
                13914887797453940944,
                3330433690892295817,
                1229183470191017903,
            ])),
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
        Fq2 {
            c0: Fq(BigInt([
                4480897313486445265,
                4797496051193971075,
                4046559893315008306,
                10569151167044009496,
                2123814803385151673,
                852749317591686856,
            ])),
            c1: Fq(BigInt([
                8921533702591418330,
                15859389534032789116,
                3389114680249073393,
                15116930867080254631,
                3288288975085550621,
                1021049300055853010,
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
            c0: Fq(BigInt([
                3974078172982593132,
                8947176549131943536,
                11547238222321620130,
                17244701004083237929,
                42144715806745195,
                208134170135164893,
            ])),
            c1: Fq(BigInt([
                9428352843095270463,
                11709709036094816655,
                14335180424952013185,
                8441381030041026197,
                5369959062663957099,
                1665664447512374973,
            ])),
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
                12895611875574011462,
                6359822009455181036,
                14936352902570693524,
                13914887797453940944,
                3330433690892295817,
                1229183470191017903,
            ])),
            c1: Fq(BigInt([
                506819140503852133,
                14297063575771579155,
                10946065744702939791,
                11771194236670323182,
                2081670087578406477,
                644615147456521963,
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
            c0: Fq(BigInt([
                8921533702591418330,
                15859389534032789116,
                3389114680249073393,
                15116930867080254631,
                3288288975085550621,
                1021049300055853010,
            ])),
            c1: Fq(BigInt([
                4480897313486445265,
                4797496051193971075,
                4046559893315008306,
                10569151167044009496,
                2123814803385151673,
                852749317591686856,
            ])),
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
                9428352843095270463,
                11709709036094816655,
                14335180424952013185,
                8441381030041026197,
                5369959062663957099,
                1665664447512374973,
            ])),
            c1: Fq(BigInt([
                3974078172982593132,
                8947176549131943536,
                11547238222321620130,
                17244701004083237929,
                42144715806745195,
                208134170135164893,
            ])),
        },
    ];

    // Fp6 * Fp2
    fn multiply_frobenius_coeff(c: &mut Self::BaseField, power: usize) {
        c.c0 = c.c0 * Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
        c.c1 = c.c1 * Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
        c.c2 = c.c2 * Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
    }
}

pub type Fq12 = QuadraticExtension<NUM_LIMBS, Fq12Config>;

impl CyclotomicGroup<NUM_LIMBS> for Fq12 {
    // since Fq12 is quadratic extension of Fq6, so inversion would be easily obtained with conjugation
    //
    fn cyclotomic_inverse(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }
        Some(Self {
            c0: self.c0,
            c1: -self.c1,
        })
    }

    // referenced from 5.3.2(P140) "Guide to Pairing-based Cryptography"
    //
    // for the convenience of utilization superior property of cyclotomic subgroup,
    // we need to reframe field tower of Fp12 first, 2/3/2 -> 3/2/2,
    // transform Fp2 -> Fp6 -> Fp12 into Fp2 -> Fp4 -> Fp12
    // namely transform Fp12 = Fp6[w] / X^2 - v, Fp6 = Fp2[v] / X^3 - (u + 1)
    //        into      Fp12 = Fp4[y] / X^3 - x, Fp4 = Fp2[x] / X^2 - (u + 1)
    //
    // alpha = (c00 + c01 v + c02 v^2) + (c10 + c11 v + c12 v^2)w
    //       = (c00 + c11 x) + (c02 + c10 x) * y  + (c01 + c12 x) * y^2
    //       = a + b * y + c * y^2
    // where alpha \in Fp12, a, b, c \in Fp4
    //
    // alpha^2 = A + B * y + C * y^2
    //         = (3 * a^2 - 2 * a') + (3 * c^2 * (u + 1) + 2 * b') * y + (3 * b^2 - 2c') * y^2
    //
    // a^2 = ((c00 - c11) * (c00 - (u + 1) * c11) + c00 * c11 + (u + 1) * c00 * c11) + 2 * (c00 * c11) * x
    fn cyclotomic_square(&self) -> Self {
        let (a0, a1) = (self.c0.c0, self.c1.c1);
        let (b0, b1) = (self.c0.c2, self.c1.c0);
        let (c0, c1) = (self.c0.c1, self.c1.c2);

        // Fp6 = Fp2[v] / X^3 - beta, where beta = u + 1
        let beta = Fq6Config::NON_CUBIC_RESIDUAL;

        // Since Fp4 is not implemented, we have manually do it here
        // square of a over Fp4
        let (a_square_0, a_square_1) = {
            let (v0, v3, v2) = (a0 - a1, a0 - beta * a1, a0 * a1);
            ((v0 * v3 + v2) + beta * v2, v2 + v2)
        };
        // square of b over Fp4
        let (b_square_0, b_square_1) = {
            let (v0, v3, v2) = (b0 - b1, b0 - beta * b1, b0 * b1);
            ((v0 * v3 + v2) + beta * v2, v2 + v2)
        };
        // square of c over Fp4
        let (mut c_square_0, mut c_square_1) = {
            let (v0, v3, v2) = (c0 - c1, c0 - beta * c1, c0 * c1);
            ((v0 * v3 + v2) + beta * v2, v2 + v2)
        };

        // A = 3 * a^2 - 2 * a' = a^2 + 2 * (a^2 - a')
        let (a_conj_0, a_conj_1) = (a0, -a1);
        let (A0, A1) = {
            let (tmp_0, tmp_1) = (a_square_0 - a_conj_0, a_square_1 - a_conj_1);
            (a_square_0 + tmp_0 + tmp_0, a_square_1 + tmp_1 + tmp_1)
        };

        // B = 3 * beta * c^2 + 2 * b' = beta * c^2 + 2 * (beta * c^2 + b')
        let (b_conj_0, b_conj_1) = (b0, -b1);
        (c_square_0, c_square_1) = (beta * c_square_0, beta * c_square_1);
        let (B0, B1) = {
            let (tmp_0, tmp_1) = (c_square_0 + b_conj_0, c_square_1 + b_conj_1);
            (c_square_0 + tmp_0 + tmp_0, c_square_1 + tmp_1 + tmp_1)
        };

        // C = 3 * b^2 - 2 * c' = b^2 + 2 * (b^2 - c')
        let (c_conj_0, c_conj_1) = (c0, -c1);
        let (C0, C1) = {
            let (tmp_0, tmp_1) = (b_square_0 - c_conj_0, b_square_1 - c_conj_1);
            (b_square_0 + tmp_0 + tmp_0, b_square_1 + tmp_1 + tmp_1)
        };

        Self {
            c0: Fq6 {
                c0: A0,
                c1: C0,
                c2: B0,
            },
            c1: Fq6 {
                c0: B1,
                c1: A1,
                c2: C1,
            },
        }
    }

    // trival implementation for cyclotomic group
    fn cyclotomic_exponentiation(&self, e: BigInt<MAX_SCALAR_LIMBS>) -> Self {
        let n_bits: Vec<u8> = e.to_bits();
        // make sure the highest bit is ONE
        assert_eq!(*n_bits.last().unwrap(), 1_u8);
        let (mut y, x) = (Fq12::ONE(), *self);
        for i in (0..n_bits.len()).rev() {
            y = y.cyclotomic_square();
            if n_bits[i] == 1 {
                y = y * x;
            }
        }
        y
    }
}

mod tests {
    use crate::finite_field_arithmetic::pairing_friendly::field::Field;

    use super::*;
    #[test]
    fn test_cycotomic() {}
}
