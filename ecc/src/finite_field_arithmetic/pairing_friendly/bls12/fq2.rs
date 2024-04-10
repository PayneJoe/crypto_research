use crate::finite_field_arithmetic::bigint::BigInt;
use crate::finite_field_arithmetic::pairing_friendly::bls12::fq::Fq;
use crate::finite_field_arithmetic::pairing_friendly::field::{Field, PrimeField};
use crate::finite_field_arithmetic::pairing_friendly::quadratic_extension::{
    QuadraticExtension, QuadraticExtensionConfig,
};

const NUM_LIMBS: usize = 6;

#[derive(Copy, Clone)]
pub struct Fq2Config;

impl QuadraticExtensionConfig<NUM_LIMBS> for Fq2Config {
    type BasePrimeField = Fq<NUM_LIMBS>;
    type BaseField = Fq<NUM_LIMBS>;
    type FrobCoeff = Fq<NUM_LIMBS>;

    // absolute degree of current (extension) field
    const DEGREE_OVER_BASE_PRIME_FIELD: usize = 2;

    // Fq[X] / X^2 - alpha, where alpha = -1 in Fq2/Fq over bls12
    const NON_QUADRATIC_RESIDUAL: Fq<NUM_LIMBS> = Fq(BigInt([
        4897101644811774638,
        3654671041462534141,
        569769440802610537,
        17053147383018470266,
        17227549637287919721,
        291242102765847046,
    ]));

    // coefficients of frobenius map over Fp2
    // alpha^{(p^d - 1)/2}
    // [1, -1] for bls12-381
    const FROBENIUS_COEFF_C1: &'static [Fq<NUM_LIMBS>] = &[
        Fq(BigInt([
            8505329371266088957,
            17002214543764226050,
            6865905132761471162,
            8632934651105793861,
            6631298214892334189,
            1582556514881692819,
        ])),
        Fq(BigInt([
            4897101644811774638,
            3654671041462534141,
            569769440802610537,
            17053147383018470266,
            17227549637287919721,
            291242102765847046,
        ])),
    ];

    fn multiply_frobenius_coeff(c: &mut Self::BaseField, power: usize) {
        *c = *c * Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
    }
}

// fq2 is an abstract type of quadratic extension over base prime field Fq
// F_q^2 = F_q[X] / X^2 - alpha
pub type Fq2 = QuadraticExtension<NUM_LIMBS, Fq2Config>;
impl Fq2 {}

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
