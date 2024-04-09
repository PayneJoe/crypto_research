use crate::finite_field_arithmetic::pairing_friendly::field::{Field, PrimeField};

pub trait CubicExtensionConfig<const N: usize>: Copy + Clone + Sized + 'static {
    type BasePrimeField: PrimeField<N>;
    type BaseField: Field<N, BasePrimeField = Self::BasePrimeField>;
    type FrobCoeff: Field<N>;

    // absolute degree of current (extension) field
    const DEGREE_OVER_BASE_PRIME_FIELD: usize;

    // beta in Fq[X] / X^3 - beta
    const NON_CUBIC_RESIDUAL: Self::BaseField;

    // frobenius coefficients
    const FROBENIUS_COEFF_C1: [Self::FrobCoeff; 3];
}

#[derive(Copy, Clone)]
pub struct CubicExtension<const N: usize, Config: CubicExtensionConfig<N>> {
    pub c0: Config::BaseField,
    pub c1: Config::BaseField,
    pub c2: Config::BaseField,
}
