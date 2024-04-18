use crate::finite_field_arithmetic::bigint::BigInt;
use crate::finite_field_arithmetic::pairing_friendly::{
    bls12::fq12::Fq12, cyclotomic::CyclotomicGroup,
};

const MAX_SCALAR_LIMBS: usize = 8;

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct Gt(Fq12);

impl Gt {
    fn square(&self) -> Self {
        Self(self.0.cyclotomic_square())
    }
    fn power(&self, e: BigInt<MAX_SCALAR_LIMBS>) -> Self {
        Self(self.0.cyclotomic_exponentiation(e))
    }
}
