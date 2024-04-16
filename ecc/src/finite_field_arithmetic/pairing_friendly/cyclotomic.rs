use crate::finite_field_arithmetic::bigint::BigInt;
use crate::finite_field_arithmetic::pairing_friendly::field::Field;

const MAX_SCALAR_LIMBS: usize = 8 as usize;

//// cyclotomic multiplicative subgroup of Fpk
pub trait CyclotomicGroup<const N: usize>: Field<N> {
    fn cyclotomic_inverse(&self) -> Option<Self>;
    fn cyclotomic_square(&self) -> Self;
    fn cyclotomic_exponentiation(&self, power: BigInt<MAX_SCALAR_LIMBS>) -> Self;
}
