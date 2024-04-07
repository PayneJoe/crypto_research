use crate::finite_field_arithmetic::bigint::BigInt;

use std::fmt::Debug;
use std::{
    ops::{Add, Div, Mul, Neg, Shl, Shr, Sub},
    str::FromStr,
};

#[derive(Debug, PartialEq, Eq)]
pub enum LegendreSymbol {
    Zero = 0,
    QuadraticResidue = 1,
    QuadraticNonResidue = -1,
}

// trait of base field F_q
pub trait PrimeField<const N: usize>:
    FromStr
    + Debug
    + PartialEq
    + Eq
    + From<BigInt<N>>
    + Into<BigInt<N>>
    + From<[u64; N]>
    + Into<[u64; N]>
    + Clone
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Neg<Output = Self>
    + Mul<Self, Output = Self>
    + for<'a> std::iter::Sum<&'a Self>
    + std::iter::Sum<Self>
    + Copy
{
    // finite field modulus
    const MODULUS: BigInt<N>;
    // W^s % MODULUS, for Montgomery reduce
    const R: BigInt<N>;
    // W^2s % MODULUS, for Montgomery reduce
    const R2: BigInt<N>;
    // W^3s % MODULUS, for Montgomery reduce
    const R3: BigInt<N>;
    // inversion of least significant word of modulus, also convenient for Montgomery reduce
    const M0: u64;
    // for the convenient of square root (Tonelli and Shanks Algorithm)
    const E: u64;
    const RODD: BigInt<N>;
    const N: BigInt<N>;

    fn ONE() -> Self;
    fn ZERO() -> Self;

    // legendre symbol
    // fn legendre_symbol(self) -> bool;

    fn is_quadratic_residual(self) -> bool;
    fn is_zero(self) -> bool;

    // reduce a bigint into the specified range [0, modulus)
    fn reduce(u: &BigInt<N>, inv: Option<bool>) -> Self;
    // convert a reduced number into a unreduced bigint
    fn rev_reduce(&self) -> BigInt<N>;

    //////////////////////////// basic operations on field
    /// F^{1/2}
    fn sqrt(&self) -> Self;
    // F^-1
    fn inv(&self) -> Self;
    // F^2
    fn square_inplace(&mut self);
    fn square(&self) -> Self;
    // F^e
    fn pow(&self, e: BigInt<N>) -> Self;
    // F * F
    fn mul_reduce(lft: &BigInt<N>, rht: &BigInt<N>) -> Self;

    //
    fn to_bytes(self) -> Vec<u8>;
    fn random() -> Self;

    fn to_string(self) -> String;
}

// trait of extension field F_q^k over F_q
pub trait Field<const N: usize>:
    Sized
    + Copy
    + Clone
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    // + for<'a> Mul<&'a Self, Output = Self>
    + Mul<Self, Output = Self>
    + Div<Self, Output = Self>
    + Neg<Output = Self>
{
    // base prime field F_q
    type BasePrimeField: PrimeField<N>;

    // coefficients with length of k
    type BasePrimeFieldIter: Iterator<Item = Self::BasePrimeField>;

    // specially for k, absolute degree of extension field F_q^k
    fn extension_degree() -> u64;

    // output coefficients of extension field
    fn to_base_prime_field_elements(&self) -> Self::BasePrimeFieldIter;

    // construct extension field from base prime field
    fn from_base_prime_field_elems(
        elems: impl IntoIterator<Item = Self::BasePrimeField>,
    ) -> Option<Self>;

    // whether is capable of squaring within current (extension) field F_q^k or not
    fn legendre(&self) -> LegendreSymbol;

    // frobenius map, a <- a^q
    fn powers_frobenius_map(&self, power: usize) -> Self;

    // arithemtics on extension field
    fn square(&self) -> Self;
    fn sqrt(&self) -> Option<Self>;
    fn inverse(&self) -> Option<Self>;
    fn pow(&self, e: BigInt<N>) -> Self;
}
