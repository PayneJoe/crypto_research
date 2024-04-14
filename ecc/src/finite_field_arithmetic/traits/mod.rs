pub mod weierstrass_field;

use crate::finite_field_arithmetic::bigint::BigInt;
use std::{
    ops::{Add, Div, Mul, Shl, Shr, Sub},
    str::FromStr,
};

pub trait Field<const N: usize>:
    FromStr
    + From<BigInt<N>>
    + Into<BigInt<N>>
    + From<[u8; N]>
    + Into<[u8; N]>
    + Clone
    + Add<Self, Output = Self>
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
    const M0: u8;
    // for the convenient of square root (Tonelli and Shanks Algorithm)
    // const E: u8;
    // const RODD: u8;
    // const N: BI<N>;

    fn ONE() -> Self;
    fn ZERO() -> Self;

    // legendre symbol
    fn legendre_symbol(self) -> bool;

    fn is_quadratic_residual(self) -> bool;

    // reduce a bigint into the specified range [0, modulus)
    fn reduce(u: &BigInt<N>, inv: Option<bool>) -> Self;
    // convert a reduced number into a unreduced bigint
    fn rev_reduce(&self) -> BigInt<N>;

    //////////////////////////// basic operations on field
    // F^-1
    fn inv(&self) -> Self;
    // F^2
    fn square_inplace(&mut self);
    fn square(&self) -> Self;
    // F^e
    fn pow(&self, e: BigInt<N>) -> Self;
    // F * F
    fn mul_reduce(lft: &BigInt<N>, rht: &BigInt<N>) -> Self;
}
