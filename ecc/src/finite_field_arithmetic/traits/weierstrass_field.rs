use crate::finite_field_arithmetic::bigint64::BigInt;
use std::marker::PhantomData;

use std::{
    ops::{Add, Div, Mul, Shl, Shr, Sub},
    str::FromStr,
};

// pub trait PrimeFieldConfig<const N: usize> {
//     const MODULUS: BigInt<N>;
//     const R: BigInt<N>;
//     const R2: BigInt<N>;
//     const R3: BigInt<N>;
//     const M0: u64;
//     const E: u64;
//     const RODD: u64;
//     const N: BigInt<N>;

//     fn inv(f: &PrimeField<Self, N>) -> PrimeField<Self, N>;
//     fn square(f: &PrimeField<Self, N>) -> PrimeField<Self, N>;
//     fn square_in_place(f: &mut PrimeField<Self, N>);
//     fn pow(f: &PrimeField<Self, N>, exp: &BigInt<N>) -> PrimeField<Self, N>;

//     fn add(lft: &PrimeField<Self, N>, rht: &PrimeField<Self, N>) -> PrimeField<Self, N>;
//     fn mul(lft: &PrimeField<Self, N>, rht: &PrimeField<Self, N>) -> PrimeField<Self, N>;
// }

// #[derive(Clone, Debug)]
// pub struct PrimeField<P: PrimeFieldConfig<N>, const N: usize>(pub BigInt<N>, pub PhantomData<P>);

// impl<const N: usize> PrimeFieldConfig<N> for PrimeField<Self, N> {

// }

// impl<P: PrimeFieldConfig, const N: usize> Add for PrimeField<P, N> {
//     type Output = PrimeField<P, N>;
//     fn add(self, other: PrimeField<P, N>) -> Self::Output {
//         Self::add(&self, &other)
//     }
// }

// impl<P: PrimeFieldConfig, const N: usize> Mul for PrimeField<P, N> {
//     type Output = PrimeField<P, N>;
//     fn mul(self, other: PrimeField<P, N>) -> Self::Output {
//         Self::mul(&self, &other)
//     }
// }

pub trait PrimeField<const N: usize>:
    FromStr
    + From<BigInt<N>>
    + Into<BigInt<N>>
    + From<[u64; N]>
    + Into<[u64; N]>
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
    const M0: u64;
    // for the convenient of square root (Tonelli and Shanks Algorithm)
    const E: u64;
    const RODD: u64;
    const N: BigInt<N>;

    fn ONE() -> Self;
    fn ZERO() -> Self;

    // legendre symbol
    // fn legendre_symbol(self) -> bool;

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
