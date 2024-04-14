use crate::finite_field_arithmetic::bigint::BigInt;
use crate::finite_field_arithmetic::pairing_friendly::field::Field;

use std::marker::PhantomData;
use std::ops::{Add, Mul, Neg, Sub};

use std::str::FromStr;

const WINDOW_SIZE: usize = 6;
const BASE_NUM_LIMBS: usize = 6;
const SCALAR_NUM_LIMBS: usize = 4;
const WORD_SIZE: usize = 64;
type Word = u64;

// type BigInteger = BigInt<NUM_LIMBS>;

#[derive(PartialEq, Eq, Clone, Debug, Copy, Hash)]
pub struct AffinePoint<F: Field<BASE_NUM_LIMBS>, R: Field<SCALAR_NUM_LIMBS>, C: Curve<F, R>> {
    pub x: F,
    pub y: F,
    pub _p1: PhantomData<R>,
    pub _p2: PhantomData<C>,
}

impl<F: Field<BASE_NUM_LIMBS>, R: Field<SCALAR_NUM_LIMBS>, C: Curve<F, R>> AffinePoint<F, R, C> {
    pub fn new(x: F, y: F) -> Self {
        Self {
            x,
            y,
            _p1: Default::default(),
            _p2: Default::default(),
        }
    }

    pub fn IDENTITY() -> Self {
        C::IDENTITY
    }

    pub fn from_bigint(x: F, y: F) -> Self {
        Self::new(F::from(x), F::from(y))
    }

    pub fn is_on_curve(&self) -> bool {
        C::is_on_curve(&self)
    }
}

impl<F: Field<BASE_NUM_LIMBS>, R: Field<SCALAR_NUM_LIMBS>, C: Curve<F, R>> Add
    for AffinePoint<F, R, C>
{
    type Output = AffinePoint<F, R, C>;

    fn add(self, other: AffinePoint<F, R, C>) -> Self::Output {
        C::addition(&self, &other)
    }
}

impl<'a, 'b, F: Field<BASE_NUM_LIMBS>, R: Field<SCALAR_NUM_LIMBS>, C: Curve<F, R>>
    Add<&'b AffinePoint<F, R, C>> for &'a AffinePoint<F, R, C>
{
    type Output = AffinePoint<F, R, C>;

    fn add(self, other: &'b AffinePoint<F, R, C>) -> Self::Output {
        C::addition(self, other)
    }
}

impl<'a, 'b, F: Field<BASE_NUM_LIMBS>, R: Field<SCALAR_NUM_LIMBS>, C: Curve<F, R>> Sub
    for AffinePoint<F, R, C>
{
    type Output = AffinePoint<F, R, C>;

    fn sub(self, other: AffinePoint<F, R, C>) -> Self::Output {
        C::addition(&self, &-other)
    }
}

impl<'a, 'b, F: Field<BASE_NUM_LIMBS>, R: Field<SCALAR_NUM_LIMBS>, C: Curve<F, R>>
    Sub<&'b AffinePoint<F, R, C>> for &'a AffinePoint<F, R, C>
{
    type Output = AffinePoint<F, R, C>;

    fn sub(self, other: &'b AffinePoint<F, R, C>) -> Self::Output {
        C::addition(self, &-other)
    }
}

impl<F: Field<BASE_NUM_LIMBS>, R: Field<SCALAR_NUM_LIMBS>, C: Curve<F, R>> Mul<R>
    for AffinePoint<F, R, C>
{
    type Output = AffinePoint<F, R, C>;

    fn mul(self, other: R) -> Self::Output {
        C::scalar_mul(&self, &other)
    }
}

impl<F: Field<BASE_NUM_LIMBS>, R: Field<SCALAR_NUM_LIMBS>, C: Curve<F, R>> Neg
    for AffinePoint<F, R, C>
{
    type Output = AffinePoint<F, R, C>;

    fn neg(self) -> Self::Output {
        C::neg(&self)
    }
}

impl<F: Field<BASE_NUM_LIMBS>, R: Field<SCALAR_NUM_LIMBS>, C: Curve<F, R>> Neg
    for &AffinePoint<F, R, C>
{
    type Output = AffinePoint<F, R, C>;

    fn neg(self) -> Self::Output {
        C::neg(self)
    }
}

pub trait Curve<BaseField: Field<BASE_NUM_LIMBS>, ScalarField: Field<SCALAR_NUM_LIMBS>>:
    Sized + Copy + Clone
{
    // parameters of standard weierstrass curve
    const a4: BaseField;
    const a6: BaseField;
    // identity point on curve
    const IDENTITY: AffinePoint<BaseField, ScalarField, Self>;
    // generator
    const GENERATOR: AffinePoint<BaseField, ScalarField, Self>;
    // order
    const ORDER: BigInt<SCALAR_NUM_LIMBS>;

    // referenced from Definition 13.2 of "handbook of elliptic and hyperelliptic curve cryptography"
    fn is_nonsingular() -> bool {
        // b2 = a1^2 + 4 * a2, b4 = a1 * a3 + 2 * a4, b6 = a3^2 + 4 * a6, b8 = a1^2 * a6 - a1 * a3 * a4 + 4 * a2 * a6 + a2 * a3^2 - a4^2
        let (b2, b4, b6, b8) = (
            BaseField::ZERO(),
            Self::a4 * BaseField::from_small_base_prime_field_str("2"),
            Self::a6 * BaseField::from_small_base_prime_field_str("4"),
            -Self::a4 * Self::a4,
        );
        // delta = -b2^2 * b8 - 8 * b4^3 - 27 * b6^2 + 9 * b2 * b4 * b6
        let delta = -b2 * b2 * b8
            - b4 * b4 * b4 * BaseField::from_small_base_prime_field_str("8")
            - b6 * b6 * BaseField::from_small_base_prime_field_str("27")
            + b2 * b4 * b6 * BaseField::from_small_base_prime_field_str("9");

        println!(
            "b2 = {:?}, b4 = {:?}, b6 = {:?}, b8 = {:?}, delta = {:?}",
            b2, b4, b6, b8, delta
        );
        delta != BaseField::ZERO()
    }

    // check one point is on curve or not
    fn is_on_curve(p: &AffinePoint<BaseField, ScalarField, Self>) -> bool;

    // evaluate y according x, it's not easy for original weierstrass curve equation
    fn to_y(x: &BaseField) -> BaseField;

    // negate ops on point
    fn neg(
        p: &AffinePoint<BaseField, ScalarField, Self>,
    ) -> AffinePoint<BaseField, ScalarField, Self>;

    // check whether two point negate each other or not
    fn is_negate(
        p1: &AffinePoint<BaseField, ScalarField, Self>,
        p2: &AffinePoint<BaseField, ScalarField, Self>,
    ) -> bool;

    // addition of two point
    fn addition(
        p1: &AffinePoint<BaseField, ScalarField, Self>,
        p2: &AffinePoint<BaseField, ScalarField, Self>,
    ) -> AffinePoint<BaseField, ScalarField, Self>;

    // addition of scalar mul
    fn scalar_mul(
        base: &AffinePoint<BaseField, ScalarField, Self>,
        scalar: &ScalarField,
    ) -> AffinePoint<BaseField, ScalarField, Self>;
}
