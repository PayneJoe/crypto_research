use crate::finite_field_arithmetic::bigint::BigInt;
use crate::finite_field_arithmetic::traits::weierstrass_field::{
    Field, LegendreSymbol, PrimeField,
};
use std::iter::Chain;
use std::ops::{Add, Div, Mul, Neg, Sub};

pub trait QuadraticExtensionConfig<const N: usize>: Copy + Clone {
    type BasePrimeField: PrimeField<N>;
    type BaseField: Field<N, BasePrimeField = Self::BasePrimeField>;
    type FrobCoeff: Field<N>;

    // absolute degree of current (extension) field
    const DEGREE_OVER_BASE_PRIME_FIELD: usize;

    // alpha in Fq[X] / X^2 - alpha
    const NON_QUADRATIC_RESIDUAL: Self::BaseField;

    // frobenius coefficients
    const FROBENIUS_COEFF_C1: [Self::FrobCoeff];
}

// there are two coefficients in quadratic extension field
#[derive(Copy, Clone)]
pub struct QuadraticExtension<const N: usize, Config: QuadraticExtensionConfig<N>> {
    pub c0: Config::BaseField,
    pub c1: Config::BaseField,
}

impl<const N: usize, Config: QuadraticExtensionConfig<N>> QuadraticExtension<N, Config> {
    pub fn new(c0: Config::BaseField, c1: Config::BaseField) -> Self {
        Self { c0, c1 }
    }

    // norm = (c0 - c1 * u) * (c0 + c1 * u) = c0^2 - c1^2 * u^2 = c0^2 - c1^2 * alpha
    pub fn norm(&self) -> Config::BaseField {
        self.c0.square() - self.c1.square() * Config::NON_QUADRATIC_RESIDUAL
    }

    // frobenius map over quadratic extension field is just its conjugative one
    // more details about this come from my another note: https://hackmd.io/@70xfCGp1QViTYYJh3AMrQg/rJZ-A_M1R
    // (a_0 + a_1 u)^q
    // = a_0^q + a_1^q * u^q
    // = a_0^q + a_1^q * u^{q - 1} * u
    // = a_0^q + a_1^q * alpha^{(q - 1)/2} * u
    // = a_0^q + a_1^q * (-1) * u
    // = a_0^q - a_1^q u
    pub fn frobenius_map(&self) -> Self {
        Self::new(self.c0, -self.c1)
    }
}

type BaseFieldIter<const N: usize, Config> =
    <<Config as QuadraticExtensionConfig<N>>::BaseField as Field<N>>::BasePrimeFieldIter;

impl<const N: usize, Config: QuadraticExtensionConfig<N>> Field<N>
    for QuadraticExtension<N, Config>
{
    type BasePrimeField = Config::BasePrimeField;
    type BasePrimeFieldIter = Chain<BaseFieldIter<N, Config>, BaseFieldIter<N, Config>>;

    fn extension_degree() -> u64 {
        2 * Config::BaseField::extension_degree()
    }

    fn to_base_prime_field_elements(&self) -> Self::BasePrimeFieldIter {
        self.c0
            .to_base_prime_field_elements()
            .chain(self.c1.to_base_prime_field_elements())
    }

    fn from_base_prime_field_elems(
        elems: impl IntoIterator<Item = Self::BasePrimeField>,
    ) -> Option<Self> {
        let mut elems = elems.into_iter();
        let elems = elems.by_ref();
        let base_ext_deg = Config::BaseField::extension_degree() as usize;
        let element = Some(Self::new(
            Config::BaseField::from_base_prime_field_elems(elems.take(base_ext_deg)).unwrap(),
            Config::BaseField::from_base_prime_field_elems(elems.take(base_ext_deg)).unwrap(),
        ));
        // degree is not proper
        if elems.next().is_some() {
            None
        } else {
            element
        }
    }

    // whether current field element has square root or not depends its Norm has square root or not
    // resolve this legendre symbol recursively
    // x^{(q^2 - 1)/2}
    // = (x^{q + 1})^{(q - 1)/2}
    // = (\Phi(x) * x)^{(q - 1)/2}
    // = (x' * x)^{(q - 1)/2}
    // = Norm(x)^{(q - 1)/2}
    fn legendre(&self) -> LegendreSymbol {
        self.norm().legendre()
    }

    // apply frobenius map power times
    fn powers_frobenius_map(&self, power: usize) -> Self {
        let mut i = power;
        let mut result = self.clone();
        while i > 0 {
            result = result.frobenius_map();
            i -= 1;
        }
        result
    }

    fn square(&self) -> Self {
        todo!()
    }

    fn sqrt(&self) -> Option<Self> {
        todo!()
    }

    fn inverse(&self) -> Option<Self> {
        todo!()
    }

    fn pow(&self, e: BigInt<N>) -> Self {
        todo!()
    }
}

impl<const N: usize, Config: QuadraticExtensionConfig<N>> Neg for QuadraticExtension<N, Config> {
    type Output = QuadraticExtension<N, Config>;

    fn neg(self) -> QuadraticExtension<N, Config> {
        Self::new(-self.c0, -self.c1)
    }
}

impl<const N: usize, Config: QuadraticExtensionConfig<N>> Add for QuadraticExtension<N, Config> {
    type Output = QuadraticExtension<N, Config>;

    fn add(self, other: Self) -> QuadraticExtension<N, Config> {
        Self::new(self.c0 + other.c0, self.c1 + other.c1)
    }
}

impl<const N: usize, Config: QuadraticExtensionConfig<N>> Sub for QuadraticExtension<N, Config> {
    type Output = QuadraticExtension<N, Config>;

    fn sub(self, other: Self) -> QuadraticExtension<N, Config> {
        Self::new(self.c0 - other.c0, self.c1 - other.c1)
    }
}

// impl<'a, const N: usize, Config: QuadraticExtensionConfig<N>> Mul<&'a QuadraticExtension<N, Config>>
//     for QuadraticExtension<N, Config>
// {
//     type Output = QuadraticExtension<N, Config>;
//
//     fn mul(self, other: &Self) -> QuadraticExtension<N, Config> {
//         todo!()
//     }
// }

// referenced from Algorithm 5.16 (P133) of "Guide to Pairing-Based Cryptography"
impl<const N: usize, Config: QuadraticExtensionConfig<N>> Mul for QuadraticExtension<N, Config> {
    type Output = QuadraticExtension<N, Config>;

    fn mul(self, other: Self) -> QuadraticExtension<N, Config> {
        let v0 = self.c0 * other.c0;
        let v1 = self.c1 * other.c1;
        let c0 = v0 + Config::NON_QUADRATIC_RESIDUAL * v1;
        let c1 = (self.c0 + self.c1) * (other.c0 + other.c1) - v0 - v1;
        Self::new(c0, c1)
    }
}

impl<const N: usize, Config: QuadraticExtensionConfig<N>> Div for QuadraticExtension<N, Config> {
    type Output = QuadraticExtension<N, Config>;

    fn div(self, other: Self) -> QuadraticExtension<N, Config> {
        todo!()
    }
}
