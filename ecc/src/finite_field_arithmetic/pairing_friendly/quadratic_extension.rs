use crate::finite_field_arithmetic::bigint::BigInt;
use crate::finite_field_arithmetic::pairing_friendly::field::{Field, LegendreSymbol, PrimeField};
use std::iter::Chain;
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::str::FromStr;

pub trait QuadraticExtensionConfig<const N: usize>: Copy + Clone + Sized + 'static {
    type BasePrimeField: PrimeField<N>;
    type BaseField: Field<N, BasePrimeField = Self::BasePrimeField>;
    type FrobCoeff: Field<N>;

    // absolute degree of current (extension) field
    const DEGREE_OVER_BASE_PRIME_FIELD: usize;

    // alpha in Fq[X] / X^2 - alpha
    const NON_QUADRATIC_RESIDUAL: Self::BaseField;

    // precomputed constant coefficients of frobenius map on quadratic extension
    // namely alpha^{(p^d - 1)/2}
    const FROBENIUS_COEFF_C1: &'static [Self::FrobCoeff];

    fn multiply_frobenius_coeff(c: &mut Self::BaseField, power: usize);
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

    // construct current field with base prime field elements
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

    // construct current field with single base prime field element
    fn from_base_prime_field_elem(elem: Self::BasePrimeField) -> Self {
        let base_field_degree = Config::BaseField::extension_degree() as usize;

        let mut c0 = vec![elem];
        c0.extend(vec![Config::BasePrimeField::ZERO(); base_field_degree - 1].iter());
        let c1 = vec![Config::BasePrimeField::ZERO(); base_field_degree];

        Self::new(
            Config::BaseField::from_base_prime_field_elems(c0.into_iter()).unwrap(),
            Config::BaseField::from_base_prime_field_elems(c1.into_iter()).unwrap(),
        )
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

    // frobenius map over quadratic extension field
    // recursively apply frobenius map power times
    //
    // (a_0 + a_1 u)^q
    // = a_0^q + a_1^q * u^q
    // = a_0^q + a_1^q * u^{q - 1} * u
    // = a_0^q + a_1^q * alpha^{(q - 1)/2} * u
    // more details about this come from my another note: https://hackmd.io/@70xfCGp1QViTYYJh3AMrQg/rJZ-A_M1R
    fn powers_frobenius_map_inplace(&mut self, power: usize) {
        self.c0.powers_frobenius_map_inplace(power);
        self.c1.powers_frobenius_map_inplace(power);
        Config::multiply_frobenius_coeff(&mut self.c1, power);
    }

    // referenced from Algorithm 5.17 (P133) of "Guide to Pairing-Based Cryptography"
    fn square(&self) -> Self {
        let mut v0 = self.c0 - self.c1;
        let v3 = self.c0 - Config::NON_QUADRATIC_RESIDUAL * self.c1;
        let v2 = self.c0 * self.c1;
        v0 = (v0 * v3) + v2;
        let c1 = v2 + v2;
        let c0 = v0 + Config::NON_QUADRATIC_RESIDUAL * v2;
        Self::new(c0, c1)
    }

    fn square_inplace(&mut self) {
        *self = self.square();
    }

    // referenced from Algorithm 5.18 (P133) of "Guide to Pairing-Based Cryptography"
    fn sqrt(&self) -> Option<Self> {
        // simple case, when c1 = 0
        if self.c1.is_zero() {
            if self.c0.legendre() == LegendreSymbol::QuadraticResidue {
                return Some(Self::new(
                    self.c0.sqrt().unwrap(),
                    Config::BaseField::ZERO(),
                ));
            } else {
                return None;
            }
        }

        // not quadratic residual
        if self.legendre() == LegendreSymbol::QuadraticNonResidue {
            return None;
        }

        // quadratic residual definitely
        // norm(a) is absolutely quadratic residual
        let mut lambda = self.norm();
        lambda = lambda.sqrt().unwrap();
        let constant_2 = Config::BasePrimeField::from(BigInt::<N>::from_str("2").unwrap());
        let base_constant_2 = Config::BaseField::from_base_prime_field_elem(constant_2);
        let mut delta = (self.c0 + lambda) / base_constant_2;
        if delta.legendre() == LegendreSymbol::QuadraticNonResidue {
            delta = (self.c0 - lambda) / base_constant_2;
        }
        let (c0, c1) = (delta.sqrt().unwrap(), self.c1 / (base_constant_2 * self.c0));

        Some(Self::new(c0, c1))
    }

    // make sure that self or its norm is not zero
    fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            let nm_inv = self.norm().inverse();
            if nm_inv.is_some() {
                Some(Self::new(
                    nm_inv.unwrap() * self.c0,
                    -nm_inv.unwrap() * self.c1,
                ))
            } else {
                None
            }
        }
    }

    fn pow(&self, e: BigInt<N>) -> Self {
        let n_bits: Vec<u8> = e.to_bits();
        let (mut y, x) = (Self::ONE(), *self);
        for i in (0..n_bits.len()).rev() {
            y.square_inplace();
            if n_bits[i] == 1 {
                y = y * x;
            }
        }
        y
    }

    // is zero or not
    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    fn ZERO() -> Self {
        Self::new(Config::BaseField::ZERO(), Config::BaseField::ZERO())
    }

    fn ONE() -> Self {
        Self::new(Config::BaseField::ONE(), Config::BaseField::ZERO())
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
