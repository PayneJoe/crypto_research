use crate::finite_field_arithmetic::bigint::BigInt;
use crate::finite_field_arithmetic::pairing_friendly::field::{Field, LegendreSymbol, PrimeField};
use std::fmt::Debug;
use std::iter::Chain;
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::str::FromStr;
pub trait CubicExtensionConfig<const N: usize>:
    Copy + Clone + Sized + 'static + Debug + Eq + PartialEq
{
    type BasePrimeField: PrimeField<N>;
    type BaseField: Field<N, BasePrimeField = Self::BasePrimeField>;
    type FrobCoeff: Field<N>;

    // absolute degree of current (extension) field
    const DEGREE_OVER_BASE_PRIME_FIELD: usize;

    // beta in Fq[X] / X^3 - beta
    const NON_CUBIC_RESIDUAL: Self::BaseField;

    // precomputed constant coefficients of frobenius map on cubic extension
    // namely beta^{(p^d - 1)/3}
    // two frobenius coefficient list, one for c1, another for c2, each with length of DEGREE_OVER_BASE_PRIME_FIELD
    const FROBENIUS_COEFF_C1: &'static [Self::FrobCoeff];
    const FROBENIUS_COEFF_C2: &'static [Self::FrobCoeff];

    fn multiply_frobenius_coeff(c1: &mut Self::BaseField, c2: &mut Self::BaseField, power: usize);
}

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct CubicExtension<const N: usize, Config: CubicExtensionConfig<N>> {
    pub c0: Config::BaseField,
    pub c1: Config::BaseField,
    pub c2: Config::BaseField,
}

impl<const N: usize, Config: CubicExtensionConfig<N>> CubicExtension<N, Config> {
    pub fn new(c0: Config::BaseField, c1: Config::BaseField, c2: Config::BaseField) -> Self {
        Self { c0, c1, c2 }
    }

    // referenced from DEFINITION 2.27 (P50) of "Guide to Pairing-based Cryptography"
    // |a| = \prod_i^n a^{p^i}
    // |a| = a * a^p * a^{p^2}
    // where p = q^n
    pub fn norm(&self) -> Config::BaseField {
        let n = Config::BaseField::extension_degree() as usize;
        let mut multiplier = *self;
        let mut result = *self;
        for i in 1..3 {
            multiplier.powers_frobenius_map_inplace(i * n);
            result = result * multiplier;
        }
        assert!(result.c1.is_zero() && result.c2.is_zero());
        result.c0
    }
}

type BaseFieldIter<const N: usize, Config> =
    <<Config as CubicExtensionConfig<N>>::BaseField as Field<N>>::BasePrimeFieldIter;

impl<const N: usize, Config: CubicExtensionConfig<N>> Field<N> for CubicExtension<N, Config> {
    type BasePrimeField = Config::BasePrimeField;
    type BasePrimeFieldIter =
        Chain<BaseFieldIter<N, Config>, Chain<BaseFieldIter<N, Config>, BaseFieldIter<N, Config>>>;

    fn extension_degree() -> u64 {
        3 * Config::BaseField::extension_degree()
    }

    fn to_base_prime_field_elements(&self) -> Self::BasePrimeFieldIter {
        self.c0.to_base_prime_field_elements().chain(
            self.c1
                .to_base_prime_field_elements()
                .chain(self.c2.to_base_prime_field_elements()),
        )
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
        let c2 = vec![Config::BasePrimeField::ZERO(); base_field_degree];

        Self::new(
            Config::BaseField::from_base_prime_field_elems(c0.into_iter()).unwrap(),
            Config::BaseField::from_base_prime_field_elems(c1.into_iter()).unwrap(),
            Config::BaseField::from_base_prime_field_elems(c2.into_iter()).unwrap(),
        )
    }

    // whether current field element has square root or not depends its Norm has square root or not
    // resolve this legendre symbol recursively
    //
    // x^{(q^3 - 1)/2}
    // = (x^{q^2 + q + 1})^{(q - 1)/2}
    // = (\Phi(x)^2 * \Phi(X) * x)^{(q - 1)/2}
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
        self.c2.powers_frobenius_map_inplace(power);
        Config::multiply_frobenius_coeff(&mut self.c1, &mut self.c2, power);
    }

    // referenced from Algorithm 5.22 (P136) of "Guide to Pairing-Based Cryptography"
    fn square(&self) -> Self {
        let constant_2 = Config::BasePrimeField::from(BigInt::<N>::from_str("2").unwrap());
        let base_constant_2 = Config::BaseField::from_base_prime_field_elem(constant_2);

        let mut v4 = base_constant_2 * (self.c0 * self.c1);
        let mut v5 = self.c2.square();
        let c1 = Config::NON_CUBIC_RESIDUAL * v5 + v4;
        let v2 = v4 - v5;
        let v3 = self.c0.square();
        v4 = self.c0 - self.c1 + self.c2;
        v5 = base_constant_2 * (self.c1 * self.c2);
        v4 = v4.square();
        let c0 = Config::NON_CUBIC_RESIDUAL * v5 + v3;
        let c2 = v2 + v4 + v5 - v3;
        Self::new(c0, c1, c2)
    }

    fn square_inplace(&mut self) {
        *self = self.square();
    }

    // referenced from Algorithm 5.23 (P137) of "Guide to Pairing-Based Cryptography"
    fn inverse(&self) -> Option<Self> {
        let (v0, v1, v2) = (self.c0.square(), self.c1.square(), self.c2.square());
        let (v3, v4, v5) = (self.c0 * self.c1, self.c0 * self.c2, self.c1 * self.c2);
        let A = v0 - Config::NON_CUBIC_RESIDUAL * v5;
        let B = Config::NON_CUBIC_RESIDUAL * v2 - v3;
        let C = v1 - v4;
        let mut v6 = self.c0 * A;
        v6 = v6 + (Config::NON_CUBIC_RESIDUAL * self.c2 * B);
        v6 = v6 + (Config::NON_CUBIC_RESIDUAL * self.c1 * C);
        let F = v6.inverse().unwrap();
        let (c0, c1, c2) = (A * F, B * F, C * F);
        Some(Self::new(c0, c1, c2))
    }

    // is zero or not
    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    fn is_one(&self) -> bool {
        self.c0.is_one() && self.c1.is_zero() && self.c2.is_zero()
    }

    fn ZERO() -> Self {
        Self::new(
            Config::BaseField::ZERO(),
            Config::BaseField::ZERO(),
            Config::BaseField::ZERO(),
        )
    }

    fn ONE() -> Self {
        Self::new(
            Config::BaseField::ONE(),
            Config::BaseField::ZERO(),
            Config::BaseField::ZERO(),
        )
    }

    fn sqrt(&self) -> Option<Self> {
        todo!()
    }

    // naive implementation of exponentiation
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
}

impl<const N: usize, Config: CubicExtensionConfig<N>> Neg for CubicExtension<N, Config> {
    type Output = CubicExtension<N, Config>;

    fn neg(self) -> CubicExtension<N, Config> {
        Self::new(-self.c0, -self.c1, -self.c2)
    }
}

impl<const N: usize, Config: CubicExtensionConfig<N>> Add for CubicExtension<N, Config> {
    type Output = CubicExtension<N, Config>;

    fn add(self, other: Self) -> CubicExtension<N, Config> {
        if self.is_zero() {
            return other;
        }
        if other.is_zero() {
            return self;
        }
        Self::new(self.c0 + other.c0, self.c1 + other.c1, self.c2 + other.c2)
    }
}

impl<const N: usize, Config: CubicExtensionConfig<N>> Sub for CubicExtension<N, Config> {
    type Output = CubicExtension<N, Config>;

    fn sub(self, other: Self) -> CubicExtension<N, Config> {
        if self.is_zero() {
            return -other;
        }
        if other.is_zero() {
            return self;
        }
        Self::new(self.c0 - other.c0, self.c1 - other.c1, self.c2 - other.c2)
    }
}

// impl<'a, const N: usize, Config: CubicExtensionConfig<N>> Mul<&'a CubicExtension<N, Config>>
//     for CubicExtension<N, Config>
// {
//     type Output = CubicExtension<N, Config>;
//
//     fn mul(self, other: &Self) -> CubicExtension<N, Config> {
//         todo!()
//     }
// }

// referenced from Algorithm 5.21 (P135) of "Guide to Pairing-Based Cryptography"
impl<const N: usize, Config: CubicExtensionConfig<N>> Mul for CubicExtension<N, Config> {
    type Output = CubicExtension<N, Config>;

    fn mul(self, other: Self) -> CubicExtension<N, Config> {
        if self.is_zero() {
            return Self::ZERO();
        }
        if other.is_zero() {
            return Self::ZERO();
        }
        if self.is_one() {
            return other;
        }
        if other.is_one() {
            return self;
        }
        let v0 = self.c0 * other.c0;
        let v1 = self.c1 * other.c1;
        let v2 = self.c2 * other.c2;
        let c0 = ((self.c1 + self.c2) * (other.c1 + other.c2) - v1 - v2)
            * Config::NON_CUBIC_RESIDUAL
            + v0;
        let c1 =
            (self.c0 + self.c1) * (other.c0 + other.c1) - v0 - v1 + Config::NON_CUBIC_RESIDUAL * v2;
        let c2 = (self.c0 + self.c2) * (other.c0 + other.c2) - v0 - v2 + v1;
        Self::new(c0, c1, c2)
    }
}

impl<const N: usize, Config: CubicExtensionConfig<N>> Div for CubicExtension<N, Config> {
    type Output = CubicExtension<N, Config>;

    fn div(self, other: Self) -> CubicExtension<N, Config> {
        if other.is_one() {
            return self;
        }
        self * other.inverse().unwrap()
    }
}
