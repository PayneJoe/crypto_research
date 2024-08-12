use std::ops::{Add, Mul, Neg, Sub};

use crate::bits::BitVectorize;

pub trait BaseBinaryField: BinaryTowerField + Copy + Clone + Eq + PartialEq {}

#[allow(non_snake_case)]
pub trait BinaryTowerField:
    PartialEq
    + Eq
    + Clone
    + Copy
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Neg<Output = Self>
    + Mul<Self, Output = Self>
{
    type BaseBinaryField: BaseBinaryField;
    fn ONE() -> Self;
    fn ZERO() -> Self;
    fn GEN() -> Self;

    fn extension_degree() -> u32;
    fn square(self) -> Self;
    fn double(self) -> Self;
    fn pow(self, e: u128) -> Self;
    fn inv(self) -> Self;
}

pub trait BinaryTowerConfig: Copy + Clone + Sized + Eq + PartialEq {
    type BaseBinaryField: BaseBinaryField;
    type BaseField: BinaryTowerField;
}

#[derive(Clone, Copy, Eq, PartialEq, Debug)]
pub struct BinaryTower<Config: BinaryTowerConfig> {
    pub c0: Config::BaseField,
    pub c1: Config::BaseField,
}

impl<Config: BinaryTowerConfig> BinaryTower<Config> {
    pub fn new(c0: Config::BaseField, c1: Config::BaseField) -> Self {
        Self { c0, c1 }
    }
}

impl<Config: BinaryTowerConfig> BinaryTowerField for BinaryTower<Config> {
    type BaseBinaryField = Config::BaseBinaryField;
    fn ONE() -> Self {
        Self {
            c0: Config::BaseField::ONE(),
            c1: Config::BaseField::ZERO(),
        }
    }

    fn ZERO() -> Self {
        Self {
            c0: Config::BaseField::ZERO(),
            c1: Config::BaseField::ZERO(),
        }
    }

    fn GEN() -> Self {
        Self {
            c0: Config::BaseField::ZERO(),
            c1: Config::BaseField::ONE(),
        }
    }

    fn extension_degree() -> u32 {
        2 * Config::BaseField::extension_degree()
    }

    fn double(self) -> Self {
        Self {
            c0: Config::BaseField::ZERO(),
            c1: Config::BaseField::ZERO(),
        }
    }

    fn pow(self, e: u128) -> Self {
        if e == 0 {
            return Self::ONE();
        }
        let bits = e.to_le_bits();
        let mut result = self.clone();
        for i in 1..bits.len() {
            result = result.square();
            if bits[i] == 1 {
                result = result * self;
            }
        }
        result
    }

    /// Little-Fermat Lemma: 1 / a = a^{p - 2}
    fn inv(self) -> Self {
        let p = 1 << Self::extension_degree();
        assert!(p < (1 << 128));
        self.pow(p - 2)
    }

    fn square(self) -> Self {
        let t0 = self.c0.square();
        let t1 = self.c1.square();
        let t2 = t1 * Config::BaseField::GEN();
        let t3 = (self.c0 + self.c1).square();
        Self {
            c0: t0 + t1,
            c1: t2 + t3 + t0 + t1,
        }
    }
}

impl<Config: BinaryTowerConfig> Add for BinaryTower<Config> {
    type Output = BinaryTower<Config>;
    fn add(self, other: Self) -> BinaryTower<Config> {
        Self {
            c0: self.c0 + other.c0,
            c1: self.c1 + other.c1,
        }
    }
}

impl<Config: BinaryTowerConfig> Neg for BinaryTower<Config> {
    type Output = BinaryTower<Config>;
    fn neg(self) -> BinaryTower<Config> {
        self
    }
}

impl<Config: BinaryTowerConfig> Sub for BinaryTower<Config> {
    type Output = BinaryTower<Config>;
    fn sub(self, other: Self) -> BinaryTower<Config> {
        Self {
            c0: self.c0 + other.c0,
            c1: self.c1 + other.c1,
        }
    }
}

// a * b = (a_0 + a_1 X_n) * (b_0 + b_1 X_n)
//       = (a_0 * b_0 + a_1 * b_1) + (a_1 * b_1 * X_{n - 1} + a_0 * b_1 + a_1 * b_0) X_n
//       = (a_0 * b_0 + a_1 * b_1) + (a_1 * b_1 * X_{n - 1} + (a_0 + a_1)*(b_0 + b_1) + a_0 * b_0 + a_1 * b_1) X_n
impl<Config: BinaryTowerConfig> Mul for BinaryTower<Config> {
    type Output = BinaryTower<Config>;
    fn mul(self, other: Self) -> BinaryTower<Config> {
        let t0 = self.c0 * other.c0;
        let t1 = self.c1 * other.c1;
        let t2 = t1 * Config::BaseField::GEN();
        let t3 = (self.c0 + self.c1) * (other.c0 + other.c1);
        Self {
            c0: t0 + t1,
            c1: t2 + t3 + t0 + t1,
        }
    }
}
