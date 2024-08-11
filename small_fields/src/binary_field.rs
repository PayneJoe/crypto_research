use std::ops::{Add, BitXor, Mul, Neg, Sub};

#[allow(non_snake_case)]
pub trait SmallField:
    PartialEq
    + Eq
    + Clone
    + Copy
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Neg<Output = Self>
    + Mul<Self, Output = Self>
{
    // fn ONE() -> Self;
    // fn ZERO() -> Self;
    // fn pow() -> Self;
    // fn inv() -> Self;
}

pub trait BinaryTowerConfig {
    type BaseBinaryField: SmallField;
}

pub struct BinaryTowerField<Config: BinaryTowerConfig> {
    pub c0: Config::BaseBinaryField,
    pub c1: Config::BaseBinaryField,
}

impl<Config: BinaryTowerConfig> BinaryTowerField<Config> {
    fn new(c0: Config::BaseBinaryField, c1: Config::BaseBinaryField) -> Self {
        Self { c0, c1 }
    }
}

impl<Config: BinaryTowerConfig> Add for BinaryTowerField<Config> {
    type Output = BinaryTowerField<Config>;
    fn add(self, other: Self) -> BinaryTowerField<Config> {
        Self {
            c0: self.c0 + other.c0,
            c1: self.c1 + other.c1,
        }
    }
}

impl<Config: BinaryTowerConfig> Neg for BinaryTowerField<Config> {
    type Output = BinaryTowerField<Config>;
    fn neg(self) -> BinaryTowerField<Config> {
        self
    }
}

impl<Config: BinaryTowerConfig> Sub for BinaryTowerField<Config> {
    type Output = BinaryTowerField<Config>;
    fn sub(self, other: Self) -> BinaryTowerField<Config> {
        Self {
            c0: self.c0 + other.c0,
            c1: self.c1 + other.c1,
        }
    }
}

impl<Config: BinaryTowerConfig> Mul for BinaryTowerField<Config> {
    type Output = BinaryTowerField<Config>;
    fn mul(self, other: Self) -> BinaryTowerField<Config> {
        let a0b0 = self.c0 * other.c0;
        let a1b1 = self.c1 * other.c1;

        todo!()
    }
}
