use super::binary_tower::{BaseBinaryField, BinaryTower, BinaryTowerConfig, BinaryTowerField};
use std::ops::{Add, Mul, Neg, Sub};

/// Fp
#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct Fp(pub u32);

impl Fp {
    pub fn new(v: u32) -> Self {
        assert!((v == 0) || (v == 1));
        Self(v)
    }
}

impl Add for Fp {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self::new(self.0 ^ other.0)
    }
}

impl Neg for Fp {
    type Output = Fp;
    fn neg(self) -> Self {
        Self::new(1 - self.0)
    }
}

impl Sub for Fp {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        self + (-other)
    }
}

impl Mul for Fp {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Self::new(self.0 & other.0)
    }
}

impl BaseBinaryField for Fp {}

impl BinaryTowerField for Fp {
    type BaseBinaryField = Self;
    fn ONE() -> Self {
        Self::new(1)
    }
    fn ZERO() -> Self {
        Self::new(0)
    }
    fn GEN() -> Self {
        Self::new(1)
    }
    fn extension_degree() -> u32 {
        1
    }
    fn double(self) -> Self {
        Self::new(0)
    }
    fn square(self) -> Self {
        self
    }
    #[allow(unused_variables)]
    fn pow(self, e: u128) -> Self {
        if e == 0 {
            return Self::new(1);
        }
        self
    }
    fn inv(self) -> Self {
        Self::new(1)
    }
}

/// Fp2
#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct Fp2Config;

impl BinaryTowerConfig for Fp2Config {
    type BaseBinaryField = Fp;
    type BaseField = Fp;
}

pub type Fp2 = BinaryTower<Fp2Config>;

impl Fp2 {}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_fp2() {
        let a = Fp2::new(Fp(0), Fp(1));
        let b = Fp2::new(Fp(1), Fp(1));

        assert_eq!(a + b, Fp2::new(Fp(1), Fp(0)));
        assert_eq!(a * b, Fp2::new(Fp(1), Fp(0)));
        assert_eq!(a.pow(17_u128), Fp2::new(Fp(1), Fp(1)));
        assert_eq!(a.pow(12_u128), Fp2::new(Fp(1), Fp(0)));
        assert_eq!(a.inv(), Fp2::new(Fp(1), Fp(1)));
        assert_eq!(a * a.inv(), Fp2::ONE());
    }
}
