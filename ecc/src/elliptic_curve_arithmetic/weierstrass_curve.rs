use std::str::FromStr;

use crate::finite_field_arithmetic::field::{Field, Foo, BI};
use crate::finite_field_arithmetic::{BigInt, PrimeField};
use std::ops::{Add, Mul};

pub trait Curve {
    const a1: PrimeField;
    const a3: PrimeField;
    const a2: PrimeField;
    const a4: PrimeField;
    const a6: PrimeField;
    const IDENTITY: AffinePoint;

    // referenced from Definition 13.2 of "handbook of elliptic and hyperelliptic curve cryptography"
    fn is_nonsingular() -> bool {
        // b2 = a1^2 + 4 * a2, b4 = a1 * a3 + 2 * a4, b6 = a3^2 + 4 * a6, b8 = a1^2 * a6 - a1 * a3 * a4 + 4 * a2 * a6 + a2 * a3^2 - a4^2
        let (b2, b4, b6, b8) = (
            Self::a1 * Self::a1 + Self::a2 * (4 as u8),
            Self::a1 * Self::a3 + Self::a4 * (2 as u8),
            Self::a3 * Self::a3 + Self::a6 * (4 as u8),
            Self::a1 * Self::a1 * Self::a6 - Self::a1 * Self::a3 * Self::a4
                + Self::a2 * Self::a6 * (4 as u8)
                + Self::a2 * Self::a3 * Self::a3
                - Self::a4 * Self::a4,
        );
        // delta = -b2^2 * b8 - 8 * b4^3 - 27 * b6^2 + 9 * b2 * b4 * b6
        let delta =
            -b2 * b2 * b8 - b4 * b4 * b4 * (8 as u8) - b6 * b6 * (27 as u8) + b2 * b4 * b6 * 9;

        println!(
            "b2 = {:?}, b4 = {:?}, b6 = {:?}, b8 = {:?}, delta = {:?}",
            b2, b4, b6, b8, delta
        );
        delta != PrimeField::ZERO()
    }

    // check one point is on curve or not
    fn is_on_curve(p: &AffinePoint) -> bool;

    // evaluate y according x, it's not easy for original weierstrass curve equation
    fn to_y(x: &PrimeField) -> PrimeField;

    // negate ops on point
    fn neg(p: &AffinePoint) -> AffinePoint;

    // check whether two point negate each other or not
    fn is_negate(p1: &AffinePoint, p2: &AffinePoint) -> bool;

    // addition of two point
    fn addition(p1: &AffinePoint, p2: &AffinePoint) -> AffinePoint;

    // addition of scalar mul
    fn scalar_mul(base: &AffinePoint, scalar: &BigInt) -> AffinePoint;
}

// custom curve instance
pub struct Bar;

// implementation of custom curve instance
impl Curve for Bar {
    // y^2 + 2xy + 8y = x^3 + 5x^2 + 1136x + 531
    // for the purpose of research, we fill these constant parameters through mannual computation
    // actually these constant parameters need to be determined dynamicly at compile-time
    const a1: PrimeField = Foo(BI([217, 4]));
    const a3: PrimeField = Foo(BI([99, 6]));
    const a2: PrimeField = Foo(BI([158, 5]));
    const a4: PrimeField = Foo(BI([165, 9]));
    const a6: PrimeField = Foo(BI([43, 6]));
    const IDENTITY: AffinePoint = AffinePoint {
        x: Foo(BI([0, 0])),
        y: Foo(BI([0, 0])),
    };

    fn is_on_curve(p: &AffinePoint) -> bool {
        let (x, y) = (p.x, p.y);
        let lft = y * y + Self::a1 * x * y + Self::a3 * y;
        let rht = x * x * x + Self::a2 * x * x + Self::a4 * x + Self::a6;
        lft == rht
    }

    fn to_y(x: &PrimeField) -> PrimeField {
        unimplemented!()
    }

    fn neg(p: &AffinePoint) -> AffinePoint {
        AffinePoint {
            x: p.x,
            y: -p.y - Self::a1 * p.x - Self::a3,
        }
    }

    fn is_negate(p1: &AffinePoint, p2: &AffinePoint) -> bool {
        (p1.x == p2.x) && (p1.y + p2.y == -Self::a1 * p1.x - Self::a3)
    }

    fn addition(p1: &AffinePoint, p2: &AffinePoint) -> AffinePoint {
        if *p1 == Self::IDENTITY {
            return p2.clone();
        }
        if *p2 == Self::IDENTITY {
            return p1.clone();
        }
        if Self::is_negate(p1, p2) {
            return Self::IDENTITY;
        }
        let (x1, y1, x2, y2) = (p1.x, p1.y, p2.x, p2.y);
        // chord and tangle method
        let (denominator, nominator) = if p1 == p2 {
            (
                x1 * x1 * (3 as u8) + Self::a2 * x1 * (2 as u8) + Self::a4 - Self::a1 * y1,
                y1 * (2 as u8) + Self::a1 * x1 + Self::a3,
            )
        } else {
            (y1 - y2, x1 - x2)
        };
        let lambda = denominator * nominator.inv();
        let x3 = lambda * lambda + Self::a1 * lambda - Self::a2 - x1 - x2;
        AffinePoint {
            x: x3,
            y: lambda * (x1 - x3) - y1 - Self::a1 * x3 - Self::a3,
        }
    }

    fn scalar_mul(base: &AffinePoint, scalar: &BigInt) -> AffinePoint {
        unimplemented!()
    }
}

#[derive(PartialEq, Eq, Clone, Debug)]
pub struct AffinePoint {
    x: PrimeField,
    y: PrimeField,
}

impl AffinePoint {
    fn from_bigint(x: BigInt, y: BigInt) -> Self {
        Self {
            x: PrimeField::from(x),
            y: PrimeField::from(y),
        }
    }
    fn from_str(x: &str, y: &str) -> Self {
        Self {
            x: PrimeField::from_str(x).unwrap(),
            y: PrimeField::from_str(y).unwrap(),
        }
    }
    fn is_on_curve(&self) -> bool {
        Bar::is_on_curve(&self)
    }
}

impl Add for AffinePoint {
    type Output = AffinePoint;

    fn add(self, other: AffinePoint) -> Self::Output {
        Bar::addition(&self, &other)
    }
}

impl Mul<BigInt> for AffinePoint {
    type Output = AffinePoint;

    fn mul(self, other: BigInt) -> Self::Output {
        Bar::scalar_mul(&self, &other)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_nonsingular() {
        assert_eq!(Bar::is_nonsingular(), true)
    }

    #[test]
    fn test_is_on_curve() {
        let a = ("", "");
        let lft = AffinePoint::from_str(a.0.to_string().as_str(), a.1.to_string().as_str());
        assert_eq!(lft.is_on_curve(), true);
    }

    #[test]
    fn test_addition() {
        let (a, b, c) = (("", ""), ("", ""), ("", ""));
        let lft = AffinePoint::from_str(a.0.to_string().as_str(), a.1.to_string().as_str());
        let rht = AffinePoint::from_str(b.0.to_string().as_str(), b.1.to_string().as_str());
        let result = AffinePoint::from_str(c.0.to_string().as_str(), c.1.to_string().as_str());
        assert_eq!(lft + rht, result);
    }

    #[test]
    fn test_scalar_mul() {
        let (a, b, c) = (("", ""), "", ("", ""));
        let lft = AffinePoint::from_str(a.0.to_string().as_str(), a.1.to_string().as_str());
        let rht = BigInt::from_str(b.to_string().as_str()).unwrap();
        let result = AffinePoint::from_str(c.0.to_string().as_str(), c.1.to_string().as_str());
        assert_eq!(lft * rht, result);
    }
}
