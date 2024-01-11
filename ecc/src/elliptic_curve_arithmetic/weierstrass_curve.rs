use std::str::FromStr;

use crate::finite_field_arithmetic::field::{Exponentiation, Field, Foo, BI};
use crate::finite_field_arithmetic::{BigInt, PrimeField};

pub trait Curve {
    const a1: PrimeField;
    const a3: PrimeField;
    const a2: PrimeField;
    const a4: PrimeField;
    const a6: PrimeField;

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

        delta != PrimeField::ZERO()
    }

    // check one point is on curve or not
    fn is_on_curve(p: &AffinePoint) -> bool;

    // evaluate y according x, it's not easy for original weierstrass curve equation
    fn to_y(x: &PrimeField) -> PrimeField;
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

    fn is_on_curve(p: &AffinePoint) -> bool {
        let (x, y) = (p.x, p.y);
        let lft = y * y + Self::a1 * x * y + Self::a3 * y;
        let rht = x * x * x + Self::a2 * x * x + Self::a4 * x + Self::a6;
        lft == rht
    }

    fn to_y(x: &PrimeField) -> PrimeField {
        unimplemented!()
    }
}

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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_nonsingular() {
        assert_eq!(Bar::is_nonsingular(), true)
    }
}
