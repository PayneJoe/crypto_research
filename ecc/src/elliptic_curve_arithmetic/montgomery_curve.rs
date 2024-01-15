use std::ops::Mul;
use std::str::FromStr;

use crate::finite_field_arithmetic::field_mont_friendly::{Field, Foo, BI};

type BigInt = BI<2>;
type PrimeField = Foo<2>;

pub trait Curve {
    // Montgomery form: By^2 = x^3 + Ax^2 + x
    // Weierstrass form: y^2 = x^3 + a_4x + a_6
    // Projective Coordinate of Weierstrass form: Y^2 Z = X^3 + a_4 X Z^2 + a_6 Z^3
    //
    // Conversion from Montgomery to Weierstrass form:
    // a_4 = 1 / B^2 − A^2 / 3B^2
    // a6 = −A^3 / 27B^3 − a^4 (A / 3B)
    const B: PrimeField;
    const A: PrimeField;
    // identity point on curve
    const IDENTITY: ProjectivePoint;
    // generator
    const GENERATOR: ProjectivePoint;
    // order
    const ORDER: BigInt;

    // check one point is on curve or not
    fn is_on_curve(p: &AffinePoint) -> bool;
    fn addition(
        p1: &ProjectivePoint,
        p2: &ProjectivePoint,
        base: &ProjectivePoint,
    ) -> ProjectivePoint;
    fn doubling(p: &ProjectivePoint) -> ProjectivePoint;
    fn scalar_mul(base: &ProjectivePoint, scalar: &BigInt) -> ProjectivePoint;
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AffinePoint {
    x: PrimeField,
    y: PrimeField,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProjectivePoint {
    x: PrimeField,
    y: PrimeField,
    z: PrimeField,
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

impl ProjectivePoint {
    fn from_affine_point(p: &AffinePoint) -> Self {
        Self {
            x: p.x,
            y: p.y,
            z: PrimeField::ONE(),
        }
    }
    fn from_bigint(x: BigInt, y: BigInt, z: BigInt) -> Self {
        Self {
            x: PrimeField::from(x),
            y: PrimeField::from(y),
            z: PrimeField::from(z),
        }
    }
    fn from_str(x: &str, y: &str, z: &str) -> Self {
        Self {
            x: PrimeField::from_str(x).unwrap(),
            y: PrimeField::from_str(y).unwrap(),
            z: PrimeField::from_str(z).unwrap(),
        }
    }
}

impl Mul<BigInt> for ProjectivePoint {
    type Output = ProjectivePoint;

    fn mul(self, other: BigInt) -> Self::Output {
        Bar::scalar_mul(&self, &other)
    }
}

#[derive(Debug, Clone)]
pub struct Bar;

impl Curve for Bar {
    // Weierstrass Curve: y^2 = x^3 + 1132x + 278
    // Montgomery Curve: 899 y^2 = x^3 + 1421 x^2 + x
    // identity = (0, 1, 0), generator = (271 : 1864 : 1), order = 1956
    //
    // Since MODULUS = 2003, where exists alpha = 1702 satisfying:
    // 1. alpha^3 + 1132 alpha + 278 = 0
    // 2. m = 3 alpha^2 + a_4 = 527, is quadratic residual number relate to MODULUS 2003
    //
    // Therefore there is a isomorphism between Weierstrass curve and Montgomery curve:
    // (x, y) = (899(x - 1702), 899y)
    // for instance, point p = (1120, 1391) on Weierstrass curve, will be mapped to point p_e = (1568, 637) on Montgomery curve
    const B: PrimeField = Foo(BI([110, 2]));
    const A: PrimeField = Foo(BI([153, 4]));
    const IDENTITY: ProjectivePoint = ProjectivePoint {
        x: Foo(BI([0, 0])),
        y: Foo(BI([160, 5])),
        z: Foo(BI([0, 0])),
    };
    const GENERATOR: ProjectivePoint = ProjectivePoint {
        x: Foo(BI([122, 6])),
        y: Foo(BI([140, 0])),
        z: Foo(BI([160, 5])),
    };
    const ORDER: BigInt = BI([164, 7]);

    fn is_on_curve(p: &AffinePoint) -> bool {
        unimplemented!()
    }

    // X_m+n = Z_m-n * ((X_m - Z_m) * (X_n + Z_n) + (X_m + Z_m) * (X_n - Z_n))^2
    // Z_m+n = X_m-n * ((X_m - Z_m) * (X_n + Z_n) - (X_m + Z_m) * (X_n - Z_n))^2
    fn addition(
        p1: &ProjectivePoint,
        p2: &ProjectivePoint,
        base: &ProjectivePoint,
    ) -> ProjectivePoint {
        // 4M + 2S
        let root_lft = (p1.x - p1.z) * (p2.x + p2.z);
        let root_rht = (p1.x + p1.z) * (p2.x - p2.z);
        let addition_square = (root_lft + root_rht).square();
        let minus_square = (root_lft - root_rht).square();
        ProjectivePoint {
            x: base.z * addition_square,
            y: PrimeField::ZERO(),
            z: base.x * minus_square,
        }
    }

    // 4X_n Z_n = (X_n + Z_n)^2 - (X_n - Z_n)^2
    // X_2n = (X_n + Z_n)^2 * (X_n - Z_n)^2
    // Z_2n = 4X_n Z_n * ((X_n - Z_n)^2 + ((A + 2) / 4) * (4 X_n Z_n))
    fn doubling(p: &ProjectivePoint) -> ProjectivePoint {
        // 3M + 2S
        let (xn_add_zn_square, xn_sub_zn_square) = ((p.x + p.z).square(), (p.x - p.z).square());
        let xn_zn_4times = xn_add_zn_square - xn_sub_zn_square;
        let constant_z = (Self::A + PrimeField::from(BI([2, 0]))) / PrimeField::from(BI([4, 0]));
        ProjectivePoint {
            x: xn_add_zn_square * xn_sub_zn_square,
            y: PrimeField::ZERO(),
            z: xn_zn_4times * (xn_sub_zn_square + constant_z * xn_zn_4times),
        }
    }

    // Montgomery Ladder, referenced from Algorithm 13.35 of "handbook of elliptic and hyperelliptic curve cryptography"
    fn scalar_mul(base: &ProjectivePoint, scalar: &BigInt) -> ProjectivePoint {
        let n = scalar.to_bits();
        if n.len() == 1 {
            let result = if n[0] == 0 {
                base.clone()
            } else {
                Self::IDENTITY
            };
            return result;
        }
        let (mut p1, mut p2) = (base.clone(), Self::doubling(base));
        for i in (0..(n.len() - 1)).rev() {
            (p1, p2) = if n[i] == 0_u8 {
                (Self::doubling(&p1), Self::addition(&p1, &p2, base))
            } else {
                (Self::addition(&p1, &p2, base), Self::doubling(&p2))
            };
        }
        return p1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scalar_mul_simple() {
        let (a, b, c) = (("1568", "637", "1"), "2", ("35", "0", "1887"));
        let lft = ProjectivePoint::from_str(
            a.0.to_string().as_str(),
            a.1.to_string().as_str(),
            a.2.to_string().as_str(),
        );
        let rht = BigInt::from_str(b.to_string().as_str()).unwrap();
        let result = ProjectivePoint::from_str(
            c.0.to_string().as_str(),
            c.1.to_string().as_str(),
            c.2.to_string().as_str(),
        );
        assert_eq!(lft * rht, result);
    }

    #[test]
    fn test_scalar_mul_bigger() {
        let (a, b, c) = (("1568", "637", "1"), "763", ("568", "0", "746"));
        let lft = ProjectivePoint::from_str(
            a.0.to_string().as_str(),
            a.1.to_string().as_str(),
            a.2.to_string().as_str(),
        );
        let rht = BigInt::from_str(b.to_string().as_str()).unwrap();
        let result = ProjectivePoint::from_str(
            c.0.to_string().as_str(),
            c.1.to_string().as_str(),
            c.2.to_string().as_str(),
        );
        assert_eq!(lft * rht, result);
    }
}
