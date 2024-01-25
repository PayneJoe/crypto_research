use std::str::FromStr;

use crate::finite_field_arithmetic::weierstrass_field::{Field, Foo, BI};
use crate::finite_field_arithmetic::{WSBigInt, WSField};
use std::ops::{Add, Mul, Neg};

pub trait Curve {
    // parameters of standard weierstrass curve
    const a1: WSField;
    const a3: WSField;
    const a2: WSField;
    const a4: WSField;
    const a6: WSField;
    // identity point on curve
    const IDENTITY: AffinePoint;
    // generator
    const GENERATOR: AffinePoint;
    // order
    const ORDER: WSBigInt;

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
        delta != WSField::ZERO()
    }

    // check one point is on curve or not
    fn is_on_curve(p: &AffinePoint) -> bool;

    // evaluate y according x, it's not easy for original weierstrass curve equation
    fn to_y(x: &WSField) -> WSField;

    // negate ops on point
    fn neg(p: &AffinePoint) -> AffinePoint;

    // check whether two point negate each other or not
    fn is_negate(p1: &AffinePoint, p2: &AffinePoint) -> bool;

    // addition of two point
    fn addition(p1: &AffinePoint, p2: &AffinePoint) -> AffinePoint;

    // addition of scalar mul
    fn scalar_mul(base: &AffinePoint, scalar: &WSBigInt) -> AffinePoint;
}

// custom curve instance
pub struct Bar;

// implementation of custom curve instance
impl Curve for Bar {
    // y^2 + 2xy + 8y = x^3 + 5x^2 + 1136x + 531
    // for the purpose of research, we fill these constant parameters through mannual computation
    // actually these constant parameters need to be determined dynamicly at compile-time
    const a1: WSField = Foo(BI([217, 4]));
    const a3: WSField = Foo(BI([99, 6]));
    const a2: WSField = Foo(BI([158, 5]));
    const a4: WSField = Foo(BI([165, 9]));
    const a6: WSField = Foo(BI([43, 6]));
    // faked identity representing the point at infinity
    const IDENTITY: AffinePoint = AffinePoint {
        x: Foo(BI([0, 0])),
        y: Foo(BI([0, 0])),
    };
    // generator of this elliptic curve, g = (2535, 2835)
    const GENERATOR: AffinePoint = AffinePoint {
        x: Foo(BI([15, 0])),
        y: Foo(BI([254, 11])),
    };
    // order of this elliptic curve is 3295, which is not a prime number
    const ORDER: WSBigInt = BI([223, 12]);

    fn is_on_curve(p: &AffinePoint) -> bool {
        let (x, y) = (p.x, p.y);
        let lft = y * y + Self::a1 * x * y + Self::a3 * y;
        let rht = x * x * x + Self::a2 * x * x + Self::a4 * x + Self::a6;
        lft == rht
    }

    fn to_y(x: &WSField) -> WSField {
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

    // referenced from P.270 of "handbook of elliptic and hyperelliptic curve cryptography"
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
        // let tmp_nom: BI<2> = nominator.into();
        // let tmp_inv: BI<2> = nominator.inv().into();
        // println!("******* inv({:?}) = {:?}", tmp_nom, tmp_inv);
        let lambda = denominator * nominator.inv();
        let x3 = lambda * lambda + Self::a1 * lambda - Self::a2 - x1 - x2;
        AffinePoint {
            x: x3,
            y: lambda * (x1 - x3) - y1 - Self::a1 * x3 - Self::a3,
        }
    }

    // referenced from Algorithm 13.6 of "handbook of elliptic and hyperelliptic curve cryptography"
    fn scalar_mul(base: &AffinePoint, scalar: &WSBigInt) -> AffinePoint {
        // k < 8, make sure u8 is big enough for store precomputated points
        let k = 4;
        let n = scalar.to_bits();

        // precomputation table
        let mut table = vec![base.clone()];
        let double_base = base + base;
        for i in 1..(1 << (k - 1)) {
            table.push(&table[i - 1] + &double_base)
        }
        // for i in 0..table.len() {
        //     let tmp_x: BI<2> = table[i].x.into();
        //     let tmp_y: BI<2> = table[i].y.into();
        //     print!("[{}]p = ({:?}, {:?}) \t", 2 * i + 1, tmp_x, tmp_y);
        // }
        // println!("n = {:?} \n\n", n);

        let (mut q, mut i) = (Self::IDENTITY, n.len() - 1);
        while i != 0 {
            // left shift skipping zeros, doubling
            if n[i] == 0 {
                (q, i) = (&q + &q, i - 1);
                // let tmp_x_3: BI<2> = q.x.into();
                // let tmp_y_3: BI<2> = q.y.into();
                // println!(
                //     "**** Doubling: i = {}, after doubliing q = ({:?}, {:?})",
                //     i, tmp_x_3, tmp_y_3
                // );
            } else {
                // left shift, doubling
                let mut s = std::cmp::max(i as i32 - k as i32 + 1, 0) as usize;
                let right_bound = s;
                while n[s] == 0 {
                    s = s + 1;
                }
                let left_bound = s;
                for _ in 0..(i - s + 1) {
                    q = &q + &q;
                }
                // let tmp_x_1: BI<2> = q.x.into();
                // let tmp_y_1: BI<2> = q.y.into();
                // println!(
                //     "**** Doubling: i = {}, after doubling q = ({:?}, {:?}), {}-{}",
                //     i, tmp_x_1, tmp_y_1, left_bound, right_bound
                // );

                // then addition with precomputated table
                let u = bits_to_u8(&n[s..(i + 1)]);
                q = &q + &table[((u - 1) / 2) as usize];
                // let tmp_x_2: BI<2> = q.x.into();
                // let tmp_y_2: BI<2> = q.y.into();
                // println!(
                //     "**** Addition: i = {}, after addition q = ({:?}, {:?}), u = {}",
                //     i, tmp_x_2, tmp_y_2, u
                // );
                i = if s >= 1 { s - 1 } else { 0 }
            }
        }

        q
    }
}

// little-endian
fn bits_to_u8(bits: &[u8]) -> u8 {
    assert!(bits.len() < 8);
    let mut result = 0;
    for i in 0..bits.len() {
        if bits[i] == 1_u8 {
            result = result + (1 << i);
        }
    }
    result
}

#[derive(PartialEq, Eq, Clone, Debug)]
pub struct AffinePoint {
    x: WSField,
    y: WSField,
}

impl AffinePoint {
    fn from_bigint(x: WSBigInt, y: WSBigInt) -> Self {
        Self {
            x: WSField::from(x),
            y: WSField::from(y),
        }
    }
    fn from_str(x: &str, y: &str) -> Self {
        Self {
            x: WSField::from_str(x).unwrap(),
            y: WSField::from_str(y).unwrap(),
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

impl<'a, 'b> Add<&'b AffinePoint> for &'a AffinePoint {
    type Output = AffinePoint;

    fn add(self, other: &'b AffinePoint) -> Self::Output {
        Bar::addition(self, other)
    }
}

impl Mul<WSBigInt> for AffinePoint {
    type Output = AffinePoint;

    fn mul(self, other: WSBigInt) -> Self::Output {
        Bar::scalar_mul(&self, &other)
    }
}

impl Neg for AffinePoint {
    type Output = AffinePoint;

    fn neg(self) -> Self::Output {
        Bar::neg(&self)
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
        // 2 * g
        let a = ("2251", "2594");
        let lft = AffinePoint::from_str(a.0.to_string().as_str(), a.1.to_string().as_str());
        assert_eq!(lft.is_on_curve(), true);
    }

    #[test]
    fn test_addition() {
        // 2 * g + 3 * g == 5 * g
        // let (a, b, c) = (("2251", "2594"), ("2820", "805"), ("70", "1903"));

        // 10 * g + 4 * g == 14 * g
        let (a, b, c) = (("2158", "2165"), ("149", "130"), ("1889", "431"));
        let lft = AffinePoint::from_str(a.0.to_string().as_str(), a.1.to_string().as_str());
        let rht = AffinePoint::from_str(b.0.to_string().as_str(), b.1.to_string().as_str());
        let result = AffinePoint::from_str(c.0.to_string().as_str(), c.1.to_string().as_str());
        println!("lft = {:?}, rht = {:?}, result = {:?}", lft, rht, result);
        assert_eq!(lft + rht, result);
    }

    #[test]
    fn test_is_negate() {
        let (a, c) = (("2251", "2594"), ("2251", "2883"));
        let lft = AffinePoint::from_str(a.0.to_string().as_str(), a.1.to_string().as_str());
        let result = AffinePoint::from_str(c.0.to_string().as_str(), c.1.to_string().as_str());
        assert_eq!(-lft, result);
    }

    #[test]
    fn test_scalar_mul_simple() {
        // (2 * g) * 3 == 6 * g
        let (a, b, c) = (("2251", "2594"), "3", ("1323", "1896"));
        let lft = AffinePoint::from_str(a.0.to_string().as_str(), a.1.to_string().as_str());
        let rht = WSBigInt::from_str(b.to_string().as_str()).unwrap();
        let result = AffinePoint::from_str(c.0.to_string().as_str(), c.1.to_string().as_str());
        assert_eq!(lft * rht, result);
    }

    #[test]
    fn test_scalar_mul_bigger() {
        // (2 * g) * 325 == 650 * g
        let (a, b, c) = (("2251", "2594"), "325", ("1103", "2591"));
        let lft = AffinePoint::from_str(a.0.to_string().as_str(), a.1.to_string().as_str());
        let rht = WSBigInt::from_str(b.to_string().as_str()).unwrap();
        let result = AffinePoint::from_str(c.0.to_string().as_str(), c.1.to_string().as_str());
        assert_eq!(lft * rht, result);
    }
}
