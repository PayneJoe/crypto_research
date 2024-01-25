/// Practical Implementation for Scalar Field (Fr) of Pallas Curve
///
///
use crate::finite_field_arithmetic::bigint::BigInt;
use crate::finite_field_arithmetic::traits::weierstrass_field::PrimeField;

use std::iter::Sum;
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::str::FromStr;

#[derive(Debug, PartialEq, Eq)]
pub struct FieldParseErr;

#[derive(Debug, PartialEq, Eq)]
pub struct NoneQuadraticResidualErr;

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub struct Fr<const N: usize>(pub BigInt<N>);

impl PrimeField<4> for Fr<4> {
    // fabricated precomputable parameters of custom finite field,
    // these constant parameters need to determined at compile time

    // MODULUS = 28948022309329048855892746252171976963363056481941647379679742748393362948097
    // W = 2^64
    // r = W^4 = 115792089237316195423570985008687907853269984665640564039457584007913129639936
    // R = r % MODULUS = 28948022309329048855892746252171976963180815219815621900418355762733040795645
    // R2 = r * r % MODULUS = 4263855311957679929489659445116329028194309752796460188622876710448966664207
    // R3 = r * r * r % MODULUS = 3557709630315679472311684181007729646594247341237824434526702614836137537100
    // M0 = (-MODULUS^{-1} % r) % W = 10108024940646105087
    //
    // let: MODULUS = 2^E * ROOD + 1
    // E = 32
    // RODD = 6739986666787659948666753771754907668419893943225417141728043264801
    // N = 5/7/...., is sampled number whose legendre symbol is -1 (is definitely quadratic residual)
    const MODULUS: BigInt<4> = BigInt([
        10108024940646105089,
        2469829653919213789,
        0,
        4611686018427387904,
    ]);
    const R: BigInt<4> = BigInt([
        6569413325480787965,
        11037255111951910247,
        18446744073709551615,
        4611686018427387903,
    ]);
    const R2: BigInt<4> = BigInt([
        18200867980676431887,
        7474641938123724515,
        9200329640471491984,
        679271340771891881,
    ]);
    const R3: BigInt<4> = BigInt([
        39197710403612236,
        16229805722976916262,
        9871554806900181859,
        566775843421393608,
    ]);
    const M0: u64 = 10108024940646105087 as u64;
    const E: u64 = 32 as u64;
    const RODD: BigInt<4> = BigInt([690362312389225249, 575052028, 0, 1073741824]);
    const N: BigInt<4> = BigInt([5, 0, 0, 0]);

    #[inline(always)]
    fn ONE() -> Self {
        Self(Self::R)
    }

    #[inline(always)]
    fn ZERO() -> Self {
        Self(BigInt::<4>::ZERO())
    }

    // referenced from Algorithm 11.12 of "handbook of elliptic and hyperelliptic curve cryptography"
    // (aR)^{-1} % N <- ((aR)^{-1} * R^2) * R^{-1} % N
    fn inv(&self) -> Self {
        let (mut r, mut s, mut t, mut v) = (
            self.0,
            BigInt::<4>::ONE(),
            Self::MODULUS,
            BigInt::<4>::ZERO(),
        );

        /////////////////////////////////////////////// STAGE ONE
        let mut k = 0 as usize;
        let mut tmp_carrier = 0 as u64;
        while r.is_zero() == false {
            if t.is_even() {
                // (t, s) = ((t >> 1).0, (s << 1).0);
                t = (t >> 1).0;
                (s, tmp_carrier) = s << 1;
                // FOR DEBUG
                assert!(tmp_carrier == 0);
            } else if r.is_even() {
                // (r, v) = ((r >> 1).0, (v << 1).0);
                r = (r >> 1).0;
                (v, tmp_carrier) = v << 1;
                // FOR DEBUG
                assert!(tmp_carrier == 0);
            } else if t > r {
                // (t, v, s) = (((t - r).0 >> 1).0, (v + s).0, (s << 1).0);
                t = ((t - r).0 >> 1).0;
                (v, tmp_carrier) = v + s;
                s = (s << 1).0;
                // FOR DEBUG
                assert!(tmp_carrier == 0);
            } else {
                // (r, s, v) = (((r - t).0 >> 1).0, (s + v).0, (v << 1).0);
                r = ((r - t).0 >> 1).0;
                (s, tmp_carrier) = s + v;
                v = (v << 1).0;
                // FOR DEBUG
                assert!(tmp_carrier == 0);
            }
            k = k + 1;
        }

        // v <- v mod MODULUS, make sure v is within range [0, MODULUS)
        if v >= Self::MODULUS {
            v = (v - Self::MODULUS).0;
        }
        (v, tmp_carrier) = Self::MODULUS - v;
        // FOR DEBUG
        assert!(tmp_carrier == 0);

        //////////////////////////////////////////// STAGE TWO
        let m = 4 * 64 as usize;
        // println!("------ k = {}, m = {}, v = {:?}", k, m, v);
        if k < m {
            (v, k) = (Self::mul_reduce(&v, &Self::R2).0, k + m);
        }

        // 2^{2m - k}
        let mut h_bit = 2 * m - k;
        assert!(h_bit < m);
        let x_tmp = (0..4).map(|i| i * 64 as usize).map(|v| {
            if (h_bit > v) && (h_bit - v < 64) {
                (1 as u64) << (h_bit - v)
            } else {
                0 as u64
            }
        });
        let x: Vec<u64> = (0..self.0 .0.len())
            .map(|i| {
                if h_bit >= 64 {
                    h_bit -= 64;
                    0 as u64
                } else {
                    let tmp_bit = h_bit;
                    h_bit = 0;
                    if tmp_bit > 0 {
                        1 << tmp_bit
                    } else {
                        0 as u64
                    }
                }
            })
            .collect();

        // println!("------- before reduce R2: v = {:?}", v);
        // REDC(v * R2)
        v = Self::mul_reduce(&v, &Self::R2).0;
        // println!(
        //     "------- before reduce 2^2m-k: v = {:?}, x = {:?}, k = {}",
        //     v, x, k
        // );
        // REDC(v * 2^{2m - k})
        let res = BigInt(x.try_into().unwrap());
        if res.is_zero() == false {
            v = Self::mul_reduce(&v, &res).0;
        }

        Self(v)
    }

    // Tonelli and Shanks Algorithm adapted for a special modulus p, satisfying: p MOD 8 == 1
    // referenced from Algorithm 11.23 of "handbook of elliptic and hyperelliptic curve cryptography"
    fn sqrt(&self) -> Self {
        let (n, r) = (Self::reduce(&Self::N, Some(false)), Self::RODD);
        let (a, z) = (self.clone(), n.pow(r));
        let (mut y, mut s, mut x, mut b) = (
            z,
            Self::E,
            a.pow(((r - BigInt::<4>::ONE()).0 >> 1).0),
            Self::ZERO(),
        );
        let tmp_x = a * x.clone();
        (b, x) = (tmp_x.clone() * x, tmp_x);
        while b != Self::ONE() {
            let mut m = 0 as u64;
            let mut tmp_b = b.clone();
            while tmp_b != Self::ONE() {
                tmp_b.square_inplace();
                m = m + 1;
            }
            // non-quadratic residual
            assert!(m != Self::E);
            let ss = s - m - 1;
            let mut t = y.clone();
            for i in 0..ss {
                t.square_inplace();
            }
            y = t * t;
            (s, x, b) = (m, t * x, y * b);
        }
        // FOR DEBUG
        if self.is_quadratic_residual() {
            assert!(x * x == *self);
        }
        x
    }

    fn pow(&self, e: BigInt<4>) -> Self {
        let n_bits: Vec<u8> = e.into();
        let (mut y, x) = (Self::ONE().0, self.0);
        for i in (0..n_bits.len()).rev() {
            y = Self::mul_reduce(&y, &y).0;
            if n_bits[i] == 1 {
                y = Self::mul_reduce(&x, &y).0;
            }
        }
        Self(y)
    }

    fn is_quadratic_residual(self) -> bool {
        self.pow(((Self::MODULUS - BigInt::<4>::ONE()).0 >> 1).0) == Self::ONE()
    }

    fn square(&self) -> Self {
        Self::mul_reduce(&self.0, &self.0)
    }

    fn square_inplace(&mut self) {
        *self = Self::mul_reduce(&(*self).0, &(*self).0);
    }

    // a % N <- aR * R^{-1} % N
    fn rev_reduce(&self) -> BigInt<4> {
        Self::mul_reduce(&self.0, &BigInt::<4>::ONE()).0
    }

    // referenced from Algorithm 11.3 of "handbook of elliptic and hyperelliptic curve cryptography"
    // abR % N <- (aR * bR) * R^{-1} % N
    fn mul_reduce(lft: &BigInt<4>, rht: &BigInt<4>) -> Self {
        if (*lft == BigInt::ZERO()) || (*rht == BigInt::ZERO()) {
            return Self::ZERO();
        }
        let s = 4;
        let mut t = BigInt([0 as u64; 4]);
        let (mut c1, mut c2, mut overflow) = (0 as u64, 0 as u64, false);
        // println!("lft = {:?}, rht = {:?}", lft.0, rht.0);
        for i in 0..s {
            // t = t + self * other[i]
            let (mut tmp_c1, mut tmp_c2, mut ab) = (0 as u64, 0 as u64, BigInt([0 as u64; 4]));
            (ab, tmp_c1) = lft.clone() * rht.0[i];
            (t, tmp_c2) = t + ab;
            (c1, overflow) = c1.overflowing_add(tmp_c1.wrapping_add(tmp_c2));
            if overflow {
                c2 += 1_u64;
            }
            // println!("i = {}, #1: t = {:?}, c1 = {:?}, c2 = {:?}", i, t.0, c1, c2);

            // t = t + ((t[0] * N'[0]) mod W) * N
            let (mut tmp_c3, mut tmp_c4, mut mn) = (0 as u64, 0 as u64, BigInt([0 as u64; 4]));
            let m = Self::M0.wrapping_mul(t.0[0]);
            (mn, tmp_c3) = Self::MODULUS * m;
            (t, tmp_c4) = t + mn;
            (c1, overflow) = c1.overflowing_add(tmp_c3.wrapping_add(tmp_c4));
            if overflow {
                c2 += 1_u64;
            }
            // println!("i = {} #2: t = {:?}, c1 = {:?}, c2 = {:?}", i, t.0, c1, c2);

            // t >> 1
            for j in 0..(s - 1) {
                t.0[j] = t.0[j + 1];
            }
            (t.0[s - 1], c1, c2) = (c1, c2, 0 as u64);
            // println!("i = {} #3: t = {:?}, c1 = {:?}, c2 = {:?}", i, t.0, c1, c2);
        }

        Self(t)
    }

    fn reduce(u: &BigInt<4>, inv: Option<bool>) -> Self {
        if let Some(true) = inv {
            Self::mul_reduce(u, &Self::R3)
        } else {
            Self::mul_reduce(u, &Self::R2)
        }
    }
}

impl From<&BigInt<4>> for Fr<4> {
    fn from(value: &BigInt<4>) -> Self {
        Self::reduce(value, Some(false))
    }
}

impl From<BigInt<4>> for Fr<4> {
    fn from(value: BigInt<4>) -> Self {
        Self::reduce(&value, Some(false))
    }
}

impl Into<BigInt<4>> for Fr<4> {
    fn into(self) -> BigInt<4> {
        self.rev_reduce()
    }
}

// non-reduced transformation between bytes and field
impl From<[u64; 4]> for Fr<4> {
    fn from(bytes: [u64; 4]) -> Self {
        Fr(BigInt(bytes))
    }
}

impl Into<[u64; 4]> for Fr<4> {
    fn into(self) -> [u64; 4] {
        self.0 .0
    }
}

impl<'a> Sum<&'a Self> for Fr<4> {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = &'a Self>,
    {
        iter.fold(Fr::<4>::ZERO(), |a, b| a + *b)
    }
}

impl Sum<Self> for Fr<4> {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        iter.fold(Fr::<4>::ZERO(), |a, b| a + b)
    }
}

impl Add for Fr<4> {
    type Output = Fr<4>;

    fn add(self, other: Self) -> Fr<4> {
        let (w, carrier) = self.0 + other.0;
        if carrier > 0 {
            Self(((BigInt::<4>::MAX() - Self::MODULUS).0 + w).0)
        } else if w >= Self::MODULUS {
            Self((w - Self::MODULUS).0)
        } else {
            Self(w)
        }
    }
}

impl Sub for Fr<4> {
    type Output = Fr<4>;

    fn sub(self, other: Self) -> Fr<4> {
        if self.0 > other.0 {
            Self((self.0 - other.0).0)
        } else {
            Self((Self::MODULUS - (other.0 - self.0).0).0)
        }
    }
}

impl Mul for Fr<4> {
    type Output = Fr<4>;

    fn mul(self, other: Self) -> Fr<4> {
        Self::mul_reduce(&self.0, &other.0)
    }
}

impl Mul<u64> for Fr<4> {
    type Output = Fr<4>;

    fn mul(self, other: u64) -> Fr<4> {
        Self::mul_reduce(
            &self.0,
            &Self::from_str(other.to_string().as_str()).unwrap().0,
        )
    }
}

impl Div for Fr<4> {
    type Output = Fr<4>;

    fn div(self, other: Self) -> Fr<4> {
        Self::mul_reduce(&self.0, &other.inv().0)
    }
}

impl Neg for Fr<4> {
    type Output = Fr<4>;

    fn neg(self) -> Fr<4> {
        Self::ZERO() - self
    }
}

impl FromStr for Fr<4> {
    type Err = FieldParseErr;

    fn from_str(text: &str) -> Result<Self, Self::Err> {
        let num = BigInt::<4>::from_str(text)
            .map_err(|_| FieldParseErr)
            .unwrap();
        Ok(Self::reduce(&num, Some(false)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() {
        let a = "10828745280282393011948633936436145363160692580455384354038716315440980557097";
        let a_bigint: BigInt<4> = Fr::<4>::from_str(a).unwrap().into();
        assert_eq!(a_bigint, BigInt::<4>::from_str(a).unwrap());

        let c = "591287673063831099393508774001657241885700464961659518539444896981632376260";
        let a_field: [u64; 4] = Fr::<4>::from_str(a).unwrap().into();
        assert_eq!(a_field, BigInt::<4>::from_str(c).unwrap().0);
    }

    #[test]
    fn test_addition() {
        let (a, b, c) = (
            "10828745280282393011948633936436145363160692580455384354038716315440980557097",
            "12983359841706137811367631271868112279657727965035633245325863046483061891040",
            "23812105121988530823316265208304257642818420545491017599364579361924042448137",
        );
        assert_eq!(
            Fr::<4>::from_str(a).unwrap() + Fr::<4>::from_str(b).unwrap(),
            Fr::<4>::from_str(c).unwrap()
        );
    }

    #[test]
    fn test_substraction() {
        let a = "21426167871899790663146927947666624753013764602916297221764179118297811791824";
        let b = "9973449002990857279361928142663956117381682723201036506158891067077462904966";
        let c = "11452718868908933383784999805002668635632081879715260715605288051220348886858";
        assert_eq!(
            Fr::<4>::from_str(a).unwrap() - Fr::<4>::from_str(b).unwrap(),
            Fr::<4>::from_str(c).unwrap()
        );
    }

    #[test]
    fn test_multiplication() {
        let a = "6375934151890180205297706674895575483273624160874557874408840128846376412200";
        let b = "28506874097417334274247958454240234974963135381959712862126761049363843908851";
        let c = "2254763930843862400398034612573101569234819393442628390705778063269694551656";
        assert_eq!(
            Fr::<4>::from_str(a).unwrap() * Fr::<4>::from_str(b).unwrap(),
            Fr::<4>::from_str(c).unwrap()
        );
    }

    #[test]
    fn test_inversion() {
        let a = "20478396197580483737789014158498979272870169265875691728003315572952274387533";
        let c = "21365542069856238450486712807496128750008055205983289024270465755393517385826";
        assert_eq!(
            Fr::<4>::from_str(a).unwrap().inv(),
            Fr::<4>::from_str(c).unwrap()
        );
    }

    #[test]
    fn test_division() {
        let a = "14315174832808638358939415012100039300513997680125272328694713606701545168752";
        let b = "20478396197580483737789014158498979272870169265875691728003315572952274387533";
        let c = "12111416690556434289459724424510259843801075164687209321777184926968376934270";
        assert_eq!(
            Fr::<4>::from_str(a).unwrap() / Fr::<4>::from_str(b).unwrap(),
            Fr::<4>::from_str(c).unwrap()
        );
    }

    #[test]
    fn test_sqrt() {
        let a = "23200564883514806523406480047239983027714487432999116974664880975079440510063";
        let c = "1411005539847194286406933453517763875751600735672730633893962034603782311192";
        let lft = Fr::<4>::from_str(a).unwrap();
        let result_1 = Fr::<4>::from_str(c).unwrap();

        assert_eq!(lft.is_quadratic_residual(), true);
        assert_eq!(result_1.pow(BigInt::<4>::from_str("2").unwrap()), lft);

        let result_2 = lft.sqrt();
        assert_eq!(result_2 * result_2, lft);
        println!(
            "two solutions for square root of {:?}: result_1 = {:?}, result_2 = {:?}",
            lft.0, result_1.0, result_2.0
        );
    }

    #[test]
    fn test_pow() {
        let a = "6667848649771366462575212241344760382420686144778213847978229720783476026494";
        let b = "234";
        let c = "896017193898729263273467687789996053512353124848347700387971974010744551774";
        assert_eq!(
            Fr::<4>::from_str(a)
                .unwrap()
                .pow(BigInt::<4>::from_str(b).unwrap()),
            Fr::<4>::from_str(c).unwrap()
        );
    }
}
