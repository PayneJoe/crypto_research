use crate::finite_field_arithmetic::bigint64::BigInt;
use crate::finite_field_arithmetic::traits::weierstrass_field::PrimeField;

use std::iter::Sum;
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::str::FromStr;

#[derive(Debug)]
pub struct FieldParseErr;

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub struct Fr<const N: usize>(pub BigInt<N>);

impl PrimeField<4> for Fr<4> {
    // fabricated precomputable parameters of custom finite field,
    // these constant parameters need to determined at compile time

    // W = 256, MODULUS/M = 2003, R = 256^2 % MODULUS = 982, M0 = (-M^{-1} % 256^2) % W =
    // MODULUS % 8 \noteq 1
    const MODULUS: BigInt<4> = BigInt([211, 7, 0, 0]);
    const R: BigInt<4> = BigInt([160, 5, 0, 0]);
    const R2: BigInt<4> = BigInt([239, 1, 0, 0]);
    const R3: BigInt<4> = BigInt([199, 6, 0, 0]);
    const M0: u64 = 165 as u64;
    // 2^E * RODD = MODULUS - 1
    const E: u64 = 8 as u64;
    const RODD: u64 = 13 as u64;
    // N is a sampled non-quadratic residual number, 3/6/11/12/...
    const N: BigInt<4> = BigInt([3, 0, 0, 0]);

    // const MODULUS: BigInt<4> = BigInt::<4>::ZERO();
    // const R: BigInt<4> = BigInt::<4>::ZERO();
    // const R2: BigInt<4> = BigInt::<4>::ZERO();
    // const R3: BigInt<4> = BigInt::<4>::ZERO();
    // const M0: u64 = 165 as u64;
    // const E: u64 = 8 as u64;
    // const RODD: u64 = 13 as u64;
    // const N: BigInt<4> = BigInt::<4>::ZERO();

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
        for i in 0..s {
            // t = t + self * other[i]
            let (mut tmp_c1, mut tmp_c2, mut ab) = (0 as u64, 0 as u64, BigInt([0 as u64; 4]));
            (ab, tmp_c1) = lft.clone() * rht.0[i];
            (t, tmp_c2) = t + ab;
            (c1, overflow) = c1.overflowing_add(tmp_c1.wrapping_add(tmp_c2));
            if overflow {
                c2 += 1_u64;
            }

            // t = t + ((t[0] * N'[0]) mod W) * N
            let (mut tmp_c3, mut tmp_c4, mut mn) = (0 as u64, 0 as u64, BigInt([0 as u64; 4]));
            let m = Self::M0.wrapping_mul(t.0[0]);
            (mn, tmp_c3) = Self::MODULUS * m;
            (t, tmp_c4) = t + mn;
            (c1, overflow) = c1.overflowing_add(tmp_c3.wrapping_add(tmp_c4));
            if overflow {
                c2 += 1_u64;
            }

            // t >> 1
            for j in 0..(s - 1) {
                t.0[j] = t.0[j + 1];
            }
            (t.0[s - 1], c1, c2) = (c1, c2, 0 as u64);
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
        unimplemented!()
    }
}
