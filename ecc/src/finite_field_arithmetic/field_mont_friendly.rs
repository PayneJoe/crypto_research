use crate::finite_field_arithmetic::traits::Field;
use std::cmp::Ordering;
use std::iter::Sum;
use std::ops::{Add, Div, Mul, Neg, Rem, Shl, Shr, Sub};
use std::str::FromStr;

//////////////////////////////////// Implementation of BigInteger Specially for Finite Field
#[derive(Debug, PartialEq, Clone, Copy, Eq)]
pub struct BI<const N: usize>(pub [u8; N]);

impl PartialOrd for BI<2> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for BI<2> {
    fn cmp(&self, other: &Self) -> Ordering {
        let n = 2;
        let mut ord = Ordering::Equal;
        for i in (0..n).rev() {
            ord = self.0[i].cmp(&other.0[i]);
            if ord != Ordering::Equal {
                return ord;
            }
        }
        ord
    }
}

impl BI<2> {
    pub fn to_bits(self) -> Vec<u8> {
        let num_leading_zeros = self.0.iter().rev().fold(0, |acc, n| {
            if acc == 0 {
                acc + n.leading_zeros()
            } else if (acc % 8) == 0 {
                acc + n.leading_zeros()
            } else {
                acc
            }
        });
        let bits = (0..self.0.len())
            .map(|i| ((0..8).map(|j| ((self.0[i] >> j) & 1) as u8)).collect::<Vec<u8>>())
            .collect::<Vec<Vec<u8>>>()
            .into_iter()
            .flatten()
            .collect::<Vec<u8>>();
        let end = bits.len() - (num_leading_zeros as usize);
        bits[..end].to_vec()
    }

    fn bit_size(self) -> usize {
        let mut num_bits = 0 as usize;
        let n = 2;
        let (mut a, mut b) = (0 as usize, 0 as usize);
        let mut found = false;
        for i in (0..n).rev() {
            if self.0[i] > 0 {
                for j in (0..8).rev() {
                    if self.0[i] & (1 << j) > 0 {
                        (a, b) = (i, j);
                        found = true;
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        a * 8 + b + 1
    }

    #[inline(always)]
    fn MAX() -> Self {
        BI(vec![((1 << 8) - 1) as u8; 2].try_into().unwrap())
    }

    #[inline(always)]
    fn zero() -> Self {
        BI([0, 0])
    }

    #[inline(always)]
    fn one() -> Self {
        BI([1, 0])
    }

    #[inline(always)]
    fn is_even(self) -> bool {
        self.0[0] % 2 == 0
    }

    #[inline(always)]
    fn is_zero(self) -> bool {
        let n = self.0.len();
        let mut zero = true;
        for i in 0..n {
            zero = zero & (self.0[i] == 0);
        }
        zero
    }
}

impl FromStr for BI<2> {
    type Err = ParseStrErr;

    fn from_str(text: &str) -> Result<Self, Self::Err> {
        // parse the large number
        let mut large_number = text.parse::<u16>().map_err(|_| ParseStrErr).unwrap();

        let apply_parse = |word: u16| -> BI<2> {
            let mut large_word = word;
            let num_bits = 8 as usize;
            let mut w = vec![0_u8; 2];
            let mask = ((1 << num_bits) - 1) as u16;
            for i in 0..2 {
                w[i] = (large_word & mask) as u8;
                large_word = large_word >> 8;
            }
            BI(w.try_into().unwrap())
        };

        // apply reduce through mul_reduce
        Ok(apply_parse(large_number))
    }
}

// naive implementation of modular for bigint
impl Rem<Self> for BI<2> {
    type Output = BI<2>;
    fn rem(self, other: Self) -> Self::Output {
        if self < other {
            return other;
        }
        let (asize, bsize) = (self.bit_size(), other.bit_size());
        assert!(asize >= bsize);
        let mut prox_lft = (other << (asize - bsize) as u8).0;
        if prox_lft > self {
            prox_lft = (prox_lft >> 1).0;
        }
        let mut remainder = (self - prox_lft).0;
        while remainder >= other {
            remainder = (remainder - other).0;
        }
        remainder
    }
}

impl Shl<u8> for BI<2> {
    type Output = (BI<2>, u8);
    fn shl(self, other: u8) -> Self::Output {
        assert!(other < 8);
        let n = 2;
        let mut w = Self([0_u8; 2]);
        let mut carrier = 0_u8;
        // from the lowest word to the highest word
        for i in 0..n {
            let tmp_carrier = self.0[i] >> (8 - other);
            let remainder = (self.0[i] << other) | carrier;
            w.0[i] = remainder;
            carrier = tmp_carrier;
        }
        (w, carrier)
    }
}

impl Shr<u8> for BI<2> {
    type Output = (BI<2>, u8);
    fn shr(self, other: u8) -> Self::Output {
        assert!(other < 8);
        let n = 2;
        let mut w = Self([0_u8; 2]);
        let mut carrier = 0_u8;
        // from the lowest word to the highest word
        for i in (0..n).rev() {
            let tmp_carrier = self.0[i] & (((1 << other) - 1) as u8);
            let remainder = (self.0[i] >> other) | (carrier << (8 - other));
            w.0[i] = remainder;
            carrier = tmp_carrier;
        }
        (w, carrier)
    }
}

impl Add for BI<2> {
    type Output = (BI<2>, u8);
    fn add(self, other: BI<2>) -> Self::Output {
        let n = self.0.len();
        let lft = &self.0;
        let rht = &other.0;
        let (mut carrier, mut remainder) = (0_u8, 0_u8);
        let mut w = BI([0_u8; 2]);
        for i in 0..n {
            let (lft_w, rht_w, carrier_w) =
                (u16::from(lft[i]), u16::from(rht[i]), u16::from(carrier));
            let t = lft_w + (rht_w + carrier_w);
            (carrier, remainder) = ((t >> 8) as u8, (t & ((1 << 8) - 1)) as u8);
            w.0[i] = remainder;
        }
        (w, carrier)
    }
}

impl<'a, 'b> Add<&'b BI<2>> for &'a BI<2> {
    type Output = (BI<2>, u8);

    fn add(self, other: &'b BI<2>) -> Self::Output {
        unimplemented!()
    }
}

impl Sub for BI<2> {
    type Output = (BI<2>, u8);
    fn sub(self, other: BI<2>) -> Self::Output {
        let n = self.0.len();
        let lft = &self.0;
        let rht = &other.0;
        let (mut carrier, mut remainder) = (0_u8, 0_u8);
        let mut w = BI([0_u8; 2]);
        for i in 0..n {
            let (lft_w, rht_w, carrier_w) =
                (u16::from(lft[i]), u16::from(rht[i]), u16::from(carrier));
            (carrier, remainder) = if lft_w >= rht_w + carrier_w {
                (0_u8, (lft_w - (rht_w + carrier_w)) as u8)
            } else {
                (1_u8, (lft_w + (1 << 8) - (rht_w + carrier_w)) as u8)
            };
            w.0[i] = remainder;
        }
        (w, carrier)
    }
}

impl<'a, 'b> Sub<&'b BI<2>> for &'a BI<2> {
    type Output = (BI<2>, u8);
    fn sub(self, other: &'b BI<2>) -> Self::Output {
        unimplemented!()
    }
}

impl Mul for BI<2> {
    type Output = (BI<2>, u8);
    fn mul(self, other: BI<2>) -> Self::Output {
        unimplemented!()
    }
}

impl Mul<u8> for BI<2> {
    type Output = (BI<2>, u8);

    fn mul(self, other: u8) -> Self::Output {
        let n = self.0.len();
        let lft = &self.0;
        let (mut carrier, mut remainder) = (0_u8, 0_u8);
        let mut w = BI([0_u8; 2]);
        for i in 0..n {
            let (lft_w, rht_w, carrier_w) =
                (u16::from(lft[i]), u16::from(other), u16::from(carrier));
            let t = lft_w * rht_w + carrier_w;
            (carrier, remainder) = ((t >> 8) as u8, (t & ((1 << 8) - 1)) as u8);
            w.0[i] = remainder;
        }
        (w, carrier)
    }
}

impl<'a, 'b> Mul<&'b BI<2>> for &'a BI<2> {
    type Output = (BI<2>, u8);
    fn mul(self, other: &'b BI<2>) -> Self::Output {
        unimplemented!()
    }
}

impl Div for BI<2> {
    type Output = (BI<2>, u8);
    fn div(self, other: BI<2>) -> Self::Output {
        unimplemented!()
    }
}

impl<'a, 'b> Div<&'b BI<2>> for &'a BI<2> {
    type Output = (BI<2>, u8);
    fn div(self, other: &'b BI<2>) -> Self::Output {
        unimplemented!()
    }
}

//////////////////////////////////////// Custom Finite Field
#[derive(Debug)]
pub struct ParseStrErr;

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub struct Foo<const T: usize>(pub BI<T>);

/////////////////////////////////////// Implementation of Custom Finite Field
impl FromStr for Foo<2> {
    type Err = ParseStrErr;

    fn from_str(text: &str) -> Result<Self, Self::Err> {
        // parse the large number
        let mut large_number = text.parse::<u32>().map_err(|_| ParseStrErr).unwrap();

        let apply_parse = |word: u16| -> BI<2> {
            let mut large_word = word;
            let num_bits = 8 as usize;
            let mut w = vec![0_u8; 2];
            let mask = ((1 << num_bits) - 1) as u16;
            for i in 0..2 {
                w[i] = (large_word & mask) as u8;
                large_word = large_word >> 8;
            }
            BI(w.try_into().unwrap())
        };

        let (high_words, low_words) = (
            (large_number >> 16) as u16,
            (large_number & ((1 << 16) - 1)) as u16,
        );

        // apply reduce through mul_reduce
        let (high_bi, low_bi) = (apply_parse(high_words), apply_parse(low_words));
        if high_bi == BI::zero() {
            Ok(Self::reduce(&low_bi, Some(false)))
        } else {
            Ok(Self::reduce(&high_bi, Some(true)) + Self::reduce(&low_bi, Some(false)))
        }
    }
}

impl From<&BI<2>> for Foo<2> {
    fn from(value: &BI<2>) -> Self {
        Self::reduce(value, Some(false))
    }
}

impl From<BI<2>> for Foo<2> {
    fn from(value: BI<2>) -> Self {
        Self::reduce(&value, Some(false))
    }
}

impl Into<BI<2>> for Foo<2> {
    fn into(self) -> BI<2> {
        self.rev_reduce()
    }
}

impl<'a> Sum<&'a Self> for Foo<2> {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = &'a Self>,
    {
        iter.fold(Foo::<2>::ZERO(), |a, b| a + *b)
    }
}

impl Sum<Self> for Foo<2> {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        iter.fold(Foo::<2>::ZERO(), |a, b| a + b)
    }
}

impl Add for Foo<2> {
    type Output = Foo<2>;

    fn add(self, other: Self) -> Foo<2> {
        let (w, carrier) = self.0 + other.0;
        if carrier > 0 {
            Self(((BI::<2>::MAX() - Self::MODULUS).0 + w).0)
        } else if w >= Self::MODULUS {
            Self((w - Self::MODULUS).0)
        } else {
            Self(w)
        }
    }
}

impl Sub for Foo<2> {
    type Output = Foo<2>;

    fn sub(self, other: Self) -> Foo<2> {
        if self.0 > other.0 {
            Self((self.0 - other.0).0)
        } else {
            Self((Self::MODULUS - (other.0 - self.0).0).0)
        }
    }
}

impl Mul for Foo<2> {
    type Output = Foo<2>;

    fn mul(self, other: Self) -> Foo<2> {
        Self::mul_reduce(&self.0, &other.0)
    }
}

impl Mul<u8> for Foo<2> {
    type Output = Foo<2>;

    fn mul(self, other: u8) -> Foo<2> {
        Self::mul_reduce(
            &self.0,
            &Self::from_str(other.to_string().as_str()).unwrap().0,
        )
    }
}

impl Div for Foo<2> {
    type Output = Foo<2>;

    fn div(self, other: Self) -> Foo<2> {
        Self::mul_reduce(&self.0, &other.inv().0)
    }
}

impl Neg for Foo<2> {
    type Output = Foo<2>;

    fn neg(self) -> Foo<2> {
        Self::ZERO() - self
    }
}

pub trait Exponentiation {
    type Output;

    // outdated interface
    fn exp(&self, n: BI<2>) -> Self;

    fn exp_raw(&self, n: &[u8]) -> Self;
}

impl Exponentiation for Foo<2> {
    type Output = Foo<2>;

    // referenced from Algorithm 11.7 of "handbook of elliptic and hyperelliptic curve cryptography"
    fn exp(&self, n: BI<2>) -> Self {
        let n_bits = n.to_bits();
        self.exp_raw(n_bits.as_slice())
    }

    fn exp_raw(&self, n: &[u8]) -> Self {
        let n_bits = n;
        let (mut y, x) = (Self::ONE().0, self.0);
        for i in (0..n_bits.len()).rev() {
            y = Self::mul_reduce(&y, &y).0;
            if n_bits[i] == 1 {
                y = Self::mul_reduce(&x, &y).0;
            }
        }
        Self(y)
    }
}

pub trait Sqrt: Sized {
    type Output;

    fn sqrt(self) -> Option<Self>;
}

impl Sqrt for Foo<2> {
    type Output = Option<Foo<2>>;

    // referenced from Lemma 11.22 of "handbook of elliptic and hyperelliptic curve cryptography"
    fn sqrt(self) -> Self::Output {
        // MODULUS % 4 == 3
        let mut n = BI([0, 0]);
        if Self::MODULUS.0[0] & (4_u8 - 1) == 3 {
            n = ((Self::MODULUS + BI([1, 0])).0 >> 2).0;
            return Some(self.exp(n));
        }
        // MODULUS % 8 == 5 and a^{(p - 1)/4} == 1
        n = ((Self::MODULUS - BI([1, 0])).0 >> 2).0;
        if (Self::MODULUS.0[0] & (8_u8 - 1) == 5) && (self.exp(n) == Self::ONE()) {
            n = ((Self::MODULUS + BI([3, 0])).0 >> 3).0;
            return Some(self.exp(n));
        }
        // MODULUS % 8 == 5 and a^{(p - 1)/4} == -1
        if (Self::MODULUS.0[0] & (8_u8 - 1) == 5) && (self.exp(n) == -Self::ONE()) {
            n = ((Self::MODULUS - BI([5, 0])).0 >> 3).0;
            let a_4times = self * 4;
            return Some(a_4times.exp(n) * self * 2);
        }
        return None;
    }
}

// define custom finite field
impl Field<2> for Foo<2> {
    // fabricated precomputable parameters of custom finite field,
    // these constant parameters need to determined at compile time

    // W = 256, MODULUS/M = 2003, R = 256^2 % MODULUS = 982, M0 = (-M^{-1} % 256^2) % W =
    // MODULUS % 8 \noteq 1
    const MODULUS: BI<2> = BI([211, 7]);
    const R: BI<2> = BI([160, 5]);
    const R2: BI<2> = BI([239, 1]);
    const R3: BI<2> = BI([199, 6]);
    const M0: u8 = 165_u8;

    fn legendre_symbol(self) -> bool {
        let mut k = 1;
        let (mut p, mut a) = (Self::MODULUS.clone(), self.0.clone());
        while p != BI::<2>::one() {
            if a.is_zero() {
                return false;
            }
            let mut v = 0;
            while a.is_even() {
                a = (a >> 1).0;
                v = v + 1;
            }
            if (v % 2 == 1) && ((p.0[0] % 8 == 3) || (p.0[0] % 8 == 5)) {
                k = -k;
            }
            if (a.0[0] % 4 == 3) && (p.0[0] % 4 == 3) {
                k = -k;
            }
            // time-consuming
            (a, p) = (p % a, a)
        }
        k == 1
    }

    fn is_quadratic_residual(self) -> bool {
        self.exp(((Self::MODULUS - BI::<2>::one()).0 >> 1).0) == Self::ONE()
    }

    fn square_inplace(&mut self) {
        *self = Self::mul_reduce(&(*self).0, &(*self).0);
    }

    fn square(&self) -> Self {
        Self::mul_reduce(&self.0, &self.0)
    }

    fn reduce(u: &BI<2>, inv: Option<bool>) -> Self {
        if let Some(true) = inv {
            Self::mul_reduce(u, &Self::R3)
        } else {
            Self::mul_reduce(u, &Self::R2)
        }
    }

    #[inline(always)]
    fn ONE() -> Self {
        Self(Self::R)
    }

    #[inline(always)]
    fn ZERO() -> Self {
        Self(BI::<2>::zero())
    }

    // F^e
    fn pow(&self, e: BI<2>) -> Self {
        let n_bits = e.to_bits();
        self.exp_raw(n_bits.as_slice())
    }

    // referenced from Algorithm 11.12 of "handbook of elliptic and hyperelliptic curve cryptography"
    // (aR)^{-1} % N <- ((aR)^{-1} * R^2) * R^{-1} % N
    fn inv(&self) -> Self {
        let (mut r, mut s, mut t, mut v) = (self.0, BI::<2>::one(), Self::MODULUS, BI::<2>::zero());

        /////////////////////////////////////////////// STAGE ONE
        let mut k = 0 as usize;
        let mut tmp_carrier = 0_u8;
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
        let m = 2 * 8 as usize;
        // println!("------ k = {}, m = {}, v = {:?}", k, m, v);
        if k < m {
            (v, k) = (Self::mul_reduce(&v, &Self::R2).0, k + m);
        }

        // 2^{2m - k}
        let mut h_bit = 2 * m - k;
        let x: Vec<u8> = (0..self.0 .0.len())
            .map(|i| {
                if h_bit >= 8 {
                    h_bit -= 8;
                    0_u8
                } else {
                    let tmp_bit = h_bit;
                    h_bit = 0;
                    if tmp_bit > 0 {
                        1 << tmp_bit
                    } else {
                        0_u8
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
        let res = BI(x.try_into().unwrap());
        if res.is_zero() == false {
            v = Self::mul_reduce(&v, &res).0;
        }

        Self(v)
    }

    // a % N <- aR * R^{-1} % N
    fn rev_reduce(&self) -> BI<2> {
        Self::mul_reduce(&self.0, &BI::<2>::one()).0
    }

    // referenced from Algorithm 11.3 of "handbook of elliptic and hyperelliptic curve cryptography"
    // abR % N <- (aR * bR) * R^{-1} % N
    fn mul_reduce(lft: &BI<2>, rht: &BI<2>) -> Self {
        if (*lft == BI::zero()) || (*rht == BI::zero()) {
            return Self::ZERO();
        }
        let s = 2;
        let mut t = BI([0_u8; 2]);
        let (mut c1, mut c2, mut overflow) = (0_u8, 0_u8, false);
        for i in 0..s {
            // t = t + self * other[i]
            let (mut tmp_c1, mut tmp_c2, mut ab) = (0_u8, 0_u8, BI([0_u8; 2]));
            (ab, tmp_c1) = lft.clone() * rht.0[i];
            (t, tmp_c2) = t + ab;
            (c1, overflow) = c1.overflowing_add(tmp_c1.wrapping_add(tmp_c2));
            if overflow {
                c2 += 1_u8;
            }

            // t = t + ((t[0] * N'[0]) mod W) * N
            let (mut tmp_c3, mut tmp_c4, mut mn) = (0_u8, 0_u8, BI([0_u8; 2]));
            let m = Self::M0.wrapping_mul(t.0[0]);
            (mn, tmp_c3) = Self::MODULUS * m;
            (t, tmp_c4) = t + mn;
            (c1, overflow) = c1.overflowing_add(tmp_c3.wrapping_add(tmp_c4));
            if overflow {
                c2 += 1_u8;
            }

            // t >> 1
            for j in 0..(s - 1) {
                t.0[j] = t.0[j + 1];
            }
            (t.0[s - 1], c1, c2) = (c1, c2, 0_u8);
        }

        Self(t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // referenced from Algorithm 10.42 of "handbook of elliptic and hyperelliptic curve cryptography"
    // GCD of two single precision positive integers
    fn gcd(a: u32, b: u32) -> (u32, u32, u32, bool) {
        assert!(a < b);
        let (mut A, mut B) = (b, a);
        let (mut Ua, mut Ub) = (0, 1);
        let (mut Va, mut Vb) = (1, 0);
        let mut n_iter = 0_u16;
        while B != 0 {
            let q = A / B;
            (A, B) = (B, A - q * B);
            (Ua, Ub) = (Ub, Ua + q * Ub);
            (Va, Vb) = (Vb, Va + q * Vb);
            n_iter = n_iter + 1;
        }
        let (d, u, v) = (A, Ua, Va);
        (u, v, d, n_iter % 2 == 0)
    }

    #[test]
    fn test_gcd() {
        let (a, b) = (2003, 65536);
        let (u, v, d, sign) = gcd(a, b);
        println!("u = {}, v = {}, d = {}, sign = {}", u, v, d, sign);
        if sign {
            assert_eq!(b * v - a * u, d)
        } else {
            assert_eq!(a * u - b * v, d);
        }
    }

    #[test]
    fn test_inv_sqrt() {
        let (a, c) = (527_u32, 899_u32);
        let lft = Foo::<2>::from_str(a.to_string().as_str()).unwrap();
        let result = Foo::<2>::from_str(c.to_string().as_str()).unwrap();
        let lft_inv = lft.inv();
        assert_eq!(lft_inv.is_quadratic_residual(), true);
        // assert_eq!(lft_inv.sqrt().unwrap(), result);
    }
}
