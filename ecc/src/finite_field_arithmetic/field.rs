use std::cmp::Ordering;
use std::ops::{Add, Div, Mul, Rem, Shl, Shr, Sub};
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
    fn to_bits(self) -> Vec<u8> {
        (0..self.0.len())
            .map(|i| ((0..8).map(|j| ((self.0[i] >> j) & 1) as u8)).collect::<Vec<u8>>())
            .collect::<Vec<Vec<u8>>>()
            .into_iter()
            .flatten()
            .collect::<Vec<u8>>()
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

////////////////////////////////// Constant configurations of field including some core behavious on field
pub trait Field<const N: usize>: FromStr + From<BI<N>> + Into<BI<N>> {
    // finite field modulus
    const MODULUS: BI<N>;
    // W^s % MODULUS, for Montgomery reduce
    const R: BI<N>;
    // W^2s % MODULUS, for Montgomery reduce
    const R2: BI<N>;
    // W^3s % MODULUS, for Montgomery reduce
    const R3: BI<N>;
    // inversion of least significant word of modulus, also convenient for Montgomery reduce
    const M0: u8;
    // for the convenient of square root
    const E: u8;
    const RODD: u8;
    const N: BI<N>;

    fn ONE() -> Self;
    fn ZERO() -> Self;

    // legendre symbol
    fn legendre_symbol(self) -> bool;

    fn is_quadratic_residual(self) -> bool;

    // reduce a bigint into the specified range [0, modulus)
    fn reduce(u: &BI<N>, inv: Option<bool>) -> Self;
    // convert a reduced number into a unreduced bigint
    fn rev_reduce(&self) -> BI<N>;
    // inversion of a reduced number
    fn inv(&self) -> Self;

    //
    fn square_inplace(&mut self);

    // Montgomery Reduction
    fn mul_reduce(lft: &BI<N>, rht: &BI<N>) -> Self;
}

//////////////////////////////////////// Custom Finite Field
#[derive(Debug)]
pub struct ParseStrErr;

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Foo<const T: usize>(BI<T>);

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
            Self((other.0 - self.0).0)
        }
    }
}

impl Mul for Foo<2> {
    type Output = Foo<2>;

    fn mul(self, other: Self) -> Foo<2> {
        Self::mul_reduce(&self.0, &other.0)
    }
}

impl Div for Foo<2> {
    type Output = Foo<2>;

    fn div(self, other: Self) -> Foo<2> {
        Self::mul_reduce(&self.0, &other.inv().0)
    }
}

pub trait Exponentiation {
    type Output;

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

    // referenced from Algorithm 11.23 of "handbook of elliptic and hyperelliptic curve cryptography"
    fn sqrt(self) -> Self::Output {
        let (n, r) = (Self::reduce(&Self::N, Some(false)), BI([Self::RODD, 0_u8]));
        let (a, z) = (self.clone(), n.exp(r));
        let (mut y, mut s, mut x, mut b) = (
            z,
            Self::E,
            a.exp(((r - BI::<2>::one()).0 >> 1).0),
            Self::ZERO(),
        );
        let tmp_x = a * x.clone();
        (b, x) = (tmp_x.clone() * x, tmp_x);
        while b != Self::ONE() {
            let mut m = 0 as u8;
            let mut tmp_b = b.clone();
            while tmp_b != Self::ONE() {
                tmp_b.square_inplace();
                m = m + 1;
            }
            // non-quadratic residual
            if m == Self::E {
                return None;
            }
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
            assert!(x * x == self);
        }
        Some(x)
    }
}

// define custom finite field
impl Field<2> for Foo<2> {
    // fabricated precomputable parameters of custom finite field,
    // these constant parameters need to determined at compile time
    // W = 256, MODULUS/M = 517, R = 256^2 % M = 394, M0 = -M[0]^{-1} % W = 51
    // const MODULUS: BI<2> = BI([5, 2]);
    // const R: BI<2> = BI([138, 1]);
    // const R2: BI<2> = BI([136, 0]);
    // const R3: BI<2> = BI([77, 1]);
    // const M0: u8 = 51_u8;

    // W = 256, MODULUS/M = 3329 = 2^8 * 13 + 1, R = 256^2 % MODULUS = 982, M0 = -M[0]^{-1} % MODULUS = 1
    const MODULUS: BI<2> = BI([1, 13]);
    const R: BI<2> = BI([237, 8]);
    const R2: BI<2> = BI([73, 5]);
    const R3: BI<2> = BI([245, 8]);
    const M0: u8 = 255_u8;
    // 2^E * RODD = MODULUS - 1
    const E: u8 = 8_u8;
    const RODD: u8 = 13_u8;
    // N is a sampled non-quadratic residual number, 3/6/11/12/...
    const N: BI<2> = BI([3, 0]);

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

    fn reduce(u: &BI<2>, inv: Option<bool>) -> Self {
        if let Some(true) = inv {
            Self::mul_reduce(u, &Self::R3)
        } else {
            Self::mul_reduce(u, &Self::R2)
        }
    }

    fn ONE() -> Self {
        Self(Self::R)
    }

    fn ZERO() -> Self {
        Self(BI::<2>::zero())
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
        // println!("-----------*****---------- First stage finished!");

        //////////////////////////////////////////// STAGE TWO
        let m = 2 * 8 as usize;
        if k < m {
            (v, k) = (Self::mul_reduce(&v, &Self::R2).0, k + m);
        }
        let h_bit = 2 * m - k;

        // 2^{2m - k}
        let mut x: Vec<u8> = vec![0_u8; h_bit / 8]
            .to_vec()
            .into_iter()
            .chain([(1 << (h_bit % 8)) as u8].into_iter())
            .collect();
        let tmp_len = x.len();
        if tmp_len < self.0 .0.len() {
            x = x
                .into_iter()
                .chain(vec![0_u8; self.0 .0.len() - tmp_len].into_iter())
                .collect();
        }
        println!("------- h_bit = {}, x = {:?}, v = {:?}", h_bit, x, v);

        // REDC(v * R2)
        v = Self::mul_reduce(&v, &Self::R2).0;
        // REDC(v * 2^{2m - k})
        v = Self::mul_reduce(&v, &BI(x.try_into().unwrap())).0;

        Self(v)
    }

    // a % N <- aR * R^{-1} % N
    fn rev_reduce(&self) -> BI<2> {
        Self::mul_reduce(&self.0, &BI::<2>::one()).0
    }

    // referenced from Algorithm 11.3 of "handbook of elliptic and hyperelliptic curve cryptography"
    // abR % N <- (aR * bR) * R^{-1} % N
    fn mul_reduce(lft: &BI<2>, rht: &BI<2>) -> Self {
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
        let (a, b) = (3329, 65536);
        let (u, v, d, sign) = gcd(a, b);
        println!("u = {}, v = {}, d = {}, sign = {}", u, v, d, sign);
        if sign {
            assert_eq!(b * v - a * u, d)
        } else {
            assert_eq!(a * u - b * v, d);
        }
    }

    #[test]
    fn test_shift() {
        let a = BI([64_u8, 1_u8]);
        assert_eq!((a >> 1).0, BI([160_u8, 0_u8]));
        let b = BI([129_u8, 2_u8]);
        assert_eq!((b << 1).0, BI([2_u8, 5_u8]));
    }

    #[test]
    fn test_bitsize() {
        let b = BI([129_u8, 2_u8]);
        assert_eq!(b.bit_size(), 10);
    }

    #[test]
    fn test_fromstr() {
        let a = 259_u16;
        let result: BI<2> = Foo::<2>::from_str(a.to_string().as_str()).unwrap().into();
        let actual = BI([(a % 256) as u8, (a / 256) as u8]);
        assert_eq!(result, actual);
    }

    #[test]
    fn test_reduce() {
        let (a, b) = (4 as u16, ((4 * 256 * 256) % 3329) as u16);
        let lft = Foo::<2>::from_str(a.to_string().as_str()).unwrap();
        // let rht = Foo::<2>::from_str(b.to_string().as_str()).unwrap();
        let rht = Foo(BI::<2>::from_str(b.to_string().as_str()).unwrap());
        assert_eq!(lft, rht);
    }

    #[test]
    fn test_add() {
        let (a, b) = (259_u16, 258_u16);
        let c = u32::from(a) + u32::from(b);
        let lft = Foo::<2>::from_str(a.to_string().as_str()).unwrap();
        let rht = Foo::<2>::from_str(b.to_string().as_str()).unwrap();
        let result = Foo::<2>::from_str(c.to_string().as_str()).unwrap();
        assert_eq!(lft + rht, result);
    }

    #[test]
    fn test_mul() {
        let (a, b) = (259_u16, 258_u16);
        let c = u32::from(a) * u32::from(b);
        let lft = Foo::<2>::from_str(a.to_string().as_str()).unwrap();
        let rht = Foo::<2>::from_str(b.to_string().as_str()).unwrap();
        let result = Foo::<2>::from_str(c.to_string().as_str()).unwrap();
        assert_eq!(lft * rht, result);
    }

    #[test]
    fn test_inv() {
        let (a, M) = (259_u32, 3329_u32);
        let (mut c, _, d, sign) = gcd(a as u32, M as u32);
        assert!(d == 1);
        c = if sign { M - c } else { c };
        let lft = Foo::<2>::from_str(a.to_string().as_str()).unwrap();
        let result = Foo::<2>::from_str(c.to_string().as_str()).unwrap();
        assert_eq!(Foo::<2>::mul_reduce(&lft.0, &result.0), Foo::<2>::ONE());
        assert_eq!(lft.inv(), result);
    }

    #[test]
    fn test_exp() {
        let (a, b) = (259_u32, 3_u32);
        let (c, n) = (
            a.pow(b.into()),
            BI::<2>::from_str(b.to_string().as_str()).unwrap(),
        );
        let lft = Foo::<2>::from_str(a.to_string().as_str()).unwrap();
        let result = Foo::<2>::from_str(c.to_string().as_str()).unwrap();
        assert_eq!(lft.exp(n), result);
    }

    #[test]
    fn test_is_quadratic_residual() {
        let residual_numbers: Vec<i32> = (1..3329)
            .map(|i| {
                let lft = Foo::<2>::from_str(i.to_string().as_str()).unwrap();
                let e = ((Foo::<2>::MODULUS - BI::<2>::one()).0 >> 1).0;
                if lft.exp(e) == Foo::<2>::ONE() {
                    i
                } else {
                    0
                }
            })
            .filter(|x| x > &0)
            .collect();
        println!(
            "samples of quadratic residual number: {:?}",
            residual_numbers[..200].to_vec()
        );
        assert_eq!(residual_numbers.len(), (3329 - 1) / 2);
    }

    #[test]
    fn test_square_root() {
        let (a, b) = (180_u32, 181_u32);
        let pos = Foo::<2>::from_str(a.to_string().as_str()).unwrap();
        let neg = Foo::<2>::from_str(b.to_string().as_str()).unwrap();
        if pos.is_quadratic_residual() {
            let result = pos.sqrt().unwrap();
            assert_eq!(result * result, pos);
        }
        if neg.is_quadratic_residual() {
            assert_eq!(neg.sqrt(), None);
        }
    }
}
