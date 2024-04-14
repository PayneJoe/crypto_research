///// Practical Implementation for Base Field (Fr) of BLS12-381
// upon which the largest prime factor of E(Fq) defined
use crate::finite_field_arithmetic::bigint::BigInt;
use crate::finite_field_arithmetic::pairing_friendly::field::{Field, LegendreSymbol, PrimeField};

use std::iter;
use std::iter::Sum;
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::str::FromStr;

#[derive(Debug, PartialEq, Eq)]
pub struct FieldParseErr;

#[derive(Debug, PartialEq, Eq)]
pub struct NoneQuadraticResidualErr;

#[derive(Debug, PartialEq, Eq, Clone, Copy, Hash)]
pub struct Fr<const N: usize>(pub BigInt<N>);

// since the bit length of Fr is 381, so we need 6 words to hold
const NUM_LIMBS: usize = 4;
const WORD_SIZE: usize = 64;
const BYTE_SIZE: usize = NUM_LIMBS * WORD_SIZE;
type Word = u64;
type Byte = u8;

// MODULUS = 52435875175126190479447740508185965837690552500527637822603658699938581184513
// W = 2^64 = 18446744073709551616
// r = W^4 = 2^256 = 115792089237316195423570985008687907853269984665640564039457584007913129639936
// R = r % MODULUS = 10920338887063814464675503992315976177888879664585288394250266608035967270910
// R2 = r * r % MODULUS = 3294906474794265442129797520630710739278575682199800681788903916070560242797
// R3 = r * r * r % MODULUS = 49829253988540319354550742249276084460127446355315915089527227471280320770991
//
// M0 = (-MODULUS^{-1} % r) % W = 18446744069414584319
//
// let: MODULUS = 2^E * ROOD + 1
// E = 32
// RODD = 12208678567578594777604504606729831043093128246378069236549469339647
// N = 5/7/...., is sampled number whose legendre symbol is -1 (is definitely non-quadratic residual)

impl PrimeField<NUM_LIMBS> for Fr<NUM_LIMBS> {
    const MODULUS: BigInt<NUM_LIMBS> = BigInt([
        18446744069414584321,
        6034159408538082302,
        3691218898639771653,
        8353516859464449352,
    ]);
    const R: BigInt<NUM_LIMBS> = BigInt([
        8589934590,
        6378425256633387010,
        11064306276430008309,
        1739710354780652911,
    ]);
    const R2: BigInt<NUM_LIMBS> = BigInt([
        14526898881837571181,
        3129137299524312099,
        419701826671360399,
        524908885293268753,
    ]);
    const R3: BigInt<NUM_LIMBS> = BigInt([
        14279814937963099055,
        1963020886675057040,
        8345518043873801240,
        7938258146690806761,
    ]);
    const M0: Word = 18446744069414584319 as Word;
    const E: Word = 32 as Word;
    const RODD: BigInt<NUM_LIMBS> = BigInt([
        18446282274530918399,
        694073334983140354,
        2998690675949164552,
        1944954707,
    ]);
    const N: BigInt<NUM_LIMBS> = BigInt([5, 0, 0, 0]);

    // a % N <- aR * R^{-1} % N
    fn rev_reduce(&self) -> BigInt<NUM_LIMBS> {
        Self::mul_reduce(&self.0, &BigInt::<NUM_LIMBS>::ONE()).0
    }

    // referenced from Algorithm 11.3 of "handbook of elliptic and hyperelliptic curve cryptography"
    // abR % N <- (aR * bR) * R^{-1} % N
    fn mul_reduce(lft: &BigInt<NUM_LIMBS>, rht: &BigInt<NUM_LIMBS>) -> Self {
        if (*lft == BigInt::ZERO()) || (*rht == BigInt::ZERO()) {
            return Self::ZERO();
        }
        let mut t = BigInt([0 as Word; NUM_LIMBS]);
        let (mut c1, mut c2, mut overflow) = (0 as Word, 0 as Word, false);
        // println!("lft = {:?}, rht = {:?}", lft.0, rht.0);
        for i in 0..NUM_LIMBS {
            // t = t + self * other[i]
            let (mut tmp_c1, mut tmp_c2, mut ab) =
                (0 as Word, 0 as Word, BigInt([0 as Word; NUM_LIMBS]));
            (ab, tmp_c1) = lft.clone() * rht.0[i];
            (t, tmp_c2) = t + ab;
            (c1, overflow) = c1.overflowing_add(tmp_c1.wrapping_add(tmp_c2));
            if overflow {
                c2 += 1 as Word;
            }
            // println!("i = {}, #1: t = {:?}, c1 = {:?}, c2 = {:?}", i, t.0, c1, c2);

            // t = t + ((t[0] * N'[0]) mod W) * N
            let (mut tmp_c3, mut tmp_c4, mut mn) =
                (0 as Word, 0 as Word, BigInt([0 as Word; NUM_LIMBS]));
            let m = Self::M0.wrapping_mul(t.0[0]);
            (mn, tmp_c3) = Self::MODULUS * m;
            (t, tmp_c4) = t + mn;
            (c1, overflow) = c1.overflowing_add(tmp_c3.wrapping_add(tmp_c4));
            if overflow {
                c2 += 1 as Word;
            }
            // println!("i = {} #2: t = {:?}, c1 = {:?}, c2 = {:?}", i, t.0, c1, c2);

            // t >> 1
            for j in 0..(NUM_LIMBS - 1) {
                t.0[j] = t.0[j + 1];
            }
            (t.0[NUM_LIMBS - 1], c1, c2) = (c1, c2, 0 as Word);
            // println!("i = {} #3: t = {:?}, c1 = {:?}, c2 = {:?}", i, t.0, c1, c2);
        }
        // 1/30/2024 fixed
        if t > Self::MODULUS {
            t = (t - Self::MODULUS).0;
        }

        Self(t)
    }

    fn reduce(u: &BigInt<NUM_LIMBS>, inv: Option<bool>) -> Self {
        if let Some(true) = inv {
            Self::mul_reduce(u, &Self::R3)
        } else {
            Self::mul_reduce(u, &Self::R2)
        }
    }

    fn to_bytes(self) -> Vec<u8> {
        self.0.into()
    }

    fn random() -> Self {
        Self::from(BigInt::<NUM_LIMBS>::random())
    }

    fn to_string(self) -> String {
        self.rev_reduce().to_string()
    }
}

impl From<&BigInt<NUM_LIMBS>> for Fr<NUM_LIMBS> {
    fn from(value: &BigInt<NUM_LIMBS>) -> Self {
        Self::reduce(value, Some(false))
    }
}

impl From<BigInt<NUM_LIMBS>> for Fr<NUM_LIMBS> {
    fn from(value: BigInt<NUM_LIMBS>) -> Self {
        Self::reduce(&value, Some(false))
    }
}

impl Into<BigInt<NUM_LIMBS>> for Fr<NUM_LIMBS> {
    fn into(self) -> BigInt<NUM_LIMBS> {
        self.rev_reduce()
    }
}

// non-reduced transformation between words and field
impl From<[Word; NUM_LIMBS]> for Fr<NUM_LIMBS> {
    fn from(bytes: [Word; NUM_LIMBS]) -> Self {
        Fr(BigInt(bytes))
    }
}
impl Into<[Word; NUM_LIMBS]> for Fr<NUM_LIMBS> {
    fn into(self) -> [Word; NUM_LIMBS] {
        self.0 .0
    }
}

// non-reduced transformation between bytes and field
impl From<[Byte; BYTE_SIZE]> for Fr<NUM_LIMBS> {
    fn from(bytes: [Byte; BYTE_SIZE]) -> Self {
        Fr(BigInt::<NUM_LIMBS>::from(bytes.as_slice()))
    }
}
impl Into<[Byte; BYTE_SIZE]> for Fr<NUM_LIMBS> {
    fn into(self) -> [Byte; BYTE_SIZE] {
        let bytes: Vec<Byte> = self.0.into();
        bytes[..BYTE_SIZE].try_into().unwrap()
    }
}

impl<'a> Sum<&'a Self> for Fr<NUM_LIMBS> {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = &'a Self>,
    {
        iter.fold(Fr::<NUM_LIMBS>::ZERO(), |a, b| a + *b)
    }
}

impl Sum<Self> for Fr<NUM_LIMBS> {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        iter.fold(Fr::<NUM_LIMBS>::ZERO(), |a, b| a + b)
    }
}

impl Add for Fr<NUM_LIMBS> {
    type Output = Fr<NUM_LIMBS>;

    fn add(self, other: Self) -> Fr<NUM_LIMBS> {
        let (w, carrier) = self.0 + other.0;
        if carrier > 0 {
            Self(((BigInt::<NUM_LIMBS>::MAX() - Self::MODULUS).0 + w).0)
        } else if w >= Self::MODULUS {
            Self((w - Self::MODULUS).0)
        } else {
            Self(w)
        }
    }
}

impl Sub for Fr<NUM_LIMBS> {
    type Output = Fr<NUM_LIMBS>;

    fn sub(self, other: Self) -> Fr<NUM_LIMBS> {
        if self.0 > other.0 {
            Self((self.0 - other.0).0)
        } else {
            Self((Self::MODULUS - (other.0 - self.0).0).0)
        }
    }
}

impl Mul for Fr<NUM_LIMBS> {
    type Output = Fr<NUM_LIMBS>;

    fn mul(self, other: Self) -> Fr<NUM_LIMBS> {
        Self::mul_reduce(&self.0, &other.0)
    }
}

impl Mul<Word> for Fr<NUM_LIMBS> {
    type Output = Fr<NUM_LIMBS>;

    fn mul(self, other: Word) -> Fr<NUM_LIMBS> {
        Self::mul_reduce(
            &self.0,
            &Self::from_str(other.to_string().as_str()).unwrap().0,
        )
    }
}

impl Div for Fr<NUM_LIMBS> {
    type Output = Fr<NUM_LIMBS>;

    fn div(self, other: Self) -> Fr<NUM_LIMBS> {
        Self::mul_reduce(&self.0, &other.inverse().unwrap().0)
    }
}

impl Neg for Fr<NUM_LIMBS> {
    type Output = Fr<NUM_LIMBS>;

    fn neg(self) -> Fr<NUM_LIMBS> {
        Self::ZERO() - self
    }
}

impl FromStr for Fr<NUM_LIMBS> {
    type Err = FieldParseErr;

    fn from_str(text: &str) -> Result<Self, Self::Err> {
        let num = BigInt::<NUM_LIMBS>::from_str(text)
            .map_err(|_| FieldParseErr)
            .unwrap();
        Ok(Self::reduce(&num, Some(false)))
    }
}

impl Field<NUM_LIMBS> for Fr<NUM_LIMBS> {
    type BasePrimeField = Self;
    type BasePrimeFieldIter = iter::Once<Self::BasePrimeField>;

    fn extension_degree() -> u64 {
        1 as u64
    }

    fn to_base_prime_field_elements(&self) -> Self::BasePrimeFieldIter {
        iter::once(*self)
    }

    fn from_base_prime_field_elems(
        elems: impl IntoIterator<Item = Self::BasePrimeField>,
    ) -> Option<Self> {
        let mut elems = elems.into_iter();
        let elem = elems.next()?;
        if elems.next().is_some() {
            return None;
        }
        Some(elem)
    }

    fn from_base_prime_field_elem(elem: Self::BasePrimeField) -> Self {
        elem
    }

    fn legendre(&self) -> LegendreSymbol {
        let result = self.pow(((Self::MODULUS - BigInt::<NUM_LIMBS>::ONE()).0 >> 1).0);
        if result == Fr::<NUM_LIMBS>::from_str("1").unwrap() {
            LegendreSymbol::QuadraticResidue
        } else if result.is_zero() {
            LegendreSymbol::Zero
        } else {
            LegendreSymbol::QuadraticNonResidue
        }
    }

    // frobenius map of base prime field is itself
    fn powers_frobenius_map_inplace(&mut self, power: usize) {}

    #[inline(always)]
    fn ONE() -> Self {
        Self(Self::R)
    }

    #[inline(always)]
    fn ZERO() -> Self {
        Self(BigInt::<NUM_LIMBS>::ZERO())
    }

    // referenced from Algorithm 11.12 of "handbook of elliptic and hyperelliptic curve cryptography"
    // (aR)^{-1} % N <- ((aR)^{-1} * R^2) * R^{-1} % N
    fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }
        let (mut r, mut s, mut t, mut v) = (
            self.0,
            BigInt::<NUM_LIMBS>::ONE(),
            Self::MODULUS,
            BigInt::<NUM_LIMBS>::ZERO(),
        );

        // println!("r = {:?}", r);

        /////////////////////////////////////////////// STAGE ONE
        let mut k = 0 as usize;
        let mut tmp_carrier = 0 as Word;
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
        let m = NUM_LIMBS * WORD_SIZE as usize;
        // println!("------ k = {}, m = {}, v = {:?}", k, m, v);
        if k < m {
            (v, k) = (Self::mul_reduce(&v, &Self::R2).0, k + m);
        }

        // 2^{2m - k}
        let mut h_bit = 2 * m - k;
        assert!(h_bit < m);
        let x: Vec<Word> = (0..NUM_LIMBS)
            .map(|i| i * WORD_SIZE as usize)
            .map(|v| {
                // !!! 1/31/2024 fixed
                if (h_bit >= v) && (h_bit - v < WORD_SIZE) {
                    (1 as Word) << (h_bit - v)
                } else {
                    0 as Word
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

        Some(Self(v))
    }

    // general method of SQRT for finite field through Tonelli and Shanks Algorithm
    // referenced from Algorithm 11.23 of "handbook of elliptic and hyperelliptic curve cryptography"
    fn sqrt(&self) -> Option<Self> {
        // not quadratic residual
        if self.legendre() == LegendreSymbol::QuadraticNonResidue {
            return None;
        }

        let (n, r) = (Self::reduce(&Self::N, Some(false)), Self::RODD);
        let (a, z) = (self.clone(), n.pow(r));
        let (mut y, mut s, mut x, mut b) = (
            z,
            Self::E,
            a.pow(((r - BigInt::<NUM_LIMBS>::ONE()).0 >> 1).0),
            Self::ZERO(),
        );
        let tmp_x = a * x.clone();
        (b, x) = (tmp_x.clone() * x, tmp_x);
        while b != Self::ONE() {
            let mut m = 0 as Word;
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
        // // FOR DEBUG
        // if self.is_quadratic_residual() {
        //     assert!(x * x == *self);
        // }
        Some(x)
    }

    fn pow(&self, e: BigInt<NUM_LIMBS>) -> Self {
        // let n_bits: Vec<u8> = e.into();
        let n_bits: Vec<u8> = e.to_bits();
        let (mut y, x) = (Self::ONE().0, self.0);
        for i in (0..n_bits.len()).rev() {
            y = Self::mul_reduce(&y, &y).0;
            if n_bits[i] == 1 {
                y = Self::mul_reduce(&x, &y).0;
            }
        }
        Self(y)
    }

    fn is_zero(&self) -> bool {
        *self == Self::ZERO()
    }

    fn square(&self) -> Self {
        Self::mul_reduce(&self.0, &self.0)
    }

    fn square_inplace(&mut self) {
        *self = self.square();
    }
}
