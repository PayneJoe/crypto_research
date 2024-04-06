/// Practical Implementation for Base Field (Fq) of BLS12-381
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

#[derive(Debug, PartialEq, Eq, Clone, Copy, Hash)]
pub struct Fq<const N: usize>(pub BigInt<N>);

// since the bit length of Fq is 381, so we need 6 words to hold
const NUM_LIMBS: usize = 6;
const WORD_SIZE: usize = 64;
const BYTE_SIZE: usize = NUM_LIMBS * WORD_SIZE;
type Word = u64;
type Byte = u8;

// MODULUS = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
// W = 2^64 = 18446744073709551616
// r = W^6 = 2^384 = 39402006196394479212279040100143613805079739270465446667948293404245721771497210611414266254884915640806627990306816
// R = r % MODULUS = 3380320199399472671518931668520476396067793891014375699959770179129436917079669831430077592723774664465579537268733
// R2 = r * r % MODULUS = 2708263910654730174793787626328176511836455197166317677006154293982164122222515399004018013397331347120527951271750
// R3 = r * r * r % MODULUS = 1639067542774625894236716575548084905938753837211594095883637014582201460755008380976950835174037649440777609978336
//
// M0 = (-MODULUS^{-1} % r) % W = 18446744073709551616 - 8506173809081122819 = 9940570264628428797
//
// let: MODULUS = 2^E * ROOD + 1
// E = 1
// RODD = 2001204777610833696708894912867952078278441409969503942666029068062015825245418932221343814564507832018947136279893
// N = 5/7/...., is sampled number whose legendre symbol is -1 (is definitely non-quadratic residual)

impl PrimeField<NUM_LIMBS> for Fq<NUM_LIMBS> {
    const MODULUS: BigInt<NUM_LIMBS> = BigInt([
        13402431016077863595,
        2210141511517208575,
        7435674573564081700,
        7239337960414712511,
        5412103778470702295,
        1873798617647539866,
    ]);
    const R: BigInt<NUM_LIMBS> = BigInt([
        8505329371266088957,
        17002214543764226050,
        6865905132761471162,
        8632934651105793861,
        6631298214892334189,
        1582556514881692819,
    ]);
    const R2: BigInt<NUM_LIMBS> = BigInt([
        17644856173732828998,
        754043588434789617,
        10224657059481499349,
        7488229067341005760,
        11130996698012816685,
        1267921511277847466,
    ]);
    const R3: BigInt<NUM_LIMBS> = BigInt([
        17098105564519244256,
        3557706395579559416,
        11120290361046346205,
        3801124253586036577,
        2671430854784468776,
        767358375875140941,
    ]);
    const M0: Word = 9940570264628428797 as Word;
    const E: Word = 1 as Word;
    const RODD: BigInt<NUM_LIMBS> = BigInt([
        15924587544893707605,
        1105070755758604287,
        12941209323636816658,
        12843041017062132063,
        2706051889235351147,
        936899308823769933,
    ]);
    const N: BigInt<NUM_LIMBS> = BigInt([5, 0, 0, 0, 0, 0]);
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
    fn inv(&self) -> Self {
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

        Self(v)
    }

    // general method of SQRT for finite field through Tonelli and Shanks Algorithm
    // referenced from Algorithm 11.23 of "handbook of elliptic and hyperelliptic curve cryptography"
    fn sqrt(&self) -> Self {
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
        // FOR DEBUG
        if self.is_quadratic_residual() {
            assert!(x * x == *self);
        }
        x
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

    fn is_quadratic_residual(self) -> bool {
        self.pow(((Self::MODULUS - BigInt::<NUM_LIMBS>::ONE()).0 >> 1).0) == Self::ONE()
    }

    fn is_zero(self) -> bool {
        self == Self::ZERO()
    }

    fn square(&self) -> Self {
        Self::mul_reduce(&self.0, &self.0)
    }

    fn square_inplace(&mut self) {
        *self = Self::mul_reduce(&(*self).0, &(*self).0);
    }

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

impl From<&BigInt<NUM_LIMBS>> for Fq<NUM_LIMBS> {
    fn from(value: &BigInt<NUM_LIMBS>) -> Self {
        Self::reduce(value, Some(false))
    }
}

impl From<BigInt<NUM_LIMBS>> for Fq<NUM_LIMBS> {
    fn from(value: BigInt<NUM_LIMBS>) -> Self {
        Self::reduce(&value, Some(false))
    }
}

impl Into<BigInt<NUM_LIMBS>> for Fq<NUM_LIMBS> {
    fn into(self) -> BigInt<NUM_LIMBS> {
        self.rev_reduce()
    }
}

// non-reduced transformation between words and field
impl From<[Word; NUM_LIMBS]> for Fq<NUM_LIMBS> {
    fn from(bytes: [Word; NUM_LIMBS]) -> Self {
        Fq(BigInt(bytes))
    }
}
impl Into<[Word; NUM_LIMBS]> for Fq<NUM_LIMBS> {
    fn into(self) -> [Word; NUM_LIMBS] {
        self.0 .0
    }
}

// non-reduced transformation between bytes and field
impl From<[Byte; BYTE_SIZE]> for Fq<NUM_LIMBS> {
    fn from(bytes: [Byte; BYTE_SIZE]) -> Self {
        Fq(BigInt::<NUM_LIMBS>::from(bytes.as_slice()))
    }
}
impl Into<[Byte; BYTE_SIZE]> for Fq<NUM_LIMBS> {
    fn into(self) -> [Byte; BYTE_SIZE] {
        let bytes: Vec<Byte> = self.0.into();
        bytes[..BYTE_SIZE].try_into().unwrap()
    }
}

impl<'a> Sum<&'a Self> for Fq<NUM_LIMBS> {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = &'a Self>,
    {
        iter.fold(Fq::<NUM_LIMBS>::ZERO(), |a, b| a + *b)
    }
}

impl Sum<Self> for Fq<NUM_LIMBS> {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        iter.fold(Fq::<NUM_LIMBS>::ZERO(), |a, b| a + b)
    }
}

impl Add for Fq<NUM_LIMBS> {
    type Output = Fq<NUM_LIMBS>;

    fn add(self, other: Self) -> Fq<NUM_LIMBS> {
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

impl Sub for Fq<NUM_LIMBS> {
    type Output = Fq<NUM_LIMBS>;

    fn sub(self, other: Self) -> Fq<NUM_LIMBS> {
        if self.0 > other.0 {
            Self((self.0 - other.0).0)
        } else {
            Self((Self::MODULUS - (other.0 - self.0).0).0)
        }
    }
}

impl Mul for Fq<NUM_LIMBS> {
    type Output = Fq<NUM_LIMBS>;

    fn mul(self, other: Self) -> Fq<NUM_LIMBS> {
        Self::mul_reduce(&self.0, &other.0)
    }
}

impl Mul<Word> for Fq<NUM_LIMBS> {
    type Output = Fq<NUM_LIMBS>;

    fn mul(self, other: Word) -> Fq<NUM_LIMBS> {
        Self::mul_reduce(
            &self.0,
            &Self::from_str(other.to_string().as_str()).unwrap().0,
        )
    }
}

impl Div for Fq<NUM_LIMBS> {
    type Output = Fq<NUM_LIMBS>;

    fn div(self, other: Self) -> Fq<NUM_LIMBS> {
        Self::mul_reduce(&self.0, &other.inv().0)
    }
}

impl Neg for Fq<NUM_LIMBS> {
    type Output = Fq<NUM_LIMBS>;

    fn neg(self) -> Fq<NUM_LIMBS> {
        Self::ZERO() - self
    }
}

impl FromStr for Fq<NUM_LIMBS> {
    type Err = FieldParseErr;

    fn from_str(text: &str) -> Result<Self, Self::Err> {
        let num = BigInt::<NUM_LIMBS>::from_str(text)
            .map_err(|_| FieldParseErr)
            .unwrap();
        Ok(Self::reduce(&num, Some(false)))
    }
}

mod tests {
    use super::*;

    #[test]
    fn test_into_bigint() {
        let a = "491024377652271573700152936577795338400293278102273691478251583805812369999324468448692790549212082767413034416424";
        let a_bigint: BigInt<NUM_LIMBS> = Fq::<NUM_LIMBS>::from_str(a).unwrap().into();
        assert_eq!(a_bigint, BigInt::<NUM_LIMBS>::from_str(a).unwrap());
    }

    #[test]
    fn test_into_field() {
        let a = "491024377652271573700152936577795338400293278102273691478251583805812369999324468448692790549212082767413034416424";
        let c = "1861486499784848706405258076333327478705098756586036708891566981207117285764111680517329729099803084735393228152325";
        let a_field: [Word; NUM_LIMBS] = Fq::<NUM_LIMBS>::from_str(a).unwrap().into();
        assert_eq!(a_field, BigInt::<NUM_LIMBS>::from_str(c).unwrap().0);
    }

    #[test]
    fn test_addition() {
        let (a, b, c) = (
            "1936685924074983148182060051590223173080719226074312827735283806945363245574769831424458169765159639840258580980296",
            "1499194227591321698461389819429023209536270652192322866826113868699015445822010524616316855955332437260591959261827",
            "3435880151666304846643449871019246382616989878266635694561397675644378691396780356040775025720492077100850540242123",
        );
        assert_eq!(
            Fq::<NUM_LIMBS>::from_str(a).unwrap() + Fq::<NUM_LIMBS>::from_str(b).unwrap(),
            Fq::<NUM_LIMBS>::from_str(c).unwrap()
        );
    }

    #[test]
    fn test_substraction() {
        let (a, b, c) = (
            "1936685924074983148182060051590223173080719226074312827735283806945363245574769831424458169765159639840258580980296",
            "1499194227591321698461389819429023209536270652192322866826113868699015445822010524616316855955332437260591959261827",
            "437491696483661449720670232161199963544448573881989960909169938246347799752759306808141313809827202579666621718469",
        );
        assert_eq!(
            Fq::<NUM_LIMBS>::from_str(a).unwrap() - Fq::<NUM_LIMBS>::from_str(b).unwrap(),
            Fq::<NUM_LIMBS>::from_str(c).unwrap()
        );
    }

    #[test]
    fn test_multiplication() {
        let (a, b, c) = (
            "1936685924074983148182060051590223173080719226074312827735283806945363245574769831424458169765159639840258580980296",
            "1499194227591321698461389819429023209536270652192322866826113868699015445822010524616316855955332437260591959261827",
            "807006275751014011246110611228296902519682995735942946709205252428279255385218161465128398467483445174243869941798",
        );
        assert_eq!(
            Fq::<NUM_LIMBS>::from_str(a).unwrap() * Fq::<NUM_LIMBS>::from_str(b).unwrap(),
            Fq::<NUM_LIMBS>::from_str(c).unwrap()
        );
    }

    #[test]
    fn test_division() {
        let (a, b, c) = (
            "1936685924074983148182060051590223173080719226074312827735283806945363245574769831424458169765159639840258580980296",
            "1499194227591321698461389819429023209536270652192322866826113868699015445822010524616316855955332437260591959261827",
            "1336785341618850264772783872266072954093474702642676025692416801937723794293608117542281462319701087670997565163102",
        );
        assert_eq!(
            Fq::<NUM_LIMBS>::from_str(a).unwrap() / Fq::<NUM_LIMBS>::from_str(b).unwrap(),
            Fq::<NUM_LIMBS>::from_str(c).unwrap()
        );
    }

    #[test]
    fn test_sqrt() {
        let (a, c) = (
            "1936685924074983148182060051590223173080719226074312827735283806945363245574769831424458169765159639840258580980296",
            "916407379142843325014424957531723593182479142585202804316555682277426719956241116057215429405117420784032484668809"
        );
        assert_eq!(
            Fq::<NUM_LIMBS>::from_str(a).unwrap().sqrt(),
            Fq::<NUM_LIMBS>::from_str(c).unwrap()
        );
    }
}