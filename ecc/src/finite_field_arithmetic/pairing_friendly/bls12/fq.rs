/// Practical Implementation for Base Field (Fq) of BLS12-381
///
///
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
        Self::mul_reduce(&self.0, &other.inverse().unwrap().0)
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

impl Field<NUM_LIMBS> for Fq<NUM_LIMBS> {
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
        if result == Fq::<NUM_LIMBS>::from_str("1").unwrap() {
            LegendreSymbol::QuadraticResidue
        } else if result.is_zero() {
            LegendreSymbol::Zero
        } else {
            LegendreSymbol::QuadraticNonResidue
        }
    }

    // frobenius map of base prime field is itself
    fn powers_frobenius_map_inplace(&mut self, power: usize){}

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

mod tests {
    use super::*;
    use crate::finite_field_arithmetic::pairing_friendly::bls12::fr::Fr;

    #[test]
    fn test_frobenius_coeff() {
        let fp2_frob_coeff = [
            "1", "4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559786"
        ];
        let fp6_frob_coeff = [
            [
                ["1", "0"], 
                ["0", "4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939436"], 
                ["793479390729215512621379701633421447060886740281060493010456487427281649075476305620758731620350", "0"], 
                ["0", "1"], 
                ["4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939436", "0"], 
                ["0", "793479390729215512621379701633421447060886740281060493010456487427281649075476305620758731620350"]
            ], 
            [
                ["1", "0"], 
                ["4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939437", "0"], 
                ["4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939436", "0"], 
                ["4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559786", "0"], 
                ["793479390729215512621379701633421447060886740281060493010456487427281649075476305620758731620350", "0"], 
                ["793479390729215512621379701633421447060886740281060493010456487427281649075476305620758731620351", "0"]
            ]
        ];
        let fp12_frob_coeff = [
            ["1", "0"], 
            ["3850754370037169011952147076051364057158807420970682438676050522613628423219637725072182697113062777891589506424760", "151655185184498381465642749684540099398075398968325446656007613510403227271200139370504932015952886146304766135027"], 
            ["793479390729215512621379701633421447060886740281060493010456487427281649075476305620758731620351", "0"], 
            ["2973677408986561043442465346520108879172042883009249989176415018091420807192182638567116318576472649347015917690530", "1028732146235106349975324479215795277384839936929757896155643118032610843298655225875571310552543014690878354869257"], 
            ["793479390729215512621379701633421447060886740281060493010456487427281649075476305620758731620350", "0"], 
            ["3125332594171059424908108096204648978570118281977575435832422631601824034463382777937621250592425535493320683825557", "877076961050607968509681729531255177986764537961432449499635504522207616027455086505066378536590128544573588734230"], 
            ["4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559786", "0"], 
            ["151655185184498381465642749684540099398075398968325446656007613510403227271200139370504932015952886146304766135027", "3850754370037169011952147076051364057158807420970682438676050522613628423219637725072182697113062777891589506424760"], 
            ["4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939436", "0"], 
            ["1028732146235106349975324479215795277384839936929757896155643118032610843298655225875571310552543014690878354869257", "2973677408986561043442465346520108879172042883009249989176415018091420807192182638567116318576472649347015917690530"], 
            ["4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939437", "0"], 
            ["877076961050607968509681729531255177986764537961432449499635504522207616027455086505066378536590128544573588734230", "3125332594171059424908108096204648978570118281977575435832422631601824034463382777937621250592425535493320683825557"]
        ];

        println!("\n\n Coefficients of Frobenius Map over Fp2: ");
        let c0 = Fq::<NUM_LIMBS>::from_str(fp2_frob_coeff[0]).unwrap();
        let c1 = Fq::<NUM_LIMBS>::from_str(fp2_frob_coeff[1]).unwrap();
        println!("##{}: [{:?}, {:?}]\n", 1, c0, c1);

        println!("\n\n Coefficients of Frobenius Map over Fp6: ");
        println!("\nC1 part: ");
        for i in 0..6 {
            let c0 = Fq::<NUM_LIMBS>::from_str(fp6_frob_coeff[0][i][0]).unwrap();
            let c1 = Fq::<NUM_LIMBS>::from_str(fp6_frob_coeff[0][i][1]).unwrap();
            println!("##{}: [{:?}, {:?}]\n", i, c0, c1);
        }
        println!("\nC2 part: ");
        for i in 0..6 {
            let c0 = Fq::<NUM_LIMBS>::from_str(fp6_frob_coeff[1][i][0]).unwrap();
            let c1 = Fq::<NUM_LIMBS>::from_str(fp6_frob_coeff[1][i][1]).unwrap();
            println!("##{}: [{:?}, {:?}]\n", i, c0, c1);
        }

        println!("\n\n Coefficients of Frobenius Map over Fp12: ");
        for i in 0..12 {
            let c0 = Fq::<NUM_LIMBS>::from_str(fp12_frob_coeff[i][0]).unwrap();
            let c1 = Fq::<NUM_LIMBS>::from_str(fp12_frob_coeff[i][1]).unwrap();
            println!("##{}: [{:?}, {:?}]\n", i, c0, c1);
        }
    }

    #[test]
    fn test_g1_params() {
        let (A_str, B_str, g_x_str, g_y_str, cofactor_str, order_str) = (
            "0",
            "4",
            "2262513090815062280530798313005799329941626325687549893214867945091568948276660786250917700289878433394123885724147",
            "3165530325623507257754644679249908411459467330345960501615736676710739703656949057125324800107717061311272030899084",
            "76329603384216526031706109802092473003",
            "4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129030796414117214202539",
        );
        let (A, B, g_x, g_y, cofactor, order) = (
            Fq::<NUM_LIMBS>::from_str(A_str).unwrap(),
            Fq::<NUM_LIMBS>::from_str(B_str).unwrap(),
            Fq::<NUM_LIMBS>::from_str(g_x_str).unwrap(),
            Fq::<NUM_LIMBS>::from_str(g_y_str).unwrap(),
            BigInt::<8>::from_str(cofactor_str).unwrap(),
            BigInt::<NUM_LIMBS>::from_str(order_str).unwrap(),
        );
        println!("## Parameters for G1 group: \n ");
        println!("A = {:?} \n B = {:?} \n\n", A, B);
        println!("g = ({:?}, {:?}) \n\n", g_x, g_y);
        println!("cofactor = {:?} \n order = {:?}", cofactor, order);
    }

    #[test]
    fn test_into_bigint() {
        // let a = "491024377652271573700152936577795338400293278102273691478251583805812369999324468448692790549212082767413034416424";
        let a = "4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559786";
        let a_bigint: BigInt<NUM_LIMBS> = Fq::<NUM_LIMBS>::from_str(a).unwrap().into();
        println!("a = {:?}", Fq::<NUM_LIMBS>::from_str(a).unwrap());
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
            Fq::<NUM_LIMBS>::from_str(a).unwrap().sqrt().unwrap(),
            Fq::<NUM_LIMBS>::from_str(c).unwrap()
        );
    }
}
