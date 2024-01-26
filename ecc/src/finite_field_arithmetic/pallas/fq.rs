/// Practical Implementation for Base Field (Fq) of Pallas Curve
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
pub struct Fq<const N: usize>(pub BigInt<N>);

const NUM_LIMBS: usize = 4;
const WORD_SIZE: usize = 64;
type Word = u64;

impl PrimeField<NUM_LIMBS> for Fq<NUM_LIMBS> {
    // fabricated precomputable parameters of custom finite field,
    // these constant parameters need to determined at compile time

    // MODULUS = 28948022309329048855892746252171976963363056481941560715954676764349967630337
    // W = 2^64
    // r = W^4 = 115792089237316195423570985008687907853269984665640564039457584007913129639936
    // R = r % MODULUS = 28948022309329048855892746252171976963180815219815881891593553714863226748925
    // R2 = r * r % MODULUS = 4263855311831330276397237192126260515652039413828781833859739249380679483407
    // R3 = r * r * r % MODULUS = 19398276961315000371481654775825491914897503251658641052577562551434289746169
    // M0 = (-MODULUS^{-1} % r) % W = 7409212017489215489
    //
    // let: MODULUS = 2^E * ROOD + 1
    // E = 32
    // RODD = 6739986666787659948666753771754907668419893943225396963757154709741
    // N = 5/7/...., is sampled number whose legendre symbol is -1 (is definitely non-quadratic residual)
    const MODULUS: BigInt<NUM_LIMBS> = BigInt([
        11037532056220336129,
        2469829653914515739,
        0,
        4611686018427387904,
    ]);
    const R: BigInt<NUM_LIMBS> = BigInt([
        3780891978758094845,
        11037255111966004397,
        18446744073709551615,
        4611686018427387903,
    ]);
    const R2: BigInt<NUM_LIMBS> = BigInt([
        10122100416058490895,
        15551789045973377255,
        8617542898466512152,
        679271340751763220,
    ]);
    const R3: BigInt<NUM_LIMBS> = BigInt([
        17403498412575166713,
        17773050464821424593,
        16108549121152898092,
        3090323811697793296,
    ]);
    const M0: Word = 11037532056220336127 as Word;
    const E: Word = 32 as Word;
    const RODD: BigInt<NUM_LIMBS> = BigInt([670184341500670189, 575052028, 0, 1073741824]);
    const N: BigInt<NUM_LIMBS> = BigInt([5, 0, 0, 0]);

    #[inline(always)]
    fn ONE() -> Self {
        Self(Self::R)
    }

    #[inline(always)]
    fn ZERO() -> Self {
        Self(BigInt::<NUM_LIMBS>::ZERO())
    }

    // referenced from Algorithm 11.12 of "handbook of elliptic and hyperelliptic curve cryptography"
    // (a^{-1} * R) % N <- (a^{-1} * R^2) * R^{-1} % N
    fn inv(&self) -> Self {
        let (mut r, mut s, mut t, mut v) = (
            self.0,
            BigInt::<NUM_LIMBS>::ONE(),
            Self::MODULUS,
            BigInt::<NUM_LIMBS>::ZERO(),
        );

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
        let h_bit = 2 * m - k;
        assert!(h_bit < m);
        let x: Vec<Word> = (0..NUM_LIMBS)
            .map(|i| i * WORD_SIZE as usize)
            .map(|v| {
                if (h_bit > v) && (h_bit - v < WORD_SIZE) {
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

    // Tonelli and Shanks Algorithm adapted for a special modulus p, satisfying: p MOD 8 == 1
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
        self.pow(((Self::MODULUS - BigInt::<NUM_LIMBS>::ONE()).0 >> 1).0) == Self::ONE()
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
        let s = 4;
        let mut t = BigInt([0 as Word; 4]);
        let (mut c1, mut c2, mut overflow) = (0 as Word, 0 as Word, false);
        // println!("lft = {:?}, rht = {:?}", lft.0, rht.0);
        for i in 0..s {
            // t = t + self * other[i]
            let (mut tmp_c1, mut tmp_c2, mut ab) = (0 as Word, 0 as Word, BigInt([0 as Word; 4]));
            (ab, tmp_c1) = lft.clone() * rht.0[i];
            (t, tmp_c2) = t + ab;
            (c1, overflow) = c1.overflowing_add(tmp_c1.wrapping_add(tmp_c2));
            if overflow {
                c2 += 1 as Word;
            }
            // println!("i = {}, #1: t = {:?}, c1 = {:?}, c2 = {:?}", i, t.0, c1, c2);

            // t = t + ((t[0] * N'[0]) mod W) * N
            let (mut tmp_c3, mut tmp_c4, mut mn) = (0 as Word, 0 as Word, BigInt([0 as Word; 4]));
            let m = Self::M0.wrapping_mul(t.0[0]);
            (mn, tmp_c3) = Self::MODULUS * m;
            (t, tmp_c4) = t + mn;
            (c1, overflow) = c1.overflowing_add(tmp_c3.wrapping_add(tmp_c4));
            if overflow {
                c2 += 1 as Word;
            }
            // println!("i = {} #2: t = {:?}, c1 = {:?}, c2 = {:?}", i, t.0, c1, c2);

            // t >> 1
            for j in 0..(s - 1) {
                t.0[j] = t.0[j + 1];
            }
            (t.0[s - 1], c1, c2) = (c1, c2, 0 as Word);
            // println!("i = {} #3: t = {:?}, c1 = {:?}, c2 = {:?}", i, t.0, c1, c2);
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

// non-reduced transformation between bytes and field
impl From<[Word; NUM_LIMBS]> for Fq<NUM_LIMBS> {
    fn from(bytes: [Word; NUM_LIMBS]) -> Self {
        Fq(BigInt(bytes))
    }
}

impl Into<[Word; 4]> for Fq<NUM_LIMBS> {
    fn into(self) -> [Word; NUM_LIMBS] {
        self.0 .0
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() {
        let a = "23299583932281281605549117567410625401281091317674300085390516703976038929236";
        let a_bigint: BigInt<NUM_LIMBS> = Fq::<NUM_LIMBS>::from_str(a).unwrap().into();
        assert_eq!(a_bigint, BigInt::<NUM_LIMBS>::from_str(a).unwrap());

        let c = "23165056470023441196371639597826156042421324766918472161609664258633197933423";
        let a_field: [Word; NUM_LIMBS] = Fq::<NUM_LIMBS>::from_str(a).unwrap().into();
        assert_eq!(a_field, BigInt::<NUM_LIMBS>::from_str(c).unwrap().0);
    }

    #[test]
    fn test_pallas_curve_params() {
        // params of weierstrass model: a = 0, b = 5
        let (a, b) = ("0", "5");
        println!(
            "a = {:?}, b = {:?}",
            Fq::<NUM_LIMBS>::from_str(a).unwrap().0 .0,
            Fq::<NUM_LIMBS>::from_str(b).unwrap().0 .0
        );

        // generator: (-1, 2)
        let (gen_x, gen_y) = (
            "28948022309329048855892746252171976963363056481941560715954676764349967630336",
            "2",
        );
        println!(
            "generator: x = {:?}, y = {:?}",
            Fq::<NUM_LIMBS>::from_str(gen_x).unwrap().0 .0,
            Fq::<NUM_LIMBS>::from_str(gen_y).unwrap().0 .0,
        );

        // order of curve(group): 28948022309329048855892746252171976963363056481941647379679742748393362948097
        let order = "28948022309329048855892746252171976963363056481941647379679742748393362948097";
        println!(
            "order: {:?}",
            BigInt::<NUM_LIMBS>::from_str(order).unwrap().0
        );
    }

    #[test]
    fn test_addition() {
        let (a, b, c) = (
            "23299583932281281605549117567410625401281091317674300085390516703976038929236",
            "12728795896941878947827723216087007012404899953936068154162936317783913641026",
            "7080357519894111697484094531325655450322934789668807523598776257409984939925",
        );
        assert_eq!(
            Fq::<NUM_LIMBS>::from_str(a).unwrap() + Fq::<NUM_LIMBS>::from_str(b).unwrap(),
            Fq::<NUM_LIMBS>::from_str(c).unwrap()
        );
    }

    #[test]
    fn test_substraction() {
        let (a, b, c) = (
            "28244056179519182624834465195545514411138275528379130540542625181634794962600",
            "13583400149048208564886004143663881202478291013837410260892386282926681992841",
            "14660656030470974059948461051881633208659984514541720279650238898708112969759",
        );
        assert_eq!(
            Fq::<NUM_LIMBS>::from_str(a).unwrap() - Fq::<NUM_LIMBS>::from_str(b).unwrap(),
            Fq::<NUM_LIMBS>::from_str(c).unwrap()
        );
    }

    #[test]
    fn test_multiplication() {
        let (a, b, c) = (
            "28420585563447636701388901488730408229744882268550871135427378007593007145818",
            "10929005576730665445090409463623600211518111311797436516349444987664000120035",
            "20884984882123841859966737092202213092058578554401209905127510218020991566255",
        );
        assert_eq!(
            Fq::<NUM_LIMBS>::from_str(a).unwrap() * Fq::<NUM_LIMBS>::from_str(b).unwrap(),
            Fq::<NUM_LIMBS>::from_str(c).unwrap()
        );
    }

    #[test]
    fn test_inversion() {
        let (a, c) = (
            "25767596886874889036540765742642627840896260537012642841805800608653635566576",
            "12977814693209616701908158718117218571334991308499348109965704251722539193017",
        );
        assert_eq!(
            Fq::<NUM_LIMBS>::from_str(a).unwrap() * Fq::<NUM_LIMBS>::from_str(c).unwrap(),
            Fq::<NUM_LIMBS>::ONE()
        );
        assert_eq!(
            Fq::<NUM_LIMBS>::from_str(a).unwrap().inv(),
            Fq::<NUM_LIMBS>::from_str(c).unwrap()
        );
    }

    #[test]
    fn test_division() {
        let (a, b, c) = (
            "11094102541548518406083044310753010621275594174978757610119501797021534623670",
            "25767596886874889036540765742642627840896260537012642841805800608653635566576",
            "10184968062787876042028836773955231702507653654365225893745657195078661715473",
        );
        assert_eq!(
            Fq::<NUM_LIMBS>::from_str(a).unwrap() / Fq::<NUM_LIMBS>::from_str(b).unwrap(),
            Fq::<NUM_LIMBS>::from_str(c).unwrap()
        );
    }

    #[test]
    fn test_sqrt() {
        let (a, c) = (
            "761940212266856713371586569342150604283558917968569240208532761798026301469",
            "4602589297423635878136638089713975619925356628567217110388194878773233887829",
        );
        let lft = Fq::<NUM_LIMBS>::from_str(a).unwrap();
        let result_1 = Fq::<NUM_LIMBS>::from_str(c).unwrap();

        assert_eq!(lft.is_quadratic_residual(), true);
        assert_eq!(
            result_1.pow(BigInt::<NUM_LIMBS>::from_str("2").unwrap()),
            lft
        );

        let result_2 = lft.sqrt();
        assert_eq!(result_2 * result_2, lft);
        println!(
            "two solutions for square root of {:?}: result_1 = {:?}, result_2 = {:?}",
            lft.0, result_1.0, result_2.0
        );
    }

    #[test]
    fn test_pow() {
        let (a, b, c) = (
            "22872133923507589581064940973199724574277684002088185742903677959102317569577",
            "234",
            "27818977364773149908565484571278600710302550015293433604567758323829469042334",
        );
        assert_eq!(
            Fq::<NUM_LIMBS>::from_str(a)
                .unwrap()
                .pow(BigInt::<NUM_LIMBS>::from_str(b).unwrap()),
            Fq::<NUM_LIMBS>::from_str(c).unwrap()
        );
    }
}
