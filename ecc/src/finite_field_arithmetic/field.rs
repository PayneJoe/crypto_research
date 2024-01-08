use std::ops::{Add, Div, Mul, Sub};
use std::str::FromStr;

//////////////////////////////////// Implementation of BigInteger Specially for Finite Field
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct BI<const N: usize>(pub [u8; N]);

impl BI<2> {
    #[inline(always)]
    fn zero() -> Self {
        BI([0, 0])
    }

    #[inline(always)]
    fn one() -> Self {
        BI([1, 0])
    }
}

impl FromStr for BI<2> {
    type Err = ParseStrErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        unimplemented!()
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
            let t = lft_w - (rht_w + carrier_w);
            (carrier, remainder) = ((t >> 8) as u8, (t & ((1 << 8) - 1)) as u8);
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

// field configurations including basic behavious on fields
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

    // one
    fn one() -> Self;

    // reduce a bigint into the specified range [0, modulus)
    fn reduce(u: &BI<N>, inv: Option<bool>) -> Self;
    // convert a reduced number into a unreduced bigint
    fn rev_reduce(&self) -> BI<N>;
    // inversion of a reduced number
    fn inv(&self) -> Self;

    // Montgomery Reduction
    fn mul_reduce(lft: &BI<N>, rht: &BI<N>) -> Self;
}

//////////////////////////////////////// custom finite field
#[derive(Debug)]
pub struct ParseStrErr;

#[derive(Debug, PartialEq)]
pub struct Foo<const T: usize>(BI<T>);

impl FromStr for Foo<2> {
    type Err = ParseStrErr;

    fn from_str(text: &str) -> Result<Self, Self::Err> {
        // parse the large number
        let mut large_number = text.parse::<u32>().map_err(|_| ParseStrErr).unwrap();

        let apply_parse = |word: u16| -> BI<2> {
            // retain the lowest N words, ignore the high words
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
        // println!("\n******* parse {}", text);
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
        let (w, _) = self.0 + other.0;
        Self(w)
    }
}

impl Sub for Foo<2> {
    type Output = Foo<2>;

    fn sub(self, other: Self) -> Foo<2> {
        let (w, _) = self.0 - other.0;
        Self(w)
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
        unimplemented!()
    }
}

// define custom finite field
impl Field<2> for Foo<2> {
    // fabricated precomputable parameters of custom finite field
    // these parameters need to determined at compile time
    // W = 256, MODULUS/M = 517, R = 256^2 % M = 394, M0 = -M[0]^{-1} % W = 51
    const MODULUS: BI<2> = BI([5, 2]);
    const R: BI<2> = BI([138, 1]);
    const R2: BI<2> = BI([136, 0]);
    const R3: BI<2> = BI([77, 1]);
    const M0: u8 = 51_u8;

    fn reduce(u: &BI<2>, inv: Option<bool>) -> Self {
        if let Some(true) = inv {
            Self::mul_reduce(u, &Self::R3)
        } else {
            Self::mul_reduce(u, &Self::R2)
        }
    }

    fn one() -> Self {
        Self(BI([1_u8, 0_u8]))
    }

    // (aR)^{-1} % N
    fn inv(&self) -> Self {
        // let (mut A, mut B) = (Self::M0, self.0);
        // let (mut Ua, mut Ub) = (BI([0, 0]), BI([1, 0]));
        // let (mut Va, mut Vb) = (BI([1, 0]), BI([0, 0]));
        // let mut n_iter = 0;
        unimplemented!()
    }

    // a % N <- aR * R^{-1} % N
    fn rev_reduce(&self) -> BI<2> {
        Self::mul_reduce(&self.0, &Self::one().0).0
    }

    // abR % N <- (aR * bR) * R^{-1} % N
    fn mul_reduce(lft: &BI<2>, rht: &BI<2>) -> Self {
        // println!("------ lft = {:?}, rht = {:?}", lft, rht);

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
            // println!(
            //     "stage 1: i = {}, t = {},{},{:?}, ({:?}, {}) = {:?} * {}",
            //     i, c2, c1, t, ab, tmp_c1, lft, rht.0[i]
            // );

            // t = t + ((t[0] * N'[0]) mod W) * N
            let (mut tmp_c3, mut tmp_c4, mut mn) = (0_u8, 0_u8, BI([0_u8; 2]));
            let m = Self::M0.wrapping_mul(t.0[0]);
            (mn, tmp_c3) = Self::MODULUS * m;
            (t, tmp_c4) = t + mn;
            (c1, overflow) = c1.overflowing_add(tmp_c3.wrapping_add(tmp_c4));
            if overflow {
                c2 += 1_u8;
            }
            // println!(
            //     "stage 2: i = {}, t = {},{},{:?}, {:?} * {}",
            //     i,
            //     c2,
            //     c1,
            //     t,
            //     Self::MODULUS,
            //     m
            // );

            // t >> 1
            for j in 0..(s - 1) {
                t.0[j] = t.0[j + 1];
            }
            (t.0[s - 1], c1, c2) = (c1, c2, 0_u8);
            // println!("stage 3: i = {}, t = {},{},{:?}", i, c2, c1, t,);
        }

        Self(t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_conversion() {
        let a = 259_u16;
        let result: BI<2> = Foo::<2>::from_str(a.to_string().as_str()).unwrap().into();
        let actual = BI([(a % 256) as u8, (a / 256) as u8]);
        assert_eq!(result, actual);
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
}
