use std::cmp::Ordering;
use std::ops::{Add, Mul, Shl, Shr, Sub};
use std::str::FromStr;

#[derive(Debug)]
pub struct BigIntParseErr;

pub trait BigInteger<const N: usize>:
    FromStr + PartialOrd + Into<Vec<u8>> + Ord + Add + Mul<u64> + Shl<usize> + Shr<usize> + Sub
{
}
impl<const N: usize> BigInteger<N> for BigInt<N> {}

//////////////////////////////////// Practical Implementation of BigInteger Specially for Finite Field
#[derive(Debug, PartialEq, Clone, Copy, Eq)]
pub struct BigInt<const N: usize>(pub [u64; N]);

impl<const N: usize> BigInt<N> {
    #[inline(always)]
    pub fn MAX() -> Self {
        Self(vec![((1 << 64) - 1) as u64; N].try_into().unwrap())
    }

    #[inline(always)]
    pub fn ZERO() -> Self {
        BigInt([0 as u64; N])
    }

    #[inline(always)]
    pub fn ONE() -> Self {
        let mut data = vec![0 as u64];
        data.extend(vec![0 as u64; N - 1]);
        BigInt(data.try_into().unwrap())
    }

    #[inline(always)]
    pub fn is_even(self) -> bool {
        self.0[0] % 2 == 0
    }

    #[inline(always)]
    pub fn is_zero(self) -> bool {
        let n = self.0.len();
        let mut zero = true;
        for i in 0..n {
            zero = zero & (self.0[i] == 0 as u64);
        }
        zero
    }
}

impl<const N: usize> PartialOrd for BigInt<N> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<const N: usize> Ord for BigInt<N> {
    fn cmp(&self, other: &Self) -> Ordering {
        let n = N;
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

impl<const N: usize> Shl<usize> for BigInt<N> {
    type Output = (BigInt<N>, u64);
    fn shl(self, other: usize) -> Self::Output {
        let n = 64 as usize;
        assert!(other < n);
        let mut w = Self([0 as u64; N]);
        let mut carrier = 0 as u64;
        // from the lowest word to the highest word
        for i in 0..n {
            let tmp_carrier = self.0[i] >> (n - other);
            let remainder = (self.0[i] << other) | carrier;
            w.0[i] = remainder;
            carrier = tmp_carrier;
        }
        (w, carrier)
    }
}

impl<const N: usize> Shr<usize> for BigInt<N> {
    type Output = (BigInt<N>, u64);
    fn shr(self, other: usize) -> Self::Output {
        let n = 64 as usize;
        assert!(other < n);
        let mut w = Self([0 as u64; N]);
        let mut carrier = 0 as u64;
        // from the lowest word to the highest word
        for i in (0..n).rev() {
            let tmp_carrier = self.0[i] & (((1 << other) - 1) as u64);
            let remainder = (self.0[i] >> other) | (carrier << (n - other));
            w.0[i] = remainder;
            carrier = tmp_carrier;
        }
        (w, carrier)
    }
}

impl<const N: usize> Add for BigInt<N> {
    type Output = (BigInt<N>, u64);
    fn add(self, other: BigInt<N>) -> Self::Output {
        let n = N;
        let lft = &self.0;
        let rht = &other.0;
        let (mut carrier, mut remainder) = (0 as u64, 0 as u64);
        let mut w = BigInt([0 as u64; N]);
        for i in 0..n {
            let (lft_w, rht_w, carrier_w) =
                (u128::from(lft[i]), u128::from(rht[i]), u128::from(carrier));
            let t = lft_w + (rht_w + carrier_w);
            (carrier, remainder) = ((t >> n) as u64, (t & ((1 << n) - 1)) as u64);
            w.0[i] = remainder;
        }
        (w, carrier)
    }
}

impl<const N: usize> Sub for BigInt<N> {
    type Output = (BigInt<N>, u64);
    fn sub(self, other: BigInt<N>) -> Self::Output {
        let n = N;
        let lft = &self.0;
        let rht = &other.0;
        let (mut carrier, mut remainder) = (0 as u64, 0 as u64);
        let mut w = BigInt([0 as u64; N]);
        for i in 0..n {
            let (lft_w, rht_w, carrier_w) =
                (u128::from(lft[i]), u128::from(rht[i]), u128::from(carrier));
            (carrier, remainder) = if lft_w >= rht_w + carrier_w {
                (0_u64, (lft_w - (rht_w + carrier_w)) as u64)
            } else {
                (1_u64, (lft_w + (1 << n) - (rht_w + carrier_w)) as u64)
            };
            w.0[i] = remainder;
        }
        (w, carrier)
    }
}

impl<const N: usize> Mul<u64> for BigInt<N> {
    type Output = (BigInt<N>, u64);

    fn mul(self, other: u64) -> Self::Output {
        let n = N;
        let lft = &self.0;
        let (mut carrier, mut remainder) = (0 as u64, 0 as u64);
        let mut w = BigInt([0 as u64; N]);
        for i in 0..n {
            let (lft_w, rht_w, carrier_w) =
                (u128::from(lft[i]), u128::from(other), u128::from(carrier));
            let t = lft_w * rht_w + carrier_w;
            (carrier, remainder) = ((t >> n) as u64, (t & ((1 << n) - 1)) as u64);
            w.0[i] = remainder;
        }
        (w, carrier)
    }
}

impl<const N: usize> Into<Vec<u8>> for BigInt<N> {
    fn into(self) -> Vec<u8> {
        let n = N as u32;
        let num_leading_zeros = self.0.iter().rev().fold(0 as u32, |acc, v| {
            if acc == 0 {
                acc + v.leading_zeros()
            } else if (acc % n) == 0 {
                acc + v.leading_zeros()
            } else {
                acc
            }
        });
        let bits = (0..self.0.len())
            .map(|i| ((0..n).map(|j| ((self.0[i] >> j) & 1) as u8)).collect::<Vec<u8>>())
            .collect::<Vec<Vec<u8>>>()
            .into_iter()
            .flatten()
            .collect::<Vec<u8>>();
        let end = bits.len() - (num_leading_zeros as usize);
        bits[..end].to_vec()
    }
}

impl<const N: usize> FromStr for BigInt<N> {
    type Err = BigIntParseErr;

    fn from_str(text: &str) -> Result<Self, Self::Err> {
        // big-endian bytes array
        let bound_be_bytes = {
            let bound = u64::MAX.to_string();
            bound.as_bytes().to_vec()
        };
        let text_be_bytes = text.as_bytes().to_vec();

        // functionality for converting a byte string into an u64 integer
        let extract_one_word = |bytes: &[u8]| {
            bytes
                .iter()
                .fold(0_u64, |acc, &v| (acc * 10 + (v as u64 - '0' as u64)))
        };

        let (bound_size, text_size) = (bound_be_bytes.len() - 1, text_be_bytes.len());
        let max_bound = 10_u64.pow(bound_size as u32);

        // initialize with the highest residual word
        let residual_len = text_size % bound_size;
        let first_word_len = if residual_len == 0 {
            bound_size
        } else {
            residual_len
        };
        let mut result_words: Vec<u64> = vec![extract_one_word(text[..first_word_len].as_bytes())];
        let mut pos = first_word_len;
        while pos < text_size {
            // reserve a empty position for single-precision multiplication and addition
            if result_words.last() != Some(&0_u64) {
                result_words.push(0_u64);
            }
            // single-precision multiplication, left-shifting a max_bound
            let mut carrier = 0_u64;
            for word in result_words.iter_mut() {
                let acc = (*word as u128) * (max_bound as u128) + carrier as u128;
                *word = (acc & ((1_u128 << 64) - 1)) as u64;
                carrier = (acc >> 64) as u64;
            }
            if carrier != 0_u64 {
                return Err(BigIntParseErr);
            }
            // single-precision addition
            let cur_word = extract_one_word(text[pos..(pos + bound_size)].as_bytes());
            carrier = 0_u64;
            for (i, word) in result_words.iter_mut().enumerate() {
                let acc = if i == 0 {
                    (*word) as u128 + cur_word as u128
                } else {
                    (*word) as u128 + carrier as u128
                };
                *word = acc as u64;
                carrier = (acc >> 64) as u64;
            }
            if carrier != 0_u64 {
                return Err(BigIntParseErr);
            }
            pos = pos + bound_size;
        }

        // padding zeros
        if result_words.len() > N {
            return Err(BigIntParseErr);
        }
        let padding_words = vec![0_u64; N - result_words.len()];
        result_words.extend(padding_words);

        Ok(BigInt(result_words.try_into().unwrap()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt as TestBigInt;

    #[test]
    fn test_from_str() {
        let s = "28948022309329048855892746252171976963363056481941647379679742748393362948097";
        let result: Vec<u64> = BigInt::<4>::from_str(s).unwrap().0.try_into().unwrap();
        assert_eq!(result, TestBigInt::from_str(s).unwrap().to_u64_digits().1);
    }
}
