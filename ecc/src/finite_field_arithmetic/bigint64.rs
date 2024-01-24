use std::cmp::Ordering;
use std::ops::{Add, Div, Mul, Neg, Rem, Shl, Shr, Sub};
use std::str::FromStr;

#[derive(Debug)]
pub struct ParseStrErr;

//////////////////////////////////// Implementation of BigInteger Specially for Finite Field
#[derive(Debug, PartialEq, Clone, Copy, Eq)]
pub struct BigInt<const N: usize>(pub [u64; N]);

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

impl<const N: usize> FromStr for BigInt<N> {
    type Err = ParseStrErr;

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
                return Err(ParseStrErr);
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
                return Err(ParseStrErr);
            }
            pos = pos + bound_size;
        }

        // padding zeros
        if result_words.len() > N {
            return Err(ParseStrErr);
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
