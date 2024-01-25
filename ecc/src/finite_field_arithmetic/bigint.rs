//////////////////////////////////// Practical Implementation of BigInteger Specially for Finite Field (fixed length)
///
///
use std::cmp::Ordering;
use std::ops::{Add, Mul, Shl, Shr, Sub};
use std::str::FromStr;

#[derive(Debug)]
pub struct BigIntParseErr;

const WORD_SIZE: usize = 64;
const DOUBLE_WORD_SIZE: usize = 64 * 2;
type Word = u64;
type DoubleWord = u128;
const WORD_BASE: DoubleWord = (1 as DoubleWord) << WORD_SIZE;

pub trait BigInteger<const N: usize>:
    FromStr + PartialOrd + Into<Vec<u8>> + Ord + Add + Mul<Word> + Shl<usize> + Shr<usize> + Sub
{
}
impl<const N: usize> BigInteger<N> for BigInt<N> {}

#[derive(Debug, PartialEq, Clone, Copy, Eq)]
pub struct BigInt<const N: usize>(pub [Word; N]);

impl<const N: usize> BigInt<N> {
    #[inline(always)]
    pub fn MAX() -> Self {
        Self(
            vec![(((1 as DoubleWord) << WORD_SIZE) - 1) as Word; N]
                .try_into()
                .unwrap(),
        )
    }

    #[inline(always)]
    pub fn ZERO() -> Self {
        BigInt([0 as Word; N])
    }

    #[inline(always)]
    pub fn ONE() -> Self {
        let mut data = vec![0 as Word; N];
        data[0] = 1 as Word;
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
            zero = zero & (self.0[i] == 0 as Word);
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
    type Output = (BigInt<N>, Word);
    fn shl(self, other: usize) -> Self::Output {
        let n = N;
        assert!(other < WORD_SIZE);
        let mut w = Self([0 as Word; N]);
        let mut carrier = 0 as Word;
        // from the lowest word to the highest word
        for i in 0..n {
            let tmp_carrier = self.0[i] >> (WORD_SIZE - other);
            let remainder = (self.0[i] << other) | carrier;
            w.0[i] = remainder;
            carrier = tmp_carrier;
        }
        (w, carrier)
    }
}

impl<const N: usize> Shr<usize> for BigInt<N> {
    type Output = (BigInt<N>, Word);
    fn shr(self, other: usize) -> Self::Output {
        let n = N;
        assert!(other < WORD_SIZE);
        let mut w = Self([0 as Word; N]);
        let mut carrier = 0 as Word;
        // from the lowest word to the highest word
        for i in (0..n).rev() {
            let tmp_carrier = self.0[i] & (((1 << other) - 1) as Word);
            let remainder = (self.0[i] >> other) | (carrier << (WORD_SIZE - other));
            w.0[i] = remainder;
            carrier = tmp_carrier;
        }
        (w, carrier)
    }
}

impl<const N: usize> Add for BigInt<N> {
    type Output = (BigInt<N>, Word);
    fn add(self, other: BigInt<N>) -> Self::Output {
        let n = N;
        let lft = &self.0;
        let rht = &other.0;
        let (mut carrier, mut remainder) = (0 as Word, 0 as Word);
        let mut w = BigInt([0 as Word; N]);
        for i in 0..n {
            let (lft_w, rht_w, carrier_w) = (
                DoubleWord::from(lft[i]),
                DoubleWord::from(rht[i]),
                DoubleWord::from(carrier),
            );
            let t = lft_w + (rht_w + carrier_w);
            (carrier, remainder) = (
                (t >> WORD_SIZE) as Word,
                (t & ((1 << WORD_SIZE) - 1)) as Word,
            );
            w.0[i] = remainder;
        }
        (w, carrier)
    }
}

impl<const N: usize> Sub for BigInt<N> {
    type Output = (BigInt<N>, Word);
    fn sub(self, other: BigInt<N>) -> Self::Output {
        let n = N;
        let lft = &self.0;
        let rht = &other.0;
        let (mut carrier, mut remainder) = (0 as Word, 0 as Word);
        let mut w = BigInt([0 as Word; N]);
        for i in 0..n {
            let (lft_w, rht_w, carrier_w) = (
                DoubleWord::from(lft[i]),
                DoubleWord::from(rht[i]),
                DoubleWord::from(carrier),
            );
            (carrier, remainder) = if lft_w >= rht_w + carrier_w {
                (0 as Word, (lft_w - (rht_w + carrier_w)) as Word)
            } else {
                (
                    1 as Word,
                    (lft_w + (1 << WORD_SIZE) - (rht_w + carrier_w)) as Word,
                )
            };
            w.0[i] = remainder;
        }
        (w, carrier)
    }
}

impl<const N: usize> Mul<Word> for BigInt<N> {
    type Output = (BigInt<N>, Word);

    fn mul(self, other: Word) -> Self::Output {
        let n = N;
        let lft = &self.0;
        let (mut carrier, mut remainder) = (0 as Word, 0 as Word);
        let mut w = BigInt([0 as Word; N]);
        for i in 0..n {
            let (lft_w, rht_w, carrier_w) = (
                DoubleWord::from(lft[i]),
                DoubleWord::from(other),
                DoubleWord::from(carrier),
            );
            let t = lft_w * rht_w + carrier_w;
            (carrier, remainder) = (
                (t >> WORD_SIZE) as Word,
                (t & ((1 << WORD_SIZE) - 1)) as Word,
            );
            w.0[i] = remainder;
        }
        (w, carrier)
    }
}

impl<const N: usize> Into<Vec<u8>> for BigInt<N> {
    fn into(self) -> Vec<u8> {
        let num_leading_zeros = self.0.iter().rev().fold(0 as u32, |acc, v| {
            if acc == 0 {
                acc + v.leading_zeros()
            } else if (acc % (WORD_SIZE as u32)) == 0 {
                acc + v.leading_zeros()
            } else {
                acc
            }
        });
        let bits = (0..self.0.len())
            .map(|i| ((0..WORD_SIZE).map(|j| ((self.0[i] >> j) & 1) as u8)).collect::<Vec<u8>>())
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
            let bound = Word::MAX.to_string();
            bound.as_bytes().to_vec()
        };
        let text_be_bytes = text.as_bytes().to_vec();

        // functionality for converting a byte string into an Word integer
        let extract_one_word = |bytes: &[u8]| {
            bytes
                .iter()
                .fold(0 as Word, |acc, &v| (acc * 10 + (v as Word - '0' as Word)))
        };

        let (bound_size, text_size) = (bound_be_bytes.len() - 1, text_be_bytes.len());
        let max_bound = (10 as Word).pow(bound_size as u32);

        // initialize with the highest residual word
        let residual_len = text_size % bound_size;
        let first_word_len = if residual_len == 0 {
            bound_size
        } else {
            residual_len
        };
        let mut result_words: Vec<Word> = vec![extract_one_word(text[..first_word_len].as_bytes())];
        let mut pos = first_word_len;
        while pos < text_size {
            // reserve a empty position for single-precision multiplication and addition
            if result_words.last() != Some(&(0 as Word)) {
                result_words.push(0 as Word);
            }
            // single-precision multiplication, left-shifting a max_bound
            let mut carrier = 0 as Word;
            for word in result_words.iter_mut() {
                let acc = (*word as DoubleWord) * (max_bound as DoubleWord) + carrier as DoubleWord;
                *word = (acc & (((1 as DoubleWord) << WORD_SIZE) - 1)) as Word;
                carrier = (acc >> WORD_SIZE) as Word;
            }
            if carrier != 0 as Word {
                return Err(BigIntParseErr);
            }
            // single-precision addition
            let cur_word = extract_one_word(text[pos..(pos + bound_size)].as_bytes());
            carrier = 0 as Word;
            for (i, word) in result_words.iter_mut().enumerate() {
                let acc = if i == 0 {
                    (*word) as DoubleWord + cur_word as DoubleWord
                } else {
                    (*word) as DoubleWord + carrier as DoubleWord
                };
                *word = acc as Word;
                carrier = (acc >> WORD_SIZE) as Word;
            }
            if carrier != 0 as Word {
                return Err(BigIntParseErr);
            }
            pos = pos + bound_size;
        }

        // padding zeros
        if result_words.len() > N {
            return Err(BigIntParseErr);
        }
        let padding_words = vec![0 as Word; N - result_words.len()];
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
        let Fr = "28948022309329048855892746252171976963363056481941647379679742748393362948097";
        let R = "28948022309329048855892746252171976963180815219815621900418355762733040795645";
        let R2 = "4263855311957679929489659445116329028194309752796460188622876710448966664207";
        let R3 = "3557709630315679472311684181007729646594247341237824434526702614836137537100";

        let FrVec: Vec<Word> = BigInt::<4>::from_str(Fr).unwrap().0.try_into().unwrap();
        let RVec: Vec<Word> = BigInt::<4>::from_str(R).unwrap().0.try_into().unwrap();
        let R2Vec: Vec<Word> = BigInt::<4>::from_str(R2).unwrap().0.try_into().unwrap();
        let R3Vec: Vec<Word> = BigInt::<4>::from_str(R3).unwrap().0.try_into().unwrap();

        println!(
            "pallas fr = {:?} \n R = {:?} \n R2 = {:?} \n R3 = {:?}\n",
            FrVec, RVec, R2Vec, R3Vec
        );

        assert_eq!(FrVec, TestBigInt::from_str(Fr).unwrap().to_u64_digits().1);
        assert_eq!(RVec, TestBigInt::from_str(R).unwrap().to_u64_digits().1);
        assert_eq!(R2Vec, TestBigInt::from_str(R2).unwrap().to_u64_digits().1);
        assert_eq!(R3Vec, TestBigInt::from_str(R3).unwrap().to_u64_digits().1);
    }

    #[test]
    fn test_addition() {}

    #[test]
    fn test_substraction() {}

    #[test]
    fn test_multiplication() {
        let (a, b) = (
            "10828745280282393011948633936436145363160692580455384354038716315440980557097",
            "18200867980676431887",
        );
        let lft = BigInt::<4>::from_str(a).unwrap();
        let rht = Word::from_str(b).unwrap();
        let result = lft * rht;
        println!("{:?}, {} = {:?} * {}", result.0 .0, result.1, lft.0, rht);
    }

    #[test]
    fn test_division() {}
}
