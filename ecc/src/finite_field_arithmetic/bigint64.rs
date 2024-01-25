////// Implementation of Basic Arithmetics for Generalized (variable length) BigInteger<T>
///
use std::ops::{Add, Div, Mul, Rem, Shl, Shr, Sub};

use crate::finite_field_arithmetic::BigInteger;

const WORD_SIZE: usize = 64;
const DOUBLE_WORD_SIZE: usize = 64 * 2;
type Word = u64;
type DoubleWord = u128;
const WORD_BASE: DoubleWord = (1 as DoubleWord) << WORD_SIZE;

fn legendre_symbol(scalar: &BigInteger<Word>, modulus: &BigInteger<Word>) -> bool {
    let mut k = 1;
    let (mut p, mut a) = (modulus.clone(), scalar.clone());
    while p != BigInteger::<Word>::ONE() {
        if a.is_zero() {
            return false;
        }
        let mut v = 0;
        while a.is_even() {
            a = &a >> (1 as usize);
            v = v + 1;
        }
        if (v % 2 == 1) && ((p.data[0] % 8 == 3) || (p.data[0] % 8 == 5)) {
            k = -k;
        }
        if (a.data[0] % 4 == 3) && (p.data[0] % 4 == 3) {
            k = -k;
        }
        (a, p) = (&p % &a, a)
    }
    k == 1
}

fn div_internal(
    lft: &BigInteger<Word>,
    rht: &BigInteger<Word>,
) -> (BigInteger<Word>, BigInteger<Word>) {
    assert!(lft.basis == rht.basis);
    let b = WORD_BASE;

    let (mut u, mut v) = (lft.clone(), rht.clone());
    let (nu, nv) = (u.data.len(), v.data.len());
    assert!(nu >= nv);

    // normalization
    let mut d = 1 as Word;
    while v.data[nv - 1] < (b / 2) as Word {
        v = &v << (1 as usize);
        u = &u << (1 as usize);
        d = d * 2;
    }
    // println!("---- normalization d = {}", d);

    // reset highest digit
    u.data.push(0 as Word);

    // proximate quotient
    let mut q = BigInteger {
        data: vec![0 as Word; nu - nv + 1],
        sign: false,
        basis: lft.basis,
    };
    // println!("u = {:?}, v = {:?}", u, v);
    // i = 0,1
    for i in (0..(nu - nv + 1)).rev() {
        // two word proximate quotient
        let u_2w =
            (DoubleWord::from(u.data[nv + i]) << WORD_SIZE) + DoubleWord::from(u.data[nv + i - 1]);
        let v_1w = DoubleWord::from(v.data[nv - 1]);
        let mut q_prox = (std::cmp::min(u_2w / v_1w, (b - 1) as DoubleWord)) as Word;
        // println!("+++++ q_prox = {} = {:?}/{:?}", q_prox, u_2w, v_1w);

        // three word proximate quotient
        if v.size() >= 2 {
            let v_2w = BigInteger {
                data: v.data[nv - 2..nv].to_vec(),
                sign: v.sign,
                basis: lft.basis,
            };
            let mut u_3w_prox = &v_2w * (q_prox as usize);
            let u_3w = BigInteger {
                data: u.data[nv + i - 2..nv + i + 1].to_vec(),
                sign: u.sign,
                basis: lft.basis,
            };
            // println!("+++++ q_prox = {}, {:?} > {:?}", q_prox, u_3w_prox, u_3w);
            while (&u_3w_prox - &u_3w).is_positive() {
                q_prox = q_prox - 1;
                u_3w_prox = &v_2w * (q_prox as usize);
            }
            // println!("+++++ q_prox = {}, {:?} < {:?}", q_prox, u_3w_prox, u_3w);
        }

        // n_plus_one word proximate quotient
        let u_nplus1w_prox = &v * (q_prox as usize);
        let u_nplus1w = BigInteger {
            data: u.data[i..(i + nv + 1)].to_vec(),
            sign: u.sign,
            basis: lft.basis,
        };
        let mut diff = &u_nplus1w - &u_nplus1w_prox;
        // println!(
        //     "++++ diff = {:?} = {:?} < {:?}",
        //     diff, u_nplus1w, u_nplus1w_prox
        // );
        if diff.is_negative() {
            q_prox = q_prox - 1;
            diff = &diff + &v;
        }
        // println!("++++ diff = {:?}, q_prox = {}", diff, q_prox);
        let updated_upper: Vec<Word> = diff
            .data
            .clone()
            .into_iter()
            .chain(vec![0 as Word; nv + 1 - diff.size()].into_iter())
            .collect();
        u.data[i..(i + nv + 1)].copy_from_slice(updated_upper.as_slice());
        q.data[i] = q_prox;
    }
    q.strip_leading_zeros();

    let mut r = &u / (d as usize);
    r.strip_leading_zeros();

    (q, r)
}

fn div_single_precision_internal(lft: &BigInteger<Word>, rht: usize) -> BigInteger<Word> {
    assert!(rht != 0);
    let (b, u, nu, v) = (&lft.basis, &lft.data, lft.size(), rht as Word);
    let mut w = BigInteger {
        data: vec![0 as Word; nu - 1 + 1],
        sign: lft.sign,
        basis: *b,
    };

    let mut carrier = 0 as Word;
    for i in (0..nu).rev() {
        let (carrier_word, unit_word, nominator_word, denominator_word) = (
            DoubleWord::from(carrier),
            *b as DoubleWord,
            DoubleWord::from(u[i]),
            DoubleWord::from(v),
        );
        let t = carrier_word * unit_word + nominator_word;
        let quotient = t / denominator_word;
        carrier = (t - quotient * denominator_word) as Word;
        w.data[i] = quotient as Word;
    }
    w.strip_leading_zeros();
    w
}

fn shl_internal(lft: &BigInteger<Word>, rht: usize) -> BigInteger<Word> {
    let n = lft.size();
    let mut w = BigInteger {
        data: vec![0 as Word; n + 1],
        sign: false,
        basis: lft.basis,
    };
    let mut carrier = 0 as Word;
    // from the lowest word to the highest word
    for i in 0..n {
        let remainder = (lft.data[i] << rht) | carrier;
        carrier = lft.data[i] >> (WORD_SIZE - rht);
        w.data[i] = remainder;
    }
    if carrier > 0 {
        w.data[n] = carrier;
    }
    w.strip_leading_zeros();
    w
}

fn shr_internal(lft: &BigInteger<Word>, rht: usize) -> BigInteger<Word> {
    let n = lft.size();
    let mut w = BigInteger {
        data: vec![0 as Word; n + 1],
        sign: false,
        basis: lft.basis,
    };
    let mut carrier = 0 as Word;
    // from the hightest word to the lowest word
    for i in (0..n).rev() {
        let remainder = (lft.data[i] >> rht) | carrier;
        carrier = lft.data[i] << (WORD_SIZE - rht);
        w.data[i] = remainder;
    }
    w.strip_leading_zeros();
    w
}

fn mul_single_precision_internal(lft: &BigInteger<Word>, rht: usize) -> BigInteger<Word> {
    let (b, u, nu, v) = (&lft.basis, &lft.data, lft.size(), rht as Word);
    let mut w = BigInteger {
        data: vec![0 as Word; nu + 1],
        sign: lft.sign,
        basis: *b,
    };
    if v != 0 {
        let mut carrier = 0 as Word;
        let mut remainder = 0 as Word;
        for i in 0..nu {
            // double precision for single precision multiplication with carriers
            let (w_word, lft_word, rht_word, carrier_word) = (
                DoubleWord::from(w.data[i]),
                DoubleWord::from(u[i]),
                DoubleWord::from(v),
                DoubleWord::from(carrier),
            );
            let t = w_word + lft_word * rht_word + carrier_word;
            (carrier, remainder) = (
                (t >> WORD_SIZE) as Word,
                (t & ((1 << WORD_SIZE) - 1)) as Word,
            );

            w.data[i] = remainder;
        }
        if carrier > 0 {
            w.data[nu] = carrier;
        }
    }
    w.strip_leading_zeros();
    w
}

fn mul_internal(lft: &BigInteger<Word>, rht: &BigInteger<Word>) -> BigInteger<Word> {
    assert!(lft.basis == rht.basis);
    let b = &lft.basis;
    let (u, v) = (&lft.data, &rht.data);
    let (nu, nv) = (u.len(), v.len());

    let mut w = BigInteger {
        data: vec![0 as Word; nu + nv],
        sign: false,
        basis: *b,
    };

    // initialization
    for i in 0..nv {
        w.data[i] = 0 as Word;
    }
    // cross multiplication
    for i in 0..nv {
        let mut carrier = 0 as Word;
        let mut remainder = 0 as Word;
        if v[i] != 0 as Word {
            for j in 0..nu {
                // double precision for single precision multiplication with carriers
                let (w_word, lft_word, rht_word, carrier_word) = (
                    DoubleWord::from(w.data[i + j]),
                    DoubleWord::from(u[j]),
                    DoubleWord::from(v[i]),
                    DoubleWord::from(carrier),
                );
                let t = w_word + lft_word * rht_word + carrier_word;
                (carrier, remainder) = (
                    (t >> WORD_SIZE) as Word,
                    (t & ((1 << WORD_SIZE) - 1)) as Word,
                );

                w.data[i + j] = remainder;
            }
            if carrier > 0 {
                w.data[nu + i] = carrier;
            }
            // println!("--- {} * {:?} = {:?}", v[i], u, w);
        }
    }
    w.strip_leading_zeros();
    w
}

fn add_internal(lft: &BigInteger<Word>, rht: &BigInteger<Word>) -> BigInteger<Word> {
    assert!(lft.basis == rht.basis);
    let b = &lft.basis;
    let (u, v) = (&lft.data, &rht.data);
    let (nu, nv) = (u.len(), v.len());

    let max_n = std::cmp::max(nu, nv);
    let mut w = BigInteger {
        data: vec![0 as Word; max_n + 1],
        sign: false,
        basis: *b,
    };
    let mut carrier = 0 as Word;
    let mut remainder = 0 as Word;
    // from right to left
    for i in 0..max_n {
        let lft = if i <= nu - 1 { u[i] } else { 0 as Word };
        let rht = if i <= nv - 1 { v[i] } else { 0 as Word };

        // double precision for single precision addition with carriers
        let result = DoubleWord::from(lft) + DoubleWord::from(rht) + DoubleWord::from(carrier);
        (carrier, remainder) = (
            (result >> WORD_SIZE) as Word,
            (result & ((1 << WORD_SIZE) - 1)) as Word,
        );

        w.data[i] = remainder;
    }
    // the highest digit must be non-zero
    if carrier > 0 {
        w.data[max_n] = carrier;
    }
    w.strip_leading_zeros();
    w
}

fn sub_internal(lft: &BigInteger<Word>, rht: &BigInteger<Word>) -> BigInteger<Word> {
    assert!(lft.basis == rht.basis);
    if lft.is_zero() {
        return BigInteger {
            data: rht.data.clone(),
            sign: !rht.sign,
            basis: rht.basis,
        };
    }

    let b = &lft.basis;
    let (u, v) = (&lft.data, &rht.data);
    let (nu, nv) = (u.len(), v.len());

    let max_n = std::cmp::max(nu, nv);
    let mut w = BigInteger {
        data: vec![0 as Word; max_n],
        sign: false,
        basis: *b,
    };
    let mut carrier = 0 as Word;
    let mut remainder = 0 as Word;
    // from right to left
    for i in 0..max_n {
        let lft = if i <= nu - 1 { u[i] } else { 0 as Word };
        let rht = if i <= nv - 1 { v[i] } else { 0 as Word };

        // double precision for single precision substraction with carriers
        let result = DoubleWord::from(lft)
            .wrapping_sub(DoubleWord::from(rht))
            .wrapping_sub(DoubleWord::from(carrier));
        (carrier, remainder) = (
            ((result >> WORD_SIZE) as Word).wrapping_neg(),
            (result & ((1 << WORD_SIZE) - 1)) as Word,
        );

        w.data[i] = remainder;
    }
    w.strip_leading_zeros();
    // we need to reverse the digits where the result has carrier
    if carrier > 0 {
        for i in 0..w.size() {
            if i == 0 {
                w.data[i] = (WORD_BASE - w.data[i] as DoubleWord) as Word;
            } else {
                w.data[i] = Word::MAX - w.data[i];
            }
        }
        w.sign = true;
    }
    w.strip_leading_zeros();
    w
}

//// &BigInteger + &BigInteger
impl<'a, 'b> Add<&'b BigInteger<Word>> for &'a BigInteger<Word> {
    type Output = BigInteger<Word>;

    fn add(self, other: &'b BigInteger<Word>) -> BigInteger<Word> {
        if self.is_zero() {
            return other.clone();
        }
        if other.is_zero() {
            return self.clone();
        }

        if self.sign ^ other.sign {
            let (lft, rht) = if self.sign == true {
                (other, self)
            } else {
                (self, other)
            };
            sub_internal(lft, rht)
        } else {
            let mut result = add_internal(self, other);
            result.sign = self.sign & other.sign;
            result
        }
    }
}

impl Add<usize> for &BigInteger<Word> {
    type Output = BigInteger<Word>;

    fn add(self, other: usize) -> BigInteger<Word> {
        assert!((other as DoubleWord) < (1 as DoubleWord) << WORD_SIZE);
        self + &(BigInteger::new(vec![other as Word].as_slice(), false, self.basis, None))
    }
}

impl<'a, 'b> Sub<&'b BigInteger<Word>> for &'a BigInteger<Word> {
    type Output = BigInteger<Word>;

    fn sub(self, other: &'b BigInteger<Word>) -> BigInteger<Word> {
        let mut other_new = other.clone();
        other_new.sign = !other_new.sign;
        self.add(&other_new)
    }
}

impl<'a, 'b> Mul<&'b BigInteger<Word>> for &'a BigInteger<Word> {
    type Output = BigInteger<Word>;

    fn mul(self, other: &'b BigInteger<Word>) -> BigInteger<Word> {
        let mut result = mul_internal(self, other);
        result.sign = self.sign & other.sign;
        result
    }
}

impl Mul<usize> for &BigInteger<Word> {
    type Output = BigInteger<Word>;

    fn mul(self, other: usize) -> BigInteger<Word> {
        assert!((other as DoubleWord) < (1 as DoubleWord) << WORD_SIZE);
        let mut result = mul_single_precision_internal(self, other);
        result.sign = self.sign;
        result
    }
}

impl Shl<usize> for &BigInteger<Word> {
    type Output = BigInteger<Word>;

    fn shl(self, rhs: usize) -> Self::Output {
        assert!(rhs < WORD_SIZE);
        let mut result = shl_internal(self, rhs);
        result.sign = self.sign;
        result
    }
}

impl Shr<usize> for &BigInteger<Word> {
    type Output = BigInteger<Word>;

    fn shr(self, rhs: usize) -> Self::Output {
        assert!(rhs < WORD_SIZE);
        let mut result = shr_internal(self, rhs);
        result.sign = self.sign;
        result
    }
}

impl Div<usize> for &BigInteger<Word> {
    type Output = BigInteger<Word>;

    fn div(self, other: usize) -> BigInteger<Word> {
        assert!((other as DoubleWord) < (1 as DoubleWord) << WORD_SIZE);
        let mut result = div_single_precision_internal(self, other);
        result.sign = self.sign;
        result
    }
}

impl<'a, 'b> Div<&'b BigInteger<Word>> for &'a BigInteger<Word> {
    type Output = BigInteger<Word>;

    fn div(self, other: &'b BigInteger<Word>) -> BigInteger<Word> {
        let (mut quotient, _) = div_internal(self, other);
        quotient.sign = self.sign & other.sign;
        quotient
    }
}

impl<'a, 'b> Rem<&'b BigInteger<Word>> for &'a BigInteger<Word> {
    type Output = BigInteger<Word>;

    fn rem(self, other: &'b BigInteger<Word>) -> BigInteger<Word> {
        let (_, mut remainder) = div_internal(self, other);
        remainder.sign = self.sign & other.sign;
        remainder
    }
}

impl Default for BigInteger<Word> {
    fn default() -> Self {
        BigInteger {
            data: vec![Word::default()],
            sign: false,
            basis: WORD_SIZE,
        }
    }
}

impl From<Vec<Word>> for BigInteger<Word> {
    fn from(words: Vec<Word>) -> Self {
        Self {
            data: words,
            sign: false,
            basis: WORD_SIZE,
        }
    }
}

impl From<&str> for BigInteger<Word> {
    fn from(text: &str) -> Self {
        // big-endian bytes array
        let bound_be_bytes = {
            let bound = Word::MAX.to_string();
            bound.as_bytes().to_vec()
        };
        let text_be_bytes = text.as_bytes().to_vec();

        // functionality for converting a byte string into an u64 integer
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
                *word = (acc & (WORD_BASE - 1)) as Word;
                carrier = (acc >> WORD_SIZE) as Word;
            }
            assert!(carrier == (0 as Word));
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
            assert!(carrier == (0 as Word));
            pos = pos + bound_size;
        }

        BigInteger {
            data: result_words.try_into().unwrap(),
            sign: false,
            basis: WORD_SIZE,
        }
    }
}

impl BigInteger<Word> {
    fn ZERO() -> BigInteger<Word> {
        Self {
            data: vec![0 as Word],
            sign: false,
            basis: WORD_SIZE,
        }
    }

    fn ONE() -> BigInteger<Word> {
        Self {
            data: vec![1 as Word],
            sign: false,
            basis: WORD_SIZE,
        }
    }

    fn squre_root(&self) -> BigInteger<Word> {
        if self.is_zero() {
            return self.clone();
        }
        assert!(self.is_positive());

        // initialization
        let n = self.data.len();
        let nr = if n % 2 == 0 { n / 2 } else { (n + 1) / 2 };
        let mut init_data = vec![0 as Word; nr];
        init_data[nr - 1] = 1 as Word;
        let mut t = BigInteger {
            data: init_data,
            sign: false,
            basis: self.basis,
        };
        let mut v = BigInteger {
            data: vec![0 as Word],
            sign: false,
            basis: self.basis,
        };

        // iterate with newton method
        loop {
            v = t;
            let delta = self / &v;
            t = &(&v + &delta) / (2 as usize);
            if t == v {
                break;
            }
        }

        v
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt as TestBigInt;
    use num_bigint::BigUint as TestBigUInt;
    use std::str::FromStr;

    #[test]
    fn test_from() {
        let Fr = "28948022309329048855892746252171976963363056481941647379679742748393362948097";
        let R = "28948022309329048855892746252171976963180815219815621900418355762733040795645";
        assert_eq!(
            TestBigInt::from_str(Fr).unwrap().to_u64_digits().1,
            BigInteger::from(Fr).data
        );
    }

    #[test]
    fn test_add() {
        let Fr = "28948022309329048855892746252171976963363056481941647379679742748393362948097";
        let R = "28948022309329048855892746252171976963180815219815621900418355762733040795645";
        assert_eq!(
            (TestBigInt::from_str(Fr).unwrap() + TestBigInt::from_str(R).unwrap())
                .to_u64_digits()
                .1,
            (&BigInteger::from(Fr) + &BigInteger::from(R)).data
        );
    }

    #[test]
    fn test_mul() {
        let Fr = "28948022309329048855892746252171976963363056481941647379679742748393362948097";
        let R = "28948022309329048855892746252171976963180815219815621900418355762733040795645";
        assert_eq!(
            (TestBigInt::from_str(Fr).unwrap() * TestBigInt::from_str(R).unwrap())
                .to_u64_digits()
                .1,
            (&BigInteger::from(Fr) * &BigInteger::from(R)).data
        );
    }

    #[test]
    fn test_Fr_R() {
        let Fr = "28948022309329048855892746252171976963363056481941647379679742748393362948097";
        let R = "28948022309329048855892746252171976963180815219815621900418355762733040795645";
        assert_eq!(
            (TestBigInt::from_str(Fr).unwrap() / TestBigInt::from_str(R).unwrap())
                .to_u64_digits()
                .1,
            (&BigInteger::from(Fr) / &BigInteger::from(R)).data
        );
    }

    #[test]
    fn test_Fr_shl() {
        let Fr = "28948022309329048855892746252171976963363056481941647379679742748393362948097";
        let R = "28948022309329048855892746252171976963180815219815621900418355762733040795645";
        assert_eq!(
            (TestBigInt::from_str(Fr).unwrap() << 1 as Word)
                .to_u64_digits()
                .1,
            (&BigInteger::from(Fr) << (1 as usize)).data
        );
    }

    #[test]
    fn test_r_div_Fr() {
        let (Fr, r) = (
            "28948022309329048855892746252171976963363056481941647379679742748393362948097",
            "115792089237316195423570985008687907853269984665640564039457584007913129639936",
        );
        let actual_result = (TestBigInt::from_str(r).unwrap() / TestBigInt::from_str(Fr).unwrap())
            .to_u64_digits()
            .1;
        println!(
            "Experimental: {:?} = {:?} / {:?}",
            (&BigInteger::from(r) / &BigInteger::from(Fr)).data,
            &BigInteger::from(r),
            &BigInteger::from(Fr)
        );

        assert_eq!(
            actual_result,
            (&BigInteger::from(r) / &BigInteger::from(Fr)).data
        );
    }

    #[test]
    fn test_Fr_factoring() {
        let Fr = "28948022309329048855892746252171976963363056481941647379679742748393362948097";
        let lft = &(&BigInteger::from(Fr) - &BigInteger::from("1"));
        let tmp_candidates: Vec<Vec<BigInteger<u64>>> = (0..4)
            .map(|i| {
                let start = i * 64;
                (start..(start + 64))
                    .map(|exp| {
                        let mut r = vec![0 as u64; i + 1];
                        r[i] = (1 as u64) << (exp - start);
                        BigInteger::from(r)
                    })
                    .collect()
            })
            .collect();
        let candidates: Vec<BigInteger<u64>> = tmp_candidates.into_iter().flatten().collect();

        for i in 0..candidates.len() {
            let rht = &candidates[i];
            let (q, r) = div_internal(lft, rht);
            if r.is_zero() && (q.is_even() == false) {
                println!("e = {}, s = {:?}, r = {:?}", i, q, r);
                assert_eq!(&(&q * rht) + &BigInteger::from("1"), BigInteger::from(Fr));
            }
        }
    }

    #[test]
    fn test_legendre_symbol() {
        let Fr = "28948022309329048855892746252171976963363056481941647379679742748393362948097";
        for v in 1..10 {
            let (a, b) = (
                BigInteger::from(v.to_string().as_str()),
                BigInteger::from(Fr),
            );

            let sym = legendre_symbol(&a, &b);
            println!("{:?} quadratic residual? {}", &a, sym);
        }
    }
}
