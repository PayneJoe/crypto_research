////// Implementation of Basic Arithmetics for BigInteger<u16>
///
use std::ops::{Add, Div, Mul, Shl, Shr, Sub};

use crate::finite_field_arithmetic::BigInteger;

fn div_internal(
    lft: &BigInteger<u16>,
    rht: &BigInteger<u16>,
) -> (BigInteger<u16>, BigInteger<u16>) {
    assert!(lft.basis == rht.basis);
    let b = lft.basis;

    let (mut u, mut v) = (lft.clone(), rht.clone());
    let (nu, nv) = (u.data.len(), v.data.len());
    assert!(nu >= nv);

    // normalization
    let mut d = 1_u16;
    while v.data[nv - 1] < (b / 2) as u16 {
        v = &v << (1 as usize);
        u = &u << (1 as usize);
        d = d * 2;
    }

    // reset highest digit
    u.data.push(0 as u16);

    // proximate quotient
    let mut q = BigInteger {
        data: vec![0_u16; nu - nv + 1],
        sign: false,
        basis: b,
    };
    for i in (0..(nu - nv + 1)).rev() {
        // two word proximate quotient
        let u_2w = (u32::from(u.data[nv + i]) << 16) + u32::from(u.data[nv + i - 1]);
        let v_1w = u32::from(v.data[nv - 1]);
        let mut q_prox = (std::cmp::min(u_2w / v_1w, (b - 1) as u32)) as u16;

        // three word proximate quotient
        if v.size() >= 2 {
            let v_2w = BigInteger {
                data: v.data[nv - 2..nv].to_vec(),
                sign: v.sign,
                basis: b,
            };
            let mut u_3w_prox = &v_2w * (q_prox as usize);
            let u_3w = BigInteger {
                data: u.data[nv + i - 2..nv + i + 1].to_vec(),
                sign: u.sign,
                basis: b,
            };
            while (&u_3w_prox - &u_3w).is_positive() {
                q_prox = q_prox - 1;
                u_3w_prox = &v_2w * (q_prox as usize);
            }
        }

        // n_plus_one word proximate quotient
        let u_nplus1w_prox = &v * (q_prox as usize);
        let u_nplus1w = BigInteger {
            data: u.data[i..(i + nv + 1)].to_vec(),
            sign: u.sign,
            basis: b,
        };
        let mut diff = &u_nplus1w - &u_nplus1w_prox;
        if diff.is_negative() {
            diff = &diff + &v;
        }
        let updated_upper: Vec<u16> = diff
            .data
            .clone()
            .into_iter()
            .chain(vec![0_u16; nv + 1 - diff.size()].into_iter())
            .collect();
        u.data[i..(i + nv + 1)].copy_from_slice(updated_upper.as_slice());
        q.data[i] = q_prox;
    }
    q.strip_leading_zeros();

    let r = &u / (d as usize);
    (q, r)
}

fn div_single_precision_internal(lft: &BigInteger<u16>, rht: usize) -> BigInteger<u16> {
    assert!(rht != 0);
    let (b, u, nu, v) = (&lft.basis, &lft.data, lft.size(), rht as u16);
    let mut w = BigInteger {
        data: vec![0 as u16; nu - 1 + 1],
        sign: lft.sign,
        basis: *b,
    };

    let mut carrier = 0_u16;
    for i in (0..nu).rev() {
        let (carrier_word, unit_word, nominator_word, denominator_word) =
            (u32::from(carrier), *b as u32, u32::from(u[i]), u32::from(v));
        let t = carrier_word * unit_word + nominator_word;
        let quotient = t / denominator_word;
        carrier = (t - quotient * denominator_word) as u16;
        w.data[i] = quotient as u16;
    }
    w.strip_leading_zeros();
    w
}

fn shl_internal(lft: &BigInteger<u16>, rht: usize) -> BigInteger<u16> {
    let n = lft.size();
    let mut w = BigInteger {
        data: vec![0_u16; n + 1],
        sign: false,
        basis: lft.basis,
    };
    let mut carrier = 0_u16;
    // from the lowest word to the highest word
    for i in 0..n {
        let remainder = (lft.data[i] << rht) | carrier;
        carrier = lft.data[i] >> (16 - rht);
        w.data[i] = remainder;
    }
    if carrier > 0 {
        w.data[n] = carrier;
    }
    w.strip_leading_zeros();
    w
}

fn shr_internal(lft: &BigInteger<u16>, rht: usize) -> BigInteger<u16> {
    let n = lft.size();
    let mut w = BigInteger {
        data: vec![0_u16; n + 1],
        sign: false,
        basis: lft.basis,
    };
    let mut carrier = 0_u16;
    // from the hightest word to the lowest word
    for i in (0..n).rev() {
        let remainder = (lft.data[i] >> rht) | carrier;
        carrier = lft.data[i] << (16 - rht);
        w.data[i] = remainder;
    }
    w.strip_leading_zeros();
    w
}

fn mul_single_precision_internal(lft: &BigInteger<u16>, rht: usize) -> BigInteger<u16> {
    let (b, u, nu, v) = (&lft.basis, &lft.data, lft.size(), rht as u16);
    let mut w = BigInteger {
        data: vec![0 as u16; nu + 1],
        sign: lft.sign,
        basis: *b,
    };
    if v != 0 {
        let mut carrier = 0_u16;
        let mut remainder = 0_u16;
        for i in 0..nu {
            // double precision for single precision multiplication with carriers
            let (w_word, lft_word, rht_word, carrier_word) = (
                u32::from(w.data[i]),
                u32::from(u[i]),
                u32::from(v),
                u32::from(carrier),
            );
            let t = w_word + lft_word * rht_word + carrier_word;
            (carrier, remainder) = ((t >> 16) as u16, (t & ((1 << 16) - 1)) as u16);

            w.data[i] = remainder;
        }
        if carrier > 0 {
            w.data[nu] = carrier;
        }
    }
    w.strip_leading_zeros();
    w
}

fn mul_internal(lft: &BigInteger<u16>, rht: &BigInteger<u16>) -> BigInteger<u16> {
    assert!(lft.basis == rht.basis);
    let b = &lft.basis;
    let (u, v) = (&lft.data, &rht.data);
    let (nu, nv) = (u.len(), v.len());

    let mut w = BigInteger {
        data: vec![0 as u16; nu + nv],
        sign: false,
        basis: *b,
    };

    // initialization
    for i in 0..nv {
        w.data[i] = 0 as u16;
    }
    // cross multiplication
    for i in 0..nv {
        let mut carrier = 0 as u16;
        let mut remainder = 0 as u16;
        if v[i] != 0_u16 {
            for j in 0..nu {
                // double precision for single precision multiplication with carriers
                let (w_word, lft_word, rht_word, carrier_word) = (
                    u32::from(w.data[i + j]),
                    u32::from(u[j]),
                    u32::from(v[i]),
                    u32::from(carrier),
                );
                let t = w_word + lft_word * rht_word + carrier_word;
                (carrier, remainder) = ((t >> 16) as u16, (t & ((1 << 16) - 1)) as u16);

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

fn add_internal(lft: &BigInteger<u16>, rht: &BigInteger<u16>) -> BigInteger<u16> {
    assert!(lft.basis == rht.basis);
    let b = &lft.basis;
    let (u, v) = (&lft.data, &rht.data);
    let (nu, nv) = (u.len(), v.len());

    let max_n = std::cmp::max(nu, nv);
    let mut w = BigInteger {
        data: vec![0_u16; max_n + 1],
        sign: false,
        basis: *b,
    };
    let mut carrier = 0 as u16;
    let mut remainder = 0 as u16;
    // from right to left
    for i in 0..max_n {
        let lft = if i <= nu - 1 { u[i] } else { 0_u16 };
        let rht = if i <= nv - 1 { v[i] } else { 0_u16 };

        // double precision for single precision addition with carriers
        let result = u32::from(lft) + u32::from(rht) + u32::from(carrier);
        (carrier, remainder) = ((result >> 16) as u16, (result & ((1 << 16) - 1)) as u16);

        w.data[i] = remainder;
    }
    // the highest digit must be non-zero
    if carrier > 0 {
        w.data[max_n] = carrier;
    }
    w.strip_leading_zeros();
    w
}

fn sub_internal(lft: &BigInteger<u16>, rht: &BigInteger<u16>) -> BigInteger<u16> {
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
        data: vec![0_u16; max_n],
        sign: false,
        basis: *b,
    };
    let mut carrier = 0 as u16;
    let mut remainder = 0 as u16;
    // from right to left
    for i in 0..max_n {
        let lft = if i <= nu - 1 { u[i] } else { 0_u16 };
        let rht = if i <= nv - 1 { v[i] } else { 0_u16 };

        // double precision for single precision substraction with carriers
        let result = u32::from(lft)
            .wrapping_sub(u32::from(rht))
            .wrapping_sub(u32::from(carrier));
        (carrier, remainder) = (
            ((result >> 16) as u16).wrapping_neg(),
            (result & ((1 << 16) - 1)) as u16,
        );

        w.data[i] = remainder;
    }
    w.strip_leading_zeros();
    // we need to reverse the digits where the result has carrier
    if carrier > 0 {
        for i in 0..w.size() {
            if i == 0 {
                w.data[i] = ((*b as u32) - w.data[i] as u32) as u16;
            } else {
                w.data[i] = u16::MAX - w.data[i];
            }
        }
        w.sign = true;
    }
    w
}

//// &BigInteger + &BigInteger
impl<'a, 'b> Add<&'b BigInteger<u16>> for &'a BigInteger<u16> {
    type Output = BigInteger<u16>;

    fn add(self, other: &'b BigInteger<u16>) -> BigInteger<u16> {
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

impl Add<usize> for &BigInteger<u16> {
    type Output = BigInteger<u16>;

    fn add(self, other: usize) -> BigInteger<u16> {
        assert!(other < 1 << 16);
        self + &(BigInteger::new(vec![other as u16].as_slice(), false, self.basis, None))
    }
}

impl<'a, 'b> Sub<&'b BigInteger<u16>> for &'a BigInteger<u16> {
    type Output = BigInteger<u16>;

    fn sub(self, other: &'b BigInteger<u16>) -> BigInteger<u16> {
        let mut other_new = other.clone();
        other_new.sign = !other_new.sign;
        self.add(&other_new)
    }
}

impl<'a, 'b> Mul<&'b BigInteger<u16>> for &'a BigInteger<u16> {
    type Output = BigInteger<u16>;

    fn mul(self, other: &'b BigInteger<u16>) -> BigInteger<u16> {
        let mut result = mul_internal(self, other);
        result.sign = self.sign & other.sign;
        result
    }
}

impl Mul<usize> for &BigInteger<u16> {
    type Output = BigInteger<u16>;

    fn mul(self, other: usize) -> BigInteger<u16> {
        assert!(other < 1 << 16);
        let mut result = mul_single_precision_internal(self, other);
        result.sign = self.sign;
        result
    }
}

impl Shl<usize> for &BigInteger<u16> {
    type Output = BigInteger<u16>;

    fn shl(self, rhs: usize) -> Self::Output {
        assert!(rhs < 16);
        let mut result = shl_internal(self, rhs);
        result.sign = self.sign;
        result
    }
}

impl Shr<usize> for &BigInteger<u16> {
    type Output = BigInteger<u16>;

    fn shr(self, rhs: usize) -> Self::Output {
        assert!(rhs < 16);
        let mut result = shr_internal(self, rhs);
        result.sign = self.sign;
        result
    }
}

impl Div<usize> for &BigInteger<u16> {
    type Output = BigInteger<u16>;

    fn div(self, other: usize) -> BigInteger<u16> {
        assert!(other < 1 << 16);
        let mut result = div_single_precision_internal(self, other);
        result.sign = self.sign;
        result
    }
}

impl<'a, 'b> Div<&'b BigInteger<u16>> for &'a BigInteger<u16> {
    type Output = BigInteger<u16>;

    fn div(self, other: &'b BigInteger<u16>) -> BigInteger<u16> {
        let (mut quotient, _) = div_internal(self, other);
        quotient.sign = self.sign & other.sign;
        quotient
    }
}

impl Default for BigInteger<u16> {
    fn default() -> Self {
        BigInteger {
            data: vec![u16::default()],
            sign: false,
            basis: (1 << 16) as usize,
        }
    }
}

impl From<&str> for BigInteger<u16> {
    fn from(text: &str) -> Self {
        match text.parse::<u128>() {
            Ok(n) => {
                let mask = ((1 << 16) - 1) as u128;
                let num_words = 128 / 16;
                let mut number = n;
                let mut w =
                    BigInteger::<u16>::new(Vec::<u16>::new().as_slice(), false, 1 << 16, None);
                for _ in 0..num_words {
                    w.data.push((number & mask) as u16);
                    number = number >> 16;
                }
                w.strip_leading_zeros();
                w
            }
            Err(e) => BigInteger {
                data: vec![u16::default()],
                sign: false,
                basis: (1 << 16) as usize,
            },
        }
    }
}

impl BigInteger<u16> {
    fn squre_root(&self) -> BigInteger<u16> {
        if self.is_zero() {
            return self.clone();
        }
        assert!(self.is_positive());

        // initialization
        let n = self.data.len();
        let nr = if n % 2 == 0 { n / 2 } else { (n + 1) / 2 };
        let mut init_data = vec![0_u16; nr];
        init_data[nr - 1] = 1_u16;
        let mut t = BigInteger {
            data: init_data,
            sign: false,
            basis: self.basis,
        };
        let mut v = BigInteger {
            data: vec![0_u16],
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

    #[test]
    fn test_add() {
        let basis = (1 << 16) as usize;
        let (a_arr, b_arr, c_arr) = (vec![12_u16, 246_u16], vec![55_u16], vec![13_u16, 45_u16]);
        let (a, b, c) = (
            BigInteger::new(&a_arr, false, basis, Some(true)),
            BigInteger::new(&b_arr, false, basis, Some(true)),
            BigInteger::new(&c_arr, false, basis, Some(true)),
        );
        assert_eq!(&a + &b, c);
    }

    #[test]
    fn test_sub() {
        let basis = (1 << 16) as usize;
        let (a_arr, b_arr, c_arr) = (vec![12_u16, 46_u16], vec![55_u16], vec![11_u16, 247_u16]);
        let (a, b, c) = (
            BigInteger::new(&a_arr, false, basis, Some(true)),
            BigInteger::new(&b_arr, false, basis, Some(true)),
            BigInteger::new(&c_arr, false, basis, Some(true)),
        );
        assert_eq!(&a - &b, c);
    }

    #[test]
    fn test_mul() {
        let basis = (1 << 16) as usize;
        let (a_arr, b_arr, c_arr) = (vec![12_u16, 46_u16], vec![6_u16], vec![73_u16, 20_u16]);
        let (a, b, c) = (
            BigInteger::new(&a_arr, false, basis, Some(true)),
            BigInteger::new(&b_arr, false, basis, Some(true)),
            BigInteger::new(&c_arr, false, basis, Some(true)),
        );
        assert_eq!(&a * &b, c);
    }

    #[test]
    fn test_div() {
        let basis = (1 << 16) as usize;
        let (a_arr, b_arr, c_arr) = (
            vec![44_u16, 60_u16, 416_u16],
            vec![63_u16, 0_u16],
            vec![179_u16],
        );
        let (a, b, c) = (
            BigInteger::new(&a_arr, false, basis, Some(true)),
            BigInteger::new(&b_arr, false, basis, Some(true)),
            BigInteger::new(&c_arr, false, basis, Some(true)),
        );
        assert_eq!(&a / &b, c);
    }

    #[test]
    fn test_squre_root() {
        let basis = (1 << 16) as usize;
        let (a_arr, c_arr) = (vec![44_u16, 60_u16, 416_u16], vec![6_u16, 166_u16]);
        let (a, c) = (
            BigInteger::new(&a_arr, false, basis, Some(true)),
            BigInteger::new(&c_arr, false, basis, Some(true)),
        );
        assert_eq!(a.squre_root(), c);
    }

    #[test]
    fn test_from() {
        let a = BigInteger::from("26498041357");
        let b = BigInteger::new(
            vec![6_u16, 11112_u16, 1549_u16].as_slice(),
            false,
            1 << 16,
            Some(true),
        );
        assert_eq!(a, b);
    }
}
