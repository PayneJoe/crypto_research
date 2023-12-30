use std::ops::{Add, Div, Mul, Shl, Shr, Sub};

use crate::finite_field_arithmetic::BigInteger;

fn div_internal(lft: &BigInteger<u8>, rht: &BigInteger<u8>) -> (BigInteger<u8>, BigInteger<u8>) {
    assert!(lft.basis == rht.basis);
    let b = lft.basis;

    let (mut u, mut v) = (lft.clone(), rht.clone());
    let (nu, nv) = (u.data.len(), v.data.len());
    assert!(nu >= nv);

    // reset highest digit
    u.data.push(0 as u8);

    // normalization
    let mut d = 1_u8;
    println!(
        "**** Before normalization: u = {:?}, v = {:?}, d = {}",
        u.data, v.data, d
    );
    while v.data[nv - 1] < (b / 2) as u8 {
        v = &v << (2 as usize);
        u = &u << (2 as usize);
        d = d * 2;
    }
    println!(
        "**** After normalization: u = {:?}, v = {:?}, d = {}",
        u.data, v.data, d
    );

    // proximate quotient
    let mut q = BigInteger {
        data: vec![0_u8; nu - nv + 1],
        sign: false,
        basis: b,
    };
    for i in (0..(nu - nv + 1)).rev() {
        // two word proximate quotient
        let u_2w = u16::from(u.data[nv + i] << 8 + u.data[nv + i - 1]);
        let v_1w = u16::from(v.data[nv - 1]);
        let mut q_prox = (std::cmp::min(u_2w / v_1w, (b - 1) as u16)) as u8;

        println!(
            "**** After 2-word approximation: q_prox = {}, q_2w = {}, u_2w = {}, v[n - 1] = {}",
            q_prox,
            u_2w / v_1w,
            u_2w,
            v_1w
        );

        // three word proximate quotient
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

        println!(
            "**** After 3-words approximation: q_prox = {}, u_3w = {:?}, u_3w_prox = {:?}",
            q_prox, u_3w.data, u_3w_prox.data
        );

        // n_plus_one word proximate quotient
        let u_nplus1w_prox = &v * (q_prox as usize);
        let u_nplus1w = BigInteger {
            data: u.data[i..(i + nv + 1)].to_vec(),
            sign: u.sign,
            basis: b,
        };
        let diff = &u_nplus1w - &u_nplus1w_prox;
        if diff.is_positive() {
            u.data[i..(i + nv + 1)].copy_from_slice(diff.data.as_slice());
        } else {
            q_prox = q_prox - 1;
            u.data[i..(i + nv + 1)].copy_from_slice((&diff + &v).data.as_slice());
        }
        println!(
            "***** After last proximation: q_prox = {}, u[i..i+n-1] = {:?}",
            q_prox,
            u.data[i..(i + nv + 1)].to_vec(),
        );

        q.data[i] = q_prox;
        println!("------------- \n");
    }
    q.strip_leading_zeros();

    let r = &u / (d as usize);
    (q, r)
}

fn div_single_precision_internal(lft: &BigInteger<u8>, rht: usize) -> BigInteger<u8> {
    assert!(rht != 0);
    let (b, u, nu, v) = (&lft.basis, &lft.data, lft.size(), rht as u8);
    let mut w = BigInteger {
        data: vec![0 as u8; nu - 1 + 1],
        sign: lft.sign,
        basis: *b,
    };

    let mut carrier = 0_u8;
    for i in (0..nu).rev() {
        let (carrier_word, unit_word, nominator_word, denominator_word) =
            (u16::from(carrier), *b as u16, u16::from(u[i]), u16::from(v));
        let t = carrier_word * unit_word + nominator_word;
        let quotient = t / denominator_word;
        carrier = (t - quotient * denominator_word) as u8;
        w.data[i] = quotient as u8;
    }
    w.strip_leading_zeros();
    w
}

fn shl_internal(lft: &BigInteger<u8>, rht: usize) -> BigInteger<u8> {
    let n = lft.size();
    let mut w = BigInteger {
        data: vec![0_u8; n + 1],
        sign: false,
        basis: lft.basis,
    };
    let mut carrier = 0_u8;
    // from the lowest word to the highest word
    for i in 0..n {
        let remainder = (lft.data[i] << rht) | carrier;
        carrier = lft.data[i] >> (8 - rht);
        w.data[i] = remainder;
    }
    if (carrier > 0) {
        w.data[n] = carrier;
    }
    w.strip_leading_zeros();
    w
}

fn shr_internal(lft: &BigInteger<u8>, rht: usize) -> BigInteger<u8> {
    let n = lft.size();
    let mut w = BigInteger {
        data: vec![0_u8; n + 1],
        sign: false,
        basis: lft.basis,
    };
    let mut carrier = 0_u8;
    // from the hightest word to the lowest word
    for i in (0..n).rev() {
        let remainder = (lft.data[i] >> rht) | carrier;
        carrier = lft.data[i] << (8 - rht);
        w.data[i] = remainder;
    }
    w.strip_leading_zeros();
    w
}

fn mul_single_precision_internal(lft: &BigInteger<u8>, rht: usize) -> BigInteger<u8> {
    let (b, u, nu, v) = (&lft.basis, &lft.data, lft.size(), rht as u8);
    let mut w = BigInteger {
        data: vec![0 as u8; nu + 1],
        sign: lft.sign,
        basis: *b,
    };
    if v != 0 {
        let mut carrier = 0_u8;
        let mut remainder = 0_u8;
        for i in 0..nu {
            // double precision for single precision multiplication with carriers
            let (w_word, lft_word, rht_word, carrier_word) = (
                u16::from(w.data[i]),
                u16::from(u[i]),
                u16::from(v),
                u16::from(carrier),
            );
            let t = w_word + lft_word * rht_word + carrier_word;
            (carrier, remainder) = ((t >> 8) as u8, (t & ((1 << 8) - 1)) as u8);

            w.data[i] = remainder;
        }
        if carrier > 0 {
            w.data[nu] = carrier;
        }
    }
    w.strip_leading_zeros();
    w
}

fn mul_internal(lft: &BigInteger<u8>, rht: &BigInteger<u8>) -> BigInteger<u8> {
    assert!(lft.basis == rht.basis);
    let b = &lft.basis;
    let (u, v) = (&lft.data, &rht.data);
    let (nu, nv) = (u.len(), v.len());

    let mut w = BigInteger {
        data: vec![0 as u8; nu + nv],
        sign: false,
        basis: *b,
    };

    // initialization
    for i in 0..nv {
        w.data[i] = 0 as u8;
    }
    // cross multiplication
    for i in 0..nv {
        let mut carrier = 0 as u8;
        let mut remainder = 0 as u8;
        if v[i] != 0_u8 {
            for j in 0..nu {
                // double precision for single precision multiplication with carriers
                let (w_word, lft_word, rht_word, carrier_word) = (
                    u16::from(w.data[i + j]),
                    u16::from(u[j]),
                    u16::from(v[i]),
                    u16::from(carrier),
                );
                let t = w_word + lft_word * rht_word + carrier_word;
                (carrier, remainder) = ((t >> 8) as u8, (t & ((1 << 8) - 1)) as u8);

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

fn add_internal(lft: &BigInteger<u8>, rht: &BigInteger<u8>) -> BigInteger<u8> {
    assert!(lft.basis == rht.basis);
    let b = &lft.basis;
    let (u, v) = (&lft.data, &rht.data);
    let (nu, nv) = (u.len(), v.len());

    let max_n = std::cmp::max(nu, nv);
    let mut w = BigInteger {
        data: vec![0_u8; max_n + 1],
        sign: false,
        basis: *b,
    };
    let mut carrier = 0 as u8;
    let mut remainder = 0 as u8;
    // from right to left
    for i in 0..max_n {
        let lft = if i <= nu - 1 { u[i] } else { 0_u8 };
        let rht = if i <= nv - 1 { v[i] } else { 0_u8 };

        // double precision for single precision addition with carriers
        let result = u16::from(lft) + u16::from(rht) + u16::from(carrier);
        (carrier, remainder) = ((result >> 8) as u8, (result & ((1 << 8) - 1)) as u8);

        w.data[i] = remainder;
    }
    // the highest digit must be non-zero
    if carrier > 0 {
        w.data[max_n] = carrier;
    }
    w.strip_leading_zeros();
    w
}

fn sub_internal(lft: &BigInteger<u8>, rht: &BigInteger<u8>) -> BigInteger<u8> {
    assert!(lft.basis == rht.basis);
    let b = &lft.basis;
    let (u, v) = (&lft.data, &rht.data);
    let (nu, nv) = (u.len(), v.len());

    let max_n = std::cmp::max(nu, nv);
    let mut w = BigInteger {
        data: vec![0_u8; max_n],
        sign: false,
        basis: *b,
    };
    let mut carrier = 0 as u8;
    let mut remainder = 0 as u8;
    // from right to left
    for i in 0..max_n {
        let lft = if i <= nu - 1 { u[i] } else { 0_u8 };
        let rht = if i <= nv - 1 { v[i] } else { 0_u8 };

        // double precision for single precision substraction with carriers
        let result = u16::from(lft)
            .wrapping_sub(u16::from(rht))
            .wrapping_sub(u16::from(carrier));
        (carrier, remainder) = (
            ((result >> 8) as u8).wrapping_neg(),
            (result & ((1 << 8) - 1)) as u8,
        );

        w.data[i] = remainder;
    }
    w.strip_leading_zeros();
    // we need to reverse the digits where the result has carrier
    if carrier > 0 {
        for i in 0..w.size() {
            w.data[i] = (*b as u8) - w.data[i];
        }
        w.sign = true;
    }
    w
}

//// &BigInteger + &BigInteger
impl<'a, 'b> Add<&'b BigInteger<u8>> for &'a BigInteger<u8> {
    type Output = BigInteger<u8>;

    fn add(self, other: &'b BigInteger<u8>) -> BigInteger<u8> {
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

impl<'a, 'b> Sub<&'b BigInteger<u8>> for &'a BigInteger<u8> {
    type Output = BigInteger<u8>;

    fn sub(self, other: &'b BigInteger<u8>) -> BigInteger<u8> {
        let mut other_new = other.clone();
        other_new.sign = !other_new.sign;
        self.add(&other_new)
    }
}

impl<'a, 'b> Mul<&'b BigInteger<u8>> for &'a BigInteger<u8> {
    type Output = BigInteger<u8>;

    fn mul(self, other: &'b BigInteger<u8>) -> BigInteger<u8> {
        let mut result = mul_internal(self, other);
        result.sign = self.sign & other.sign;
        result
    }
}

impl Mul<usize> for &BigInteger<u8> {
    type Output = BigInteger<u8>;

    fn mul(self, other: usize) -> BigInteger<u8> {
        assert!(other < 1 << 8);
        let mut result = mul_single_precision_internal(self, other);
        result.sign = self.sign;
        result
    }
}

impl Shl<usize> for &BigInteger<u8> {
    type Output = BigInteger<u8>;

    fn shl(self, rhs: usize) -> Self::Output {
        assert!(rhs < 8);
        let mut result = shl_internal(self, rhs);
        result.sign = self.sign;
        result
    }
}

impl Shr<usize> for &BigInteger<u8> {
    type Output = BigInteger<u8>;

    fn shr(self, rhs: usize) -> Self::Output {
        assert!(rhs < 8);
        let mut result = shr_internal(self, rhs);
        result.sign = self.sign;
        result
    }
}

impl Div<usize> for &BigInteger<u8> {
    type Output = BigInteger<u8>;

    fn div(self, other: usize) -> BigInteger<u8> {
        assert!(other < 1 << 8);
        let mut result = div_single_precision_internal(self, other);
        result.sign = self.sign;
        result
    }
}

impl<'a, 'b> Div<&'b BigInteger<u8>> for &'a BigInteger<u8> {
    type Output = BigInteger<u8>;

    fn div(self, other: &'b BigInteger<u8>) -> BigInteger<u8> {
        let (mut quotient, _) = div_internal(self, other);
        quotient.sign = self.sign & other.sign;
        quotient
    }
}

impl BigInteger<u8> {
    fn squre_root(&self) -> BigInteger<u8> {
        if self.is_zero() {
            return self.clone();
        }
        assert!(self.is_positive());

        // initialization
        let n = self.data.len();
        let nr = if n % 2 == 0 { n / 2 } else { (n + 1) / 2 };
        let mut init_data = vec![0_u8; nr];
        init_data[nr - 1] = 1_u8;
        let mut t = BigInteger {
            data: init_data,
            sign: false,
            basis: self.basis,
        };
        let mut v = BigInteger {
            data: vec![0_u8],
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
    fn test_add_1() {
        let basis = (1 << 8) as usize;
        let mut a_arr: Vec<u8> = vec![12, 246];
        let mut b_arr: Vec<u8> = vec![55];
        let mut c_arr: Vec<u8> = vec![13, 45];
        a_arr.reverse();
        b_arr.reverse();
        c_arr.reverse();
        let (a, b, c) = (
            BigInteger {
                data: a_arr.to_vec(),
                sign: false,
                basis: basis,
            },
            BigInteger {
                data: b_arr.to_vec(),
                sign: false,
                basis: basis,
            },
            BigInteger {
                data: c_arr.to_vec(),
                sign: false,
                basis: basis,
            },
        );
        assert_eq!(&a + &b, c);
    }

    #[test]
    fn test_sub_1() {
        let basis = (1 << 8) as usize;
        let mut a_arr: Vec<u8> = vec![12, 46];
        let mut b_arr: Vec<u8> = vec![55];
        let mut c_arr: Vec<u8> = vec![11, 247];
        a_arr.reverse();
        b_arr.reverse();
        c_arr.reverse();
        let (a, b, c) = (
            BigInteger {
                data: a_arr.to_vec(),
                sign: false,
                basis: basis,
            },
            BigInteger {
                data: b_arr.to_vec(),
                sign: false,
                basis: basis,
            },
            BigInteger {
                data: c_arr.to_vec(),
                sign: false,
                basis: basis,
            },
        );
        assert_eq!(&a - &b, c);
    }

    #[test]
    fn test_mul_1() {
        let basis = (1 << 8) as usize;
        let mut a_arr: Vec<u8> = vec![12, 46];
        let mut b_arr: Vec<u8> = vec![6];
        let mut c_arr: Vec<u8> = vec![73, 20];
        a_arr.reverse();
        b_arr.reverse();
        c_arr.reverse();
        let (a, b, c) = (
            BigInteger {
                data: a_arr.to_vec(),
                sign: false,
                basis: basis,
            },
            BigInteger {
                data: b_arr.to_vec(),
                sign: false,
                basis: basis,
            },
            BigInteger {
                data: c_arr.to_vec(),
                sign: false,
                basis: basis,
            },
        );
        assert_eq!(&a * &b, c);
    }
}
