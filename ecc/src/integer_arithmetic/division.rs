///// Short division of positive multiple-precision integers
///// referenced by Algorithm 10.28 of "Handbook of Elliptic and Hyperelliptic Curve Cryptography"

use crate::integer_arithmetic::{
    multiplication::{MultiplePrecisionMultiplication, SinglePrecisionMultiplication},
    BigInteger,
};

pub trait SinglePrecision {
    fn divide_by_single_precision(self, denominator: i32) -> (BigInteger, i32);
}

pub trait MultiplePrecision {
    fn divide_by_multiple_precision(&self, denominator: &BigInteger) -> (BigInteger, BigInteger);
}

impl SinglePrecision for BigInteger {
    fn divide_by_single_precision(self, denominator: i32) -> (BigInteger, i32) {
        assert!(self.data.len() >= 1);

        let n = self.data.len();
        let b = self.basis;
        let mut r0 = 0 as i32;
        let u_ref = &self.data;
        let mut q = BigInteger {
            data: vec![],
            basis: b,
        };
        for i in (0..n).rev() {
            let t = &r0 * &b + u_ref[i];
            let quotient = t / denominator;
            r0 = t - quotient * denominator;
            if ((q.data.len() == 0) && (quotient > 0)) || (q.data.len() > 0) {
                q.data.push(quotient);
            }
        }
        (q, r0)
    }
}

impl MultiplePrecision for BigInteger {
    fn divide_by_multiple_precision(&self, denominator: &BigInteger) -> (BigInteger, BigInteger) {
        assert!(self.basis == denominator.basis);
        let b = self.basis;

        let (mut u, mut v) = (self.clone(), denominator.clone());
        let (nu, nv) = (u.data.len(), v.data.len());
        assert!(nu >= nv);

        u.data.push(0 as i32);
        let mut d = 1;

        // normalization
        while v.data[nv - 1] < b / 2 {
            v = v.multiply_single_precision(2 as i32);
            u = u.multiply_single_precision(2 as i32);
            d *= 2;
        }

        // proximate quotient
        for i in (0..nu).rev() {
            // two word proximate quotient
            let u_2w = BigInteger {
                data: u.data[nv + i - 1..nv + i + 1].to_vec(),
                basis: b,
            };
            let (q_2w, _) = u_2w.divide_by_single_precision(v.data[nv - 1]);
            let q_prox = match q_2w.data.len() {
                2 => b - 1,
                1 => q_2w.data[0],
                _ => 0,
            };
            // three word proximate quotient
            let v_2w = BigInteger {
                data: v.data[nv - 2..nv].to_vec(),
                basis: b,
            };
            let u_prox = v_2w.multiply_single_precision(q_prox);
            let u_3w = BigInteger {
                data: u.data[nv + i - 2..nv + i + 1].to_vec(),
                basis: b,
            };
        }
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_division() {
        let a = BigInteger {
            data: vec![8789],
            basis: 10,
        };
        let (q, r) = a.divide_by_single_precision(7 as i32);
        println!("q: {:?}, r: {}", q, r);
    }
}
