use crate::integer_arithmetic::{
    addition::Addition, multiplication::SinglePrecisionMultiplication, substraction::Substraction,
    BigInteger,
};

pub trait SinglePrecision {
    fn divide_by_single_precision(&self, denominator: u8) -> (BigInteger, u8);
}

pub trait MultiplePrecision {
    fn divide_by_multiple_precision(&self, denominator: &BigInteger) -> (BigInteger, BigInteger);
}

///// Short division of positive multiple-precision integers
///// referenced by Algorithm 10.28 of "Handbook of Elliptic and Hyperelliptic Curve Cryptography"
///
impl SinglePrecision for BigInteger {
    fn divide_by_single_precision(&self, denominator: u8) -> (BigInteger, u8) {
        assert!(self.data.len() >= 1);
        assert!(denominator != 0);

        let n = self.data.len();
        let b = self.basis;
        let mut carrier = 0 as u8;
        let u = &self.data;
        let mut q = BigInteger {
            data: vec![],
            basis: b,
        };
        for i in (0..n).rev() {
            let t = &carrier * &b + u[i];
            let quotient = t / denominator;
            carrier = t - quotient * denominator;
            // the highest must be non-zero
            if ((q.data.len() == 0) && (quotient > 0)) || (q.data.len() > 0) {
                q.data.push(quotient);
            }
        }
        if q.data.len() == 0 {
            q.data.push(0);
        }
        q.data.reverse();
        (q, carrier)
    }
}

//// Long division of positive multiple-precision integers
/// referenced by Algorithm 10.30 of "Handbook of Elliptic and Hyperelliptic Curve Cryptography"
///
impl MultiplePrecision for BigInteger {
    fn divide_by_multiple_precision(&self, denominator: &BigInteger) -> (BigInteger, BigInteger) {
        assert!(self.basis == denominator.basis);
        let b = self.basis;

        let (mut u, mut v) = (self.clone(), denominator.clone());
        let (nu, nv) = (u.data.len(), v.data.len());
        assert!(nu >= nv);

        // reset highest digit
        u.data.push(0 as u8);

        // normalization
        let mut d = 1;
        println!(
            "**** Before normalization: u = {:?}, v = {:?}, d = {}",
            u.data, v.data, d
        );
        while v.data[nv - 1] < b / 2 {
            v = v.multiply_single_precision(2 as u8);
            u = u.multiply_single_precision(2 as u8);
            d = d * 2;
        }
        println!(
            "**** After normalization: u = {:?}, v = {:?}, d = {}",
            u.data, v.data, d
        );

        // proximate quotient
        let mut q = BigInteger {
            data: vec![],
            basis: b,
        };
        for i in (0..(nu - nv + 1)).rev() {
            // two word proximate quotient
            let u_2w = BigInteger {
                data: u.data[nv + i - 1..nv + i + 1].to_vec(),
                basis: b,
            };
            let (q_2w, _) = &u_2w.divide_by_single_precision(v.data[nv - 1]);
            let mut q_prox = match q_2w.data.len() {
                2 => b - 1,
                1 => q_2w.data[0],
                _ => 0,
            };
            println!(
                "**** After 2-word approximation: q_prox = {}, q_2w = {:?}, u_2w = {:?}, v[n - 1] = {}",
                q_prox,
                q_2w.data,
                &u_2w.data,
                v.data[nv - 1]
            );

            // three word proximate quotient
            let v_2w = BigInteger {
                data: v.data[nv - 2..nv].to_vec(),
                basis: b,
            };
            let mut u_3w_prox = v_2w.multiply_single_precision(q_prox);
            let u_3w = BigInteger {
                data: u.data[nv + i - 2..nv + i + 1].to_vec(),
                basis: b,
            };
            while u_3w_prox.greater(&u_3w) {
                q_prox = q_prox - 1;
                u_3w_prox = v_2w.multiply_single_precision(q_prox);
            }
            println!(
                "**** After 3-words approximation: q_prox = {}, u_3w = {:?}, u_3w_prox = {:?}",
                q_prox, u_3w.data, u_3w_prox.data
            );

            // n_plus_one word proximate quotient
            let u_n1w_prox = v.multiply_single_precision(q_prox);
            let u_n1w = BigInteger {
                data: u.data[i..(i + nv + 1)].to_vec(),
                basis: b,
            };
            if u_n1w.greater(&u_n1w_prox) {
                let diff = u_n1w.substract(&u_n1w_prox);
                u.data[i..(i + nv + 1)].copy_from_slice(diff.data.as_slice());
            } else {
                q_prox = q_prox - 1;
                let diff = u_n1w.add(&v).substract(&u_n1w_prox);
                u.data[i..(i + nv + 1)].copy_from_slice(diff.data.as_slice());
            }
            println!(
                "***** After last proximation: q_prox = {}, u[i..i+n-1] = {:?}",
                q_prox,
                u.data[i..(i + nv + 1)].to_vec(),
            );

            // the highest must be non-zero
            if ((q.data.len() == 0) && (q_prox > 0)) || (q.data.len() > 0) {
                q.data.push(q_prox);
            }
            println!("------------- \n");
        }
        q.data.reverse();

        let r = u.divide_by_single_precision(d).0;
        (q, r)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_division() {
        let mut a_arr = vec![8, 7, 8, 9];
        let b = 7;
        let mut c_arr = vec![1, 2, 5, 5];
        a_arr.reverse();
        c_arr.reverse();

        let a = BigInteger {
            data: a_arr,
            basis: 10,
        };
        let c = BigInteger {
            data: c_arr,
            basis: 10,
        };
        assert_eq!(a.divide_by_single_precision(b).0, c);
    }

    #[test]
    fn test_multiple_division() {
        let mut a_arr = vec![1, 1, 5, 9, 2, 3];
        let mut b_arr = vec![3, 4, 4];
        let mut c_arr = vec![3, 3, 6];
        a_arr.reverse();
        b_arr.reverse();
        c_arr.reverse();

        let a = BigInteger {
            data: a_arr,
            basis: 10,
        };
        let b = BigInteger {
            data: b_arr,
            basis: 10,
        };
        let c = BigInteger {
            data: c_arr,
            basis: 10,
        };
        assert_eq!(a.divide_by_multiple_precision(&b).0, c);
    }

    #[test]
    fn test_pow2_basis() {
        let basis: u8 = 1 << 6;
        let mut a_arr = vec![11, 15, 12];
        let mut b_arr = vec![15, 48];
        let mut c_arr = vec![45];
        let mut d_arr = vec![10, 28];
        a_arr.reverse();
        b_arr.reverse();
        c_arr.reverse();
        d_arr.reverse();

        let a = BigInteger {
            data: a_arr,
            basis: basis,
        };
        let b = BigInteger {
            data: b_arr,
            basis: basis,
        };
        let c = BigInteger {
            data: c_arr,
            basis: basis,
        };
        let d = BigInteger {
            data: d_arr,
            basis: basis,
        };
        let result = a.divide_by_multiple_precision(&b);
        assert_eq!(result.0, c);
        assert_eq!(result.1, d);
    }
}
