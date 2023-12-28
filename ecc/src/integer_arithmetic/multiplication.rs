use crate::integer_arithmetic::BigInteger;

pub trait SinglePrecisionMultiplication {
    fn multiply_single_precision(&self, multiplier: u8) -> BigInteger;
}

pub trait MultiplePrecisionMultiplication {
    fn multiply_multiple_precision(&self, multiplier: &BigInteger) -> BigInteger;
}

impl MultiplePrecisionMultiplication for BigInteger {
    fn multiply_multiple_precision(&self, multiplier: &BigInteger) -> BigInteger {
        assert!(self.basis == multiplier.basis);
        let b = self.basis;

        let (u, v) = (&self.data, &multiplier.data);
        let (nu, nv) = (u.len(), v.len());
        assert!((nu != 0) && (nv != 0));
        if ((nu == 1) && (u[0] == 0)) || ((nv == 1) && (v[0] == 0)) {
            return BigInteger {
                data: vec![0],
                basis: b,
            };
        }
        if nu > 1 {
            assert!(u[nu - 1] != 0);
        }
        if nv > 1 {
            assert!(v[nv - 1] != 0);
        }

        let mut w = BigInteger {
            data: vec![0 as u8; nu + nv],
            basis: self.basis,
        };

        // initialization
        for i in 0..nv {
            w.data[i] = 0 as u8;
        }
        // cross multiplication
        for i in 0..nv {
            let mut carrier = 0 as u8;
            if v[i] != 0 {
                for j in 0..nu {
                    let t = w.data[i + j] + v[i] * u[j] + carrier;
                    carrier = t / b;
                    w.data[i + j] = t - carrier * b;
                }
                if carrier > 0 {
                    w.data[nu + i] = carrier;
                }
            }
        }
        // the highest must be non-zero
        if w.data[nu + nv - 1] == 0 {
            w.data = w.data[..nu + nv - 1].to_vec();
        }

        w
    }
}

impl SinglePrecisionMultiplication for BigInteger {
    fn multiply_single_precision(&self, multiplier: u8) -> BigInteger {
        assert!(multiplier < self.basis);
        let u = &self.data;
        let b = self.basis;

        if multiplier == 0 {
            return BigInteger {
                data: vec![0],
                basis: b,
            };
        }

        let n = self.data.len();
        let mut w = BigInteger {
            // data: vec![0 as u8; n + 1],
            data: vec![],
            basis: self.basis,
        };
        let mut carrier = 0 as u8;
        for i in 0..n {
            let t = /*w.data[i] + */ u[i] * multiplier + carrier;
            carrier = t / b;
            let remainder = t - carrier * b;
            w.data.push(remainder);
        }
        // the highest must be non-zero
        if carrier > 0 {
            w.data.push(carrier);
        }
        w
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multiple_precision_multiplication() {
        let mut a_arr = vec![2, 7, 3];
        let mut b_arr = vec![1, 2];
        let mut c_arr = vec![3, 2, 7, 6];
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
        assert_eq!(a.multiply_multiple_precision(&b), c);
    }

    #[test]
    fn test_single_precision_multiplication() {
        let mut a_arr = vec![2, 7, 3];
        let b = 6;
        let mut c_arr = vec![1, 6, 3, 8];
        let mut d_arr = vec![8, 1, 9];
        a_arr.reverse();
        c_arr.reverse();
        d_arr.reverse();
        let a = BigInteger {
            data: a_arr,
            basis: 10,
        };
        let c = BigInteger {
            data: c_arr,
            basis: 10,
        };
        let d = BigInteger {
            data: d_arr,
            basis: 10,
        };
        assert_eq!(a.multiply_single_precision(b), c);
        assert_eq!(a.multiply_single_precision(3), d);
    }
}
