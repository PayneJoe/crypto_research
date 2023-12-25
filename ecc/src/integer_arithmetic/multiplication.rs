use crate::integer_arithmetic::BigInteger;

pub trait SinglePrecisionMultiplication {
    fn multiply_single_precision(&self, multiplier: i32) -> BigInteger;
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
        let mut w = BigInteger {
            data: vec![0 as i32; nu + nv],
            basis: self.basis,
        };

        // initialization
        for i in 0..nv {
            w.data[i] = 0 as i32;
        }
        // cross multiplication
        for i in 0..nv {
            let carrier = 0 as i32;
            if v[i] == 0 {
                w.data[nu + i] = 0 as i32;
            } else {
                for j in 0..nu {
                    let t = w.data[i + j] + v[i] * u[j] + carrier;
                    let quotient = t / b;
                    w.data[i + j] = t - quotient * b;
                }
                w.data[nu + i] = carrier;
            }
        }
        w
    }
}

impl SinglePrecisionMultiplication for BigInteger {
    fn multiply_single_precision(&self, multiplier: i32) -> BigInteger {
        assert!(multiplier < self.basis);
        let u = &self.data;
        let b = self.basis;

        let n = self.data.len();
        let mut w = BigInteger {
            data: vec![0 as i32; n + 1],
            basis: self.basis,
        };
        let mut carrier = 0 as i32;
        for i in 0..n {
            let t = w.data[i] + u[i] * multiplier + carrier;
            carrier = t / b;
            w.data[i] = t - carrier * b;
        }
        w.data[n] = carrier;
        w
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_single_precision() {}
}
