use crate::integer_arithmetic::{
    addition::Addition,
    division::{MultiplePrecision, SinglePrecision},
    BigInteger,
};

pub trait SqureRoot {
    fn squre_root_newton_approximate(&self) -> Self;
}

/// squre root with Newton approximation method
/// referenced by Algorithm 10.55 from "Handbook of Elliptic and Hyperelliptic Curve Cryptography"
impl SqureRoot for BigInteger {
    fn squre_root_newton_approximate(&self) -> Self {
        if self.is_zero() {
            return BigInteger {
                data: vec![0_u8],
                basis: self.basis,
            };
        }

        // initialization
        let n = self.data.len();
        let nr = if n % 2 == 0 { n / 2 } else { (n + 1) / 2 };
        let mut init_data = vec![0_u8; nr];
        init_data[nr - 1] = 1_u8;
        let mut t = BigInteger {
            data: init_data,
            basis: self.basis,
        };
        let mut v = BigInteger {
            data: vec![0_u8],
            basis: self.basis,
        };

        // iterate with newton method
        loop {
            v = t;
            let (delta, _) = self.divide_by_multiple_precision(&v);
            (t, _) = v.add(&delta).divide_by_single_precision(2);
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
    fn test_newton_squre_root() {
        let mut input_data = vec![1, 6, 1];
        let mut actual_data = vec![1, 2];
        input_data.reverse();
        actual_data.reverse();

        let u = BigInteger {
            data: input_data,
            basis: 10,
        };
        let result = u.squre_root_newton_approximate();
        assert_eq!(
            result,
            BigInteger {
                data: actual_data,
                basis: 10
            }
        );
    }
}
