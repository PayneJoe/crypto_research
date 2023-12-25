use core::fmt::Debug;

#[derive(Debug, Clone, PartialEq)]
pub struct BigInteger {
    data: Vec<i32>,
    basis: i32,
}

pub trait SinglePrecision {
    fn divide_by(self, denominator: i32) -> (BigInteger, i32);
}

impl SinglePrecision for BigInteger {
    fn divide_by(self, denominator: i32) -> (BigInteger, i32) {
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
            q.data.push(quotient);
        }
        (q, r0)
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
        let (q, r) = a.divide_by(7 as i32);
        println!("q: {:?}, r: {}", q, r);
    }
}
