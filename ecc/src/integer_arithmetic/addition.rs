use crate::integer_arithmetic::BigInteger;

pub trait Addition {
    fn add(&self, other: &BigInteger) -> BigInteger;
}

impl Addition for BigInteger {
    fn add(&self, other: &BigInteger) -> BigInteger {
        assert!(self.basis == other.basis);
        let b = &self.basis;
        let (u, v) = (&self.data, &other.data);
        let (nu, nv) = (u.len(), v.len());

        let mut w = BigInteger {
            data: vec![],
            basis: *b,
        };
        let max_n = std::cmp::max(nu, nv);
        let mut carrier = 0 as u8;
        // from right to left
        for i in 0..max_n {
            let lft = if i <= nu - 1 { u[i] } else { 0 as u8 };
            let rht = if i <= nv - 1 { v[i] } else { 0 as u8 };
            let t = lft + rht + carrier;
            let remainder = if t >= *b { t - b } else { t };
            carrier = if t >= *b { 1 } else { 0 };
            w.data.push(remainder);
        }
        // the highest digit must be non-zero
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
    fn test_addition() {
        let mut a_arr = vec![9, 6, 3, 5];
        let mut b_arr = vec![8, 2, 7];
        let mut c_arr = vec![1, 0, 4, 6, 2];
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
        assert_eq!(a.add(&b), c);
    }
}
