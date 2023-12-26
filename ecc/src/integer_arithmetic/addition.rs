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
        let mut carrier = 0 as i32;
        for i in 0..max_n {
            let lft = if i <= nu - 1 { u[i] } else { 0 as i32 };
            let rht = if i <= nv - 1 { v[i] } else { 0 as i32 };
            let t = lft + rht + carrier;
            let remainder = if t > *b { t - b } else { t };
            carrier = t / b;
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
    #[test]
    fn test_addition() {}
}
