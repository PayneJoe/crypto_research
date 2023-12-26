use crate::integer_arithmetic::BigInteger;

pub trait Substraction {
    fn substract(&self, other: &BigInteger) -> BigInteger;
    fn greater(&self, other: &BigInteger) -> bool;
}

impl Substraction for BigInteger {
    fn greater(&self, other: &BigInteger) -> bool {
        assert!(self.basis == other.basis);
        let b = &self.basis;
        let (u, v) = (&self.data, &other.data);
        let (nu, nv) = (u.len(), v.len());
        assert!((self.data[nu - 1] != 0) && (other.data[nv - 1] != 0));
        if nu < nv {
            return false;
        } else if nu > nv {
            return true;
        }

        let mut i = nu - 1;
        while (i >= 0) && (u[i] == v[i]) {
            i -= 1
        }
        if (i >= 0) && (u[i] > v[i]) {
            true
        } else {
            false
        }
    }

    fn substract(&self, other: &BigInteger) -> BigInteger {
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
            let lft = if i <= nu - 1 {
                u[i] - carrier
            } else {
                -carrier as i32
            };
            let rht = if i <= nv - 1 { v[i] } else { 0 as i32 };
            let remainder = if lft > rht { lft - rht } else { lft + b - rht };
            carrier = if lft > rht { 0 } else { 1 };
            w.data.push(remainder);
        }
        w
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_substraction() {}
}
