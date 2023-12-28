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
        assert!((nu != 0) && (nv != 0));
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
        let mut carrier = 0 as u8;
        for i in 0..max_n {
            let lft = if i <= nu - 1 { u[i] } else { 0 };
            let rht = if i <= nv - 1 {
                v[i] + carrier
            } else {
                carrier as u8
            };
            let remainder = if lft >= rht { lft - rht } else { lft + b - rht };
            carrier = if lft >= rht { 0 } else { 1 };
            w.data.push(remainder);
        }
        assert!(carrier == 0);
        w
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_substraction() {
        let mut a_arr = vec![5, 8, 6, 4, 6];
        let mut b_arr = vec![2, 1, 6, 5];
        let mut c_arr = vec![5, 6, 4, 8, 1];
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
        let result = a.substract(&b);
        assert_eq!(c, result);
    }

    #[test]
    fn test_greater() {
        let mut a_arr = vec![5, 8, 6, 4, 6];
        let mut b_arr = vec![2, 1, 6, 5];
        a_arr.reverse();
        b_arr.reverse();
        let a = BigInteger {
            data: a_arr,
            basis: 10,
        };
        let b = BigInteger {
            data: b_arr,
            basis: 10,
        };
        let result = a.greater(&b);
        assert_eq!(result, true);
        assert_eq!(b.greater(&a), false);
    }
}
