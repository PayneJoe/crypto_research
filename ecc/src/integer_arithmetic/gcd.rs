use crate::integer_arithmetic::BigInteger;

/// gcd for single precision of positive integer
/// referenced by Algorithm 10.42 of "Handbook of Elliptic and Hyperelliptic Curve Cryptography"
pub trait SinglePrecisionGCD {
    fn euclid_extended_gcd(x: i32, N: i32) -> (i32, i32, i32) {
        assert!(x < N);
        let (mut A, mut B) = (N, x);
        let (mut Ua, mut Ub) = (0, 1);
        let (mut Va, mut Vb) = (1, 0);
        while B != 0 {
            let q = A / B;
            (A, B) = (B, A - q * B);
            (Ua, Ub) = (Ub, Ua - q * Ub);
            (Va, Vb) = (Vb, Va - q * Vb);
        }
        let (d, u, v) = (A, Ua, Va);
        (u, v, d)
    }
}

/// gcd for multiple precision of positive integers
/// referenced by Algorithm 10.45 and 10.46 of "Handbook of Elliptic and Hyperelliptic Curve Cryptography"
pub trait MultiplePrecisionGCD {
    fn lehmer_extended_gcd(x: BigInteger, N: BigInteger) -> (BigInteger, BigInteger, BigInteger) {
        todo!()
    }
}

impl SinglePrecisionGCD for BigInteger {}
impl MultiplePrecisionGCD for BigInteger {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_not_coprime_gcd() {
        let (u, v, d) = BigInteger::euclid_extended_gcd(54, 888);
        assert_eq!(d, 6);
    }

    #[test]
    fn test_coprime_gcd() {
        let (u, v, d) = BigInteger::euclid_extended_gcd(45, 127);
        assert_eq!(d, 1);
    }
}
