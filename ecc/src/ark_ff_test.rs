use ark_ff::fields::{Field, Fp64, MontBackend, MontConfig};
use num_bigint::*;

#[derive(MontConfig)]
#[modulus = "17"]
#[generator = "3"]
pub struct FqConfig;
pub type Fq = Fp64<MontBackend<FqConfig, 1>>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_derive() {
        let a = Fq::from(9);
        let b = Fq::from(10);

        let a_bint: num_bigint::BigUint = a.into();
        assert_eq!(a_bint.to_bytes_le()[0] & 1, 1);
        let m = num_bigint::BigUint::from(4u8);
        assert_eq!(m.to_bytes_le()[0] & 1, 0);

        assert_eq!(a, Fq::from(26)); // 26 =  9 mod 17
        assert_eq!(a - b, Fq::from(16)); // -1 = 16 mod 17
        assert_eq!(a + b, Fq::from(2)); // 19 =  2 mod 17
        assert_eq!(a * b, Fq::from(5)); // 90 =  5 mod 17
        assert_eq!(a.square(), Fq::from(13)); // 81 = 13 mod 17
        assert_eq!(b.double(), Fq::from(3)); // 20 =  3 mod 17
        assert_eq!(a / b, a * b.inverse().unwrap());

        assert_eq!(1 == 1, true);
    }
}
