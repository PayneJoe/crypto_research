use ark_crypto_primitives::sponge::poseidon::{PoseidonDefaultConfig, PoseidonDefaultConfigEntry};
use ark_ff::fields::{Field, Fp64, MontBackend, MontConfig};
use num_bigint::*;

#[derive(MontConfig)]
#[modulus = "17"]
#[generator = "3"]
pub struct FqConfig;
pub type Fq = Fp64<MontBackend<FqConfig, 1>>;

// impl PoseidonDefaultConfig<4> for MontBackend<FqConfig, 1> {
//     const PARAMS_OPT_FOR_CONSTRAINTS: [PoseidonDefaultConfigEntry; 7] = [
//         PoseidonDefaultConfigEntry::new(2, 17, 8, 31, 0),
//         PoseidonDefaultConfigEntry::new(3, 5, 8, 56, 0),
//         PoseidonDefaultConfigEntry::new(4, 5, 8, 56, 0),
//         PoseidonDefaultConfigEntry::new(5, 5, 8, 57, 0),
//         PoseidonDefaultConfigEntry::new(6, 5, 8, 57, 0),
//         PoseidonDefaultConfigEntry::new(7, 5, 8, 57, 0),
//         PoseidonDefaultConfigEntry::new(8, 5, 8, 57, 0),
//     ];
//     const PARAMS_OPT_FOR_WEIGHTS: [PoseidonDefaultConfigEntry; 7] = [
//         PoseidonDefaultConfigEntry::new(2, 257, 8, 13, 0),
//         PoseidonDefaultConfigEntry::new(3, 257, 8, 13, 0),
//         PoseidonDefaultConfigEntry::new(4, 257, 8, 13, 0),
//         PoseidonDefaultConfigEntry::new(5, 257, 8, 13, 0),
//         PoseidonDefaultConfigEntry::new(6, 257, 8, 13, 0),
//         PoseidonDefaultConfigEntry::new(7, 257, 8, 13, 0),
//         PoseidonDefaultConfigEntry::new(8, 257, 8, 13, 0),
//     ];
// }

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
