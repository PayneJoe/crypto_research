use crate::elliptic_curve_arithmetic::pairing_friendly::bls12::{g1::G1, g2::G2, gt::Gt};
use crate::finite_field_arithmetic::bigint::BigInt;

const MAX_SCALAR_LIMBS: usize = 8;

pub trait WeilPairing {
    fn miller_loop(&self, Q: Self, r: BigInt<MAX_SCALAR_LIMBS>) -> Gt;
}
