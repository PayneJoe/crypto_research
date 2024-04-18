pub mod g1;
pub mod g12;
pub mod g2;
pub mod gt;

const BASE_NUM_LIMBS: usize = 6;
const SCALAR_NUM_LIMBS: usize = 4;

use crate::elliptic_curve_arithmetic::pairing_friendly::bls12::{g1::G1, g12::G12, g2::G2};
use crate::finite_field_arithmetic::pairing_friendly::bls12::{
    fq::Fq, fq12::Fq12, fq2::Fq2, fr::Fr,
};
use crate::pairings::weil::WeilPairing;
pub struct bls12_381;

impl WeilPairing for bls12_381 {
    type G1BaseField = Fq<BASE_NUM_LIMBS>;
    type ScalarField = Fr<SCALAR_NUM_LIMBS>;
    type Gt = Fq12;
    type G1Curve = G1;
    type G2Curve = G12;
}

mod tests {
    #[test]
    fn test_weil_pairing() {
        todo!()
    }
}
