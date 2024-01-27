use crate::ipa::{Commitment, IPAProof};
use crate::polynomial::sparse::SparsePolynomial;

use ecc::elliptic_curve_arithmetic::weierstrass_model::{pallas::curve::Pallas, AffinePoint};
use ecc::finite_field_arithmetic::bigint::BigInt;
use ecc::finite_field_arithmetic::pallas::{fq::Fq, fr::Fr};
use hash::shake128::transcript::Shake128Transcript;
use std::marker::PhantomData;

const WORD_SIZE: usize = 64;
const NUM_LIMBS: usize = 4;
const STATE_SIZE: usize = NUM_LIMBS * (WORD_SIZE / 8);
type Word = u64;

type BaseField = Fq<NUM_LIMBS>;
type ScalarField = Fr<NUM_LIMBS>;
type Curve = Pallas;
type Point = AffinePoint<BaseField, ScalarField, Curve>;
type RandomOracle = Shake128Transcript<ScalarField>;

#[derive(Debug)]
pub struct IPA<const D: usize> {
    G: [Point; D],
    H: Point,
    RO: RandomOracle,
}

impl<const K: usize, const D: usize> Commitment<K, BaseField, ScalarField, Curve> for IPA<D> {
    const GENERATOR: Point = AffinePoint {
        x: Fq(BigInt([0 as Word; NUM_LIMBS])),
        y: Fq(BigInt([0 as Word; NUM_LIMBS])),
        _p1: PhantomData,
        _p2: PhantomData,
    };
    fn setup() -> Self {
        unimplemented!()
    }
    fn commit(poly: &SparsePolynomial<ScalarField>) -> Point {
        unimplemented!()
    }
    fn prove(
        poly: &SparsePolynomial<ScalarField>,
        x: ScalarField,
    ) -> IPAProof<K, BaseField, ScalarField, Curve> {
        unimplemented!()
    }
    fn verify(cm: Point, x: ScalarField, proof: IPAProof<K, BaseField, ScalarField, Curve>) {
        unimplemented!()
    }
}
