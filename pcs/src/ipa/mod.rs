pub mod ipa_over_pallas;

// Inner Product Argument based Polynomial Commitment Scheme
use crate::polynomial::sparse::SparsePolynomial;
use ecc::elliptic_curve_arithmetic::weierstrass_model::{AffinePoint, Curve};
use ecc::finite_field_arithmetic::traits::weierstrass_field::PrimeField;

const WORD_SIZE: usize = 64;
const NUM_LIMBS: usize = 4;
const STATE_SIZE: usize = NUM_LIMBS * (WORD_SIZE / 8);
type Word = u64;

pub struct IPAProof<
    const K: usize,
    F: PrimeField<NUM_LIMBS>,
    R: PrimeField<NUM_LIMBS>,
    C: Curve<F, R>,
> {
    L: [AffinePoint<F, R, C>; K],
    R: [AffinePoint<F, R, C>; K],
    a: R,
    v: R,
}

pub trait Commitment<
    const K: usize,
    F: PrimeField<NUM_LIMBS>,
    R: PrimeField<NUM_LIMBS>,
    C: Curve<F, R>,
>
{
    // for the purpose of generating random (universal) SRS
    const GENERATOR: AffinePoint<F, R, C>;
    // generate universal SRS and new a RO instance
    fn setup() -> Self;
    // commit for target polynomial (a vector)
    fn commit(poly: &SparsePolynomial<R>) -> AffinePoint<F, R, C>;
    // generate ipa proof (non-succinct) for taget polynomial
    fn prove(poly: &SparsePolynomial<R>, x: R) -> IPAProof<K, F, R, C>;
    // verify with provided ipa proof
    fn verify(cm: AffinePoint<F, R, C>, x: R, proof: IPAProof<K, F, R, C>);
}
