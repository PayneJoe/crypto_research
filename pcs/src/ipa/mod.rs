pub mod ipa_over_pallas;

// Inner Product Argument based Polynomial Commitment Scheme
use crate::polynomial::{dense::DensePolynomial, sparse::SparsePolynomial};
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

impl<const K: usize, F: PrimeField<NUM_LIMBS>, R: PrimeField<NUM_LIMBS>, C: Curve<F, R>>
    IPAProof<K, F, R, C>
{
    fn new(
        l_term: [AffinePoint<F, R, C>; K],
        r_term: [AffinePoint<F, R, C>; K],
        a_term: R,
        eval: R,
    ) -> Self {
        Self {
            L: l_term,
            R: r_term,
            a: a_term,
            v: eval,
        }
    }
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
    // commit for target sparse polynomial (a vector)
    fn commit_sparse(&self, poly: &SparsePolynomial<R>) -> AffinePoint<F, R, C>;
    // commit for target dense polynomial (a vector)
    fn commit_dense(
        poly: &DensePolynomial<R>,
        bases: &Vec<AffinePoint<F, R, C>>,
    ) -> AffinePoint<F, R, C>;
    // inner product for scalar field
    fn inner_product(lft: &Vec<R>, rht: &Vec<R>) -> R;
    // generate ipa proof (non-succinct) for taget polynomial
    fn prove(&mut self, poly: &SparsePolynomial<R>, x: R) -> IPAProof<K, F, R, C>;
    // verify with provided ipa proof
    fn verify(&mut self, aG: AffinePoint<F, R, C>, x: R, proof: &IPAProof<K, F, R, C>) -> bool;
}
