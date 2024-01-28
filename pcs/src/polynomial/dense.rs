use crate::polynomial::sparse::SparsePolynomial;
use ecc::finite_field_arithmetic::bigint::BigInt;
use ecc::finite_field_arithmetic::traits::weierstrass_field::PrimeField;
use std::collections::BTreeMap;

const WORD_SIZE: usize = 64;
const NUM_LIMBS: usize = 4;
const STATE_SIZE: usize = NUM_LIMBS * (WORD_SIZE / 8);
type Word = u64;

pub struct DensePolynomial<F: PrimeField<NUM_LIMBS>> {
    pub coefficients: Vec<F>,
}

impl<F: PrimeField<NUM_LIMBS>> DensePolynomial<F> {
    pub fn is_zero(&self) -> bool {
        self.coefficients.len() == 0
    }

    pub fn ZERO() -> Self {
        Self {
            coefficients: Vec::new(),
        }
    }

    pub fn degreee(&self) -> usize {
        self.coefficients.len() - 1
    }
}

impl<F: PrimeField<NUM_LIMBS>> From<&SparsePolynomial<F>> for DensePolynomial<F> {
    fn from(poly: &SparsePolynomial<F>) -> Self {
        let idx2coeff: BTreeMap<usize, F> = poly
            .coefficients
            .iter()
            .map(|(idx, coeff)| (*idx, *coeff))
            .collect();

        let d = poly.degreee();
        let dense_vec: Vec<F> = (0..(d + 1))
            .map(|i| {
                if idx2coeff.contains_key(&i) {
                    idx2coeff.get(&i).unwrap().clone()
                } else {
                    F::ZERO()
                }
            })
            .collect();

        DensePolynomial {
            coefficients: dense_vec,
        }
    }
}

impl<F: PrimeField<NUM_LIMBS>> From<Vec<F>> for DensePolynomial<F> {
    fn from(vec: Vec<F>) -> Self {
        Self { coefficients: vec }
    }
}
