use ecc::finite_field_arithmetic::traits::{BigInt, Field, PrimeField};
use std::ops::{Add, Mul};
use std::str::FromStr;

#[derive(Debug, Clone)]
pub struct SparsePolynomial<F: Field<2>> {
    coefficients: Vec<(usize, F)>,
}

impl<F: Field<2>> SparsePolynomial<F> {
    fn from_sparse_vec(sparse_vals: Vec<(usize, F)>) -> Self {
        let mut sparse_coefficients = sparse_vals.clone();
        sparse_coefficients.sort_by(|(ida, va), (idb, vb)| ida.cmp(idb));
        Self {
            coefficients: sparse_coefficients,
        }
    }

    fn from_dense_vec(dense_vals: Vec<F>) -> Self {
        unimplemented!()
    }

    // TODO: redundant implementation through multiple pow operations
    // one-shot (double and add) might be more efficient
    fn evaluation(self, x: F) -> F {
        self.coefficients
            .iter()
            .map(|(idx, coeff)| *coeff * x.pow(BigInt::from_str(idx.to_string().as_str()).unwrap()))
            .sum()
    }
}

impl<F: Field<2>> Add for SparsePolynomial<F> {
    type Output = SparsePolynomial<F>;
    fn add(self, other: SparsePolynomial<F>) -> Self::Output {
        unimplemented!()
    }
}

impl<F: Field<2>> Mul for SparsePolynomial<F> {
    type Output = SparsePolynomial<F>;
    fn mul(self, other: SparsePolynomial<F>) -> Self::Output {
        unimplemented!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_evaluation() {
        let sparse_vals = [(1, 10), (3, 20), (5, 30)]
            .map(|(idx, v)| {
                (
                    idx as usize,
                    PrimeField::from_str(v.to_string().as_str()).unwrap(),
                )
            })
            .to_vec();
        let x = PrimeField::from_str("3").unwrap();
        let y = PrimeField::from_str("7860").unwrap();
        let poly = SparsePolynomial::from_sparse_vec(sparse_vals);
        assert_eq!(poly.evaluation(x), y);
    }
}
