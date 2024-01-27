use ecc::finite_field_arithmetic::bigint::BigInt;
use ecc::finite_field_arithmetic::traits::weierstrass_field::PrimeField;
use std::collections::BTreeMap;
use std::ops::{Add, Mul};
use std::str::FromStr;

const WORD_SIZE: usize = 64;
const NUM_LIMBS: usize = 4;
const STATE_SIZE: usize = NUM_LIMBS * (WORD_SIZE / 8);
type Word = u64;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SparsePolynomial<F: PrimeField<NUM_LIMBS>> {
    coefficients: Vec<(usize, F)>,
}

impl<F: PrimeField<NUM_LIMBS>> SparsePolynomial<F> {
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
            .map(|(idx, coeff)| {
                *coeff * x.pow(BigInt::<NUM_LIMBS>::from_str(idx.to_string().as_str()).unwrap())
            })
            .sum()
    }

    fn is_zero(&self) -> bool {
        self.coefficients.len() == 0
    }

    fn ZERO() -> Self {
        Self {
            coefficients: Vec::new(),
        }
    }

    fn degreee(&self) -> usize {
        if self.is_zero() {
            return 0;
        }
        self.coefficients.last().unwrap().0
    }
}

impl<F: PrimeField<NUM_LIMBS>> Add for SparsePolynomial<F> {
    type Output = SparsePolynomial<F>;
    fn add(self, other: SparsePolynomial<F>) -> Self::Output {
        &self + &other
    }
}

impl<F: PrimeField<NUM_LIMBS>> Mul for SparsePolynomial<F> {
    type Output = SparsePolynomial<F>;
    fn mul(self, other: SparsePolynomial<F>) -> Self::Output {
        &self * &other
    }
}

// merge two polynomials into one
impl<'a, 'b, F: PrimeField<NUM_LIMBS>> Add<&'a SparsePolynomial<F>> for &'b SparsePolynomial<F> {
    type Output = SparsePolynomial<F>;
    fn add(self, other: &'a SparsePolynomial<F>) -> Self::Output {
        if self.is_zero() {
            return other.clone();
        }
        if other.is_zero() {
            return self.clone();
        }
        let (mut pos_a, mut pos_b) = (0, 0);
        let mut result: Vec<(usize, F)> = Vec::new();
        loop {
            if pos_a == self.coefficients.len() {
                for i in pos_b..other.coefficients.len() {
                    result.push(other.coefficients[i])
                }
                break;
            }
            if pos_b == other.coefficients.len() {
                for i in pos_a..self.coefficients.len() {
                    result.push(self.coefficients[i])
                }
                break;
            }
            let (idx_a, coeff_a) = self.coefficients[pos_a];
            let (idx_b, coeff_b) = other.coefficients[pos_b];
            if idx_a == idx_b {
                result.push((idx_a, coeff_a + coeff_b));
                (pos_a, pos_b) = (pos_a + 1, pos_b + 1);
            } else if idx_a > idx_b {
                result.push((idx_b, coeff_b));
                pos_b = pos_b + 1;
            } else if idx_a < idx_b {
                result.push((idx_a, coeff_a));
                pos_a = pos_a + 1;
            }
        }
        SparsePolynomial::<F>::from_sparse_vec(result)
    }
}

// naive implementation of polynomial multiplication
impl<'a, 'b, F: PrimeField<NUM_LIMBS>> Mul<&'a SparsePolynomial<F>> for &'b SparsePolynomial<F> {
    type Output = SparsePolynomial<F>;
    fn mul(self, other: &'a SparsePolynomial<F>) -> Self::Output {
        if self.is_zero() || other.is_zero() {
            return SparsePolynomial::<F>::ZERO();
        }
        let mut result: BTreeMap<usize, F> = BTreeMap::new();
        for i in 0..self.coefficients.len() {
            for j in 0..other.coefficients.len() {
                let mut ent = result
                    .entry(self.coefficients[i].0 + other.coefficients[j].0)
                    .or_insert(F::ZERO());
                *ent = *ent + self.coefficients[i].1 * other.coefficients[j].1;
            }
        }
        SparsePolynomial::from_sparse_vec(result.into_iter().collect::<Vec<(usize, F)>>())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ecc::finite_field_arithmetic::pallas::{fq::Fq, fr::Fr};
    use std::str::FromStr;
    type ScalarField = Fr<NUM_LIMBS>;
    type BaseField = Fq<NUM_LIMBS>;

    #[test]
    fn test_evaluation() {
        let poly_params = [(1, 10), (3, 20), (5, 30)];
        let (x, y) = (3, 7860);

        let sparse_vals = poly_params
            .map(|(idx, v)| {
                (
                    idx as usize,
                    ScalarField::from_str(v.to_string().as_str()).unwrap(),
                )
            })
            .to_vec();
        let point = ScalarField::from_str(x.to_string().as_str()).unwrap();
        let eval = ScalarField::from_str(y.to_string().as_str()).unwrap();
        let poly = SparsePolynomial::from_sparse_vec(sparse_vals);
        assert_eq!(poly.evaluation(point), eval);
    }

    #[test]
    fn test_add() {
        let (params_a, params_b, params_c) = (
            [(1, 10), (3, 20), (5, 30)],
            [(2, 10), (3, 20), (6, 30)],
            [(1, 10), (2, 10), (3, 40), (5, 30), (6, 30)],
        );
        let poly_a = SparsePolynomial::from_sparse_vec(
            params_a
                .map(|(idx, v)| {
                    (
                        idx as usize,
                        ScalarField::from_str(v.to_string().as_str()).unwrap(),
                    )
                })
                .to_vec(),
        );
        let poly_b = SparsePolynomial::from_sparse_vec(
            params_b
                .map(|(idx, v)| {
                    (
                        idx as usize,
                        ScalarField::from_str(v.to_string().as_str()).unwrap(),
                    )
                })
                .to_vec(),
        );
        let poly_c = SparsePolynomial::from_sparse_vec(
            params_c
                .map(|(idx, v)| {
                    (
                        idx as usize,
                        ScalarField::from_str(v.to_string().as_str()).unwrap(),
                    )
                })
                .to_vec(),
        );
        assert_eq!(poly_a + poly_b, poly_c);
    }
    #[test]
    fn test_mul() {
        let (params_a, params_b, params_c) = (
            [(1, 10), (3, 20), (5, 30)],
            [(2, 10), (3, 20), (6, 30)],
            [
                (3, 100),
                (4, 200),
                (7, 600),
                (5, 200),
                (6, 400),
                (9, 600),
                (8, 600),
                (11, 900),
            ],
        );
        let poly_a = SparsePolynomial::from_sparse_vec(
            params_a
                .map(|(idx, v)| {
                    (
                        idx as usize,
                        ScalarField::from_str(v.to_string().as_str()).unwrap(),
                    )
                })
                .to_vec(),
        );
        let poly_b = SparsePolynomial::from_sparse_vec(
            params_b
                .map(|(idx, v)| {
                    (
                        idx as usize,
                        ScalarField::from_str(v.to_string().as_str()).unwrap(),
                    )
                })
                .to_vec(),
        );
        let poly_c = SparsePolynomial::from_sparse_vec(
            params_c
                .map(|(idx, v)| {
                    (
                        idx as usize,
                        ScalarField::from_str(v.to_string().as_str()).unwrap(),
                    )
                })
                .to_vec(),
        );
        assert_eq!(poly_a * poly_b, poly_c);
    }
}
