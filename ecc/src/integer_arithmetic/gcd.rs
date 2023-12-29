use crate::integer_arithmetic::{
    addition::Addition, division::MultiplePrecision, multiplication::SinglePrecisionMultiplication,
    substraction::Substraction, BigInteger,
};

/// lehmer matrix for lehmer extended gcd
/// referenced by Algorithm 10.46 of "Handbook of Elliptic and Hyperelliptic Curve Cryptography"
#[derive(PartialEq, Eq)]
pub struct LehmerMatrix(u8, u8, u8, u8, bool);

impl LehmerMatrix {
    pub const IDENTITY: Self = Self(1, 0, 0, 1, true);
    pub const HALF_WORD: u8 = 1 << 7;
    pub const HALF_SIZE_WORD: u8 = 1 << 4;

    fn MixedApproximation(A: &BigInteger, B: &BigInteger) -> Self {
        assert!(A.basis == B.basis);
        assert!(A.greater(&B));
        let (na, nb) = (A.data.len(), B.data.len());

        ///// step 1: get the most significant word from A and B
        let mut one_word_idx = std::cmp::max(na - 1, 0);
        let (mut hat_A, mut hat_B) = (
            A.data[one_word_idx],
            if one_word_idx + 1 > nb {
                0_u8
            } else {
                B.data[one_word_idx]
            },
        );
        // make sure the highest words are right
        assert!(hat_A >= hat_B);
        if hat_B == 0_u8 {
            return Self::IDENTITY;
        }

        // initial lehmer matrix
        let (mut alpha, mut beta, mut alpha_new, mut beta_new) = (1_u8, 0_u8, 0_u8, 1_u8);

        ///// step 2: update lehmer matrix through single precision positive integers (the most significant words)
        // first trial iteration
        let mut n_iter = 0_u8;
        let mut q = hat_A / hat_B;
        let mut T = hat_A - q * hat_B;
        if T >= Self::HALF_SIZE_WORD {
            while true {
                // second trial iteration
                let q_new = hat_B / T;
                let T_new = hat_B - q_new * T;
                if T_new < Self::HALF_SIZE_WORD {
                    break;
                }
                // update lehmer matrix with the first quotient (is 'q', not 'q_new'), once the second trial passed
                (alpha, beta, alpha_new, beta_new) = (
                    alpha_new,
                    beta_new,
                    alpha + q * alpha_new,
                    beta + q * beta_new,
                );
                // update quotient with the new one 'q_new', and replace remainder with the new one 'T_new'
                (q, T) = (q_new, T_new);
                n_iter = n_iter + 1;
            }
        }
        // single precision failed, upper loop did not proceed more than one iteration
        if beta == 0 {
            return Self::IDENTITY;
        }

        ///// step 3: get the most two significant words from A and B
        let two_word_idx = std::cmp::max(na - 2, 0);
        let (mut hat_hat_A, mut hat_hat_B) = (
            BigInteger {
                data: A.data[two_word_idx..].to_vec(),
                basis: A.basis,
            },
            if two_word_idx + 1 > nb {
                BigInteger {
                    data: vec![0_u8],
                    basis: A.basis,
                }
            } else {
                BigInteger {
                    data: B.data[two_word_idx..].to_vec(),
                    basis: A.basis,
                }
            },
        );

        ///// step 4: update the two words with freshed lehmer matrix
        (hat_hat_A, hat_hat_B) = (
            hat_hat_A
                .multiply_single_precision(alpha)
                .add(&hat_hat_B.multiply_single_precision(beta)),
            hat_hat_A
                .multiply_single_precision(alpha_new)
                .add(&hat_hat_B.multiply_single_precision(beta_new)),
        );

        ///// step 5: get the updated most significant word after refreshed two-words
        one_word_idx = std::cmp::max(hat_hat_A.data.len() - 1, 0);
        (hat_A, hat_B) = (
            hat_hat_A.data[one_word_idx],
            if one_word_idx + 1 > nb {
                0_u8
            } else {
                hat_hat_B.data[one_word_idx]
            },
        );
        assert!(hat_A >= hat_B);
        if hat_B == 0_u8 {
            return Self::IDENTITY;
        }

        ///// step 6: update lehmer matrix through double precision positive integers
        q = hat_A / hat_B;
        T = hat_A - q * hat_B;
        if T >= Self::HALF_SIZE_WORD {
            while true {
                // second trial iteration
                let q_new = hat_B / T;
                let T_new = hat_B - q_new * T;
                if T_new < Self::HALF_SIZE_WORD {
                    break;
                }
                // update lehmer matrix with the first quotient (is 'q', not 'q_new'), once the second trial passed
                (alpha, beta, alpha_new, beta_new) = (
                    alpha_new,
                    beta_new,
                    alpha + q * alpha_new,
                    beta + q * beta_new,
                );
                // update quotient with the new one 'q_new', and replace remainder with the new one 'T_new'
                (q, T) = (q_new, T_new);
                n_iter = n_iter + 1;
            }
        }
        Self(alpha, beta, alpha_new, beta_new, n_iter % 2 == 0)
    }
}

/// gcd for multiple precision of positive integers
pub trait GCD {
    /// native method of euclid extended gcd, many more methods can be applied to achieve few iterations
    /// referenced by Algorithm 10.42 of "Handbook of Elliptic and Hyperelliptic Curve Cryptography"
    fn euclid_extended_gcd(x: u8, N: u8) -> (u8, u8, u8) {
        assert!(x < N);
        let (mut A, mut B) = (N, x);
        let (mut Ua, mut Ub) = (0, 1);
        let (mut Va, mut Vb) = (1, 0);
        while B != 0 {
            let q = A / B;
            (A, B) = (B, A - q * B);
            (Ua, Ub) = (Ub, Ua - q * Ub);
            (Va, Vb) = (Vb, Va - q * Vb);
        }
        let (d, u, v) = (A, Ua, Va);
        (u, v, d)
    }

    /// one step forwards optimal gcd:
    /// euclid extended gcd with least remainder, which is approximately 30% faster than traditional euclid extended gcd
    fn euclid_extended_gcd_least_remainder(x: u8, N: u8) -> (u8, u8, u8) {
        assert!(x < N);

        let (mut A, mut B) = (N, x);
        let (mut Ua, mut Ub) = (0, 1);
        let (mut Va, mut Vb) = (1, 0);
        while B != 0 {
            let mut q = A / B;
            let r = A - q * B;
            // skip cases (iterations) specially when quotient satisfy 'q = 1'
            // therefore the overall iterations must be fewer than before
            (A, B, q) = if r > B / 2 {
                (B, B - r, q + 1)
            } else {
                (B, r, q)
            };
            (Ua, Ub) = (Ub, Ua - q * Ub);
            (Va, Vb) = (Vb, Va - q * Vb);
        }
        let (d, u, v) = (A, Ua, Va);
        (u, v, d)
    }

    /// lehmer extended gcd for multi-precision of positive integers (big integer)
    /// referenced by Algorithm 10.45 of "Handbook of Elliptic and Hyperelliptic Curve Cryptography"
    fn lehmer_extended_gcd(x: BigInteger, N: BigInteger) -> (u8, u8, u8) {
        let (mut A, mut B) = (N, x);
        let (mut U_A, mut U_B) = (0, 1);
        let (mut V_A, mut V_B) = (1, 0);
        assert!(!B.is_zero());

        //// step 1: reduce A/B into single precisions through lehmer mixed approximation
        while B.data.len() > 1 {
            // get the optimal lehmer matrix
            let mat = LehmerMatrix::MixedApproximation(&A, &B);
            // update source input integers
            (A, B) = if mat == LehmerMatrix::IDENTITY {
                let (_, r) = A.divide_by_multiple_precision(&B);
                (B, r)
            } else {
                (
                    A.multiply_single_precision(mat.0)
                        .add(&B.multiply_single_precision(mat.1)),
                    A.multiply_single_precision(mat.2)
                        .add(&B.multiply_single_precision(mat.3)),
                )
            };
            // update coefficients of input integers
            (U_A, U_B) = (U_A * mat.0 + U_B * mat.1, U_A * mat.2 + U_B * mat.3);
            (V_A, V_B) = (V_A * mat.0 + V_B * mat.1, V_A * mat.2 + V_B * mat.3);
        }

        //// step 2: conduct euclide extended gcd
        assert!(A.data.len() == B.data.len());
        let (mut u, mut v, d) = Self::euclid_extended_gcd(A.data[0], B.data[0]);

        //// step 3: coefficients combination
        (u, v) = (u * U_A + v * U_B, u * V_A + v * V_B);
        (u, v, d)
    }
}

impl GCD for BigInteger {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_not_coprime_gcd() {
        let (u, v, d) = BigInteger::euclid_extended_gcd(54, 188);
        assert_eq!(d, 6);
    }

    #[test]
    fn test_coprime_gcd() {
        let (u, v, d) = BigInteger::euclid_extended_gcd(45, 127);
        assert_eq!(d, 1);
    }
}
