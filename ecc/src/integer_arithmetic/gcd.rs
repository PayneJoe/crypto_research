use crate::integer_arithmetic::{
    addition::Addition, division::MultiplePrecision, multiplication::SinglePrecisionMultiplication,
    substraction::Substraction, BigInteger,
};

/// lehmer matrix for lehmer extended gcd
/// referenced by Algorithm 10.46 of "Handbook of Elliptic and Hyperelliptic Curve Cryptography"
#[derive(Debug, PartialEq, Eq)]
pub struct LehmerMatrix(u8, u8, u8, u8, bool);

impl LehmerMatrix {
    pub const IDENTITY: Self = Self(1, 0, 0, 1, true);
    pub const HALF_WORD: u8 = 1 << 5;
    pub const HALF_SIZE_WORD: u8 = 1 << 3;

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
    fn euclid_extended_gcd(x: u8, N: u8) -> (u8, u8, u8, bool) {
        assert!(x < N);
        let (mut A, mut B) = (N, x);
        let (mut Ua, mut Ub) = (0, 1);
        let (mut Va, mut Vb) = (1, 0);
        let mut n_iter = 0_u8;
        while B != 0 {
            let q = A / B;
            (A, B) = (B, A - q * B);
            (Ua, Ub) = (Ub, Ua + q * Ub);
            (Va, Vb) = (Vb, Va + q * Vb);
            n_iter = n_iter + 1;
        }
        let (d, u, v) = (A, Ua, Va);
        (u, v, d, n_iter % 2 == 0)
    }

    /// one step forwards optimal gcd:
    /// euclid extended gcd with least remainder, which is approximately 30% faster than traditional euclid extended gcd
    fn euclid_extended_gcd_least_remainder(x: u8, N: u8) -> (u8, u8, u8, bool) {
        assert!(x < N);

        let (mut A, mut B) = (N, x);
        let (mut Ua, mut Ub) = (0, 1);
        let (mut Va, mut Vb) = (1, 0);
        let mut n_iter = 0_u8;
        while B != 0 {
            let mut q = A / B;
            let r = A - q * B;
            // skip cases (iterations) specially when quotient satisfy 'q = 1'
            // therefore the overall iterations must be fewer than before
            (A, B, q) = if r > B / 2 {
                n_iter = n_iter + 1;
                (B, B - r, q + 1)
            } else {
                (B, r, q)
            };
            (Ua, Ub) = (Ub, Ua + q * Ub);
            (Va, Vb) = (Vb, Va + q * Vb);
            n_iter = n_iter + 1;
        }
        let (d, u, v) = (A, Ua, Va);
        (u, v, d, n_iter % 2 == 0)
    }

    /// lehmer extended gcd for multi-precision of positive integers (big integer)
    /// referenced by Algorithm 10.45 of "Handbook of Elliptic and Hyperelliptic Curve Cryptography"
    fn lehmer_extended_gcd(x: &BigInteger, N: &BigInteger) -> (u8, u8, u8, bool) {
        let (mut A, mut B) = (N.clone(), x.clone());
        let (mut U_A, mut U_B) = (0, 1);
        let (mut V_A, mut V_B) = (1, 0);
        assert!(!B.is_zero());

        //// step 1: reduce A/B into single precisions through lehmer mixed approximation
        println!("\n Reducing begins ...");
        let mut signs: Vec<bool> = vec![];
        let mut trial = 0_u8;
        while B.data.len() > 1 {
            println!(" -------- {}th trial ---------\n", trial);
            // get the optimal lehmer matrix along with its sign
            let mat = LehmerMatrix::MixedApproximation(&A, &B);
            println!("lehmer matrix {:?}", mat);
            // update source input integers
            (A, B) = if mat == LehmerMatrix::IDENTITY {
                println!("deadlock happened ...");
                let (_, r) = A.divide_by_multiple_precision(&B);
                (B, r)
            } else {
                // sign must be considered while updating source input integers
                println!("update with lehmer matrix ...");
                if mat.4 {
                    (
                        A.multiply_single_precision(mat.0)
                            .substract(&B.multiply_single_precision(mat.1)),
                        B.multiply_single_precision(mat.3)
                            .substract(&A.multiply_single_precision(mat.2)),
                    )
                } else {
                    (
                        B.multiply_single_precision(mat.1)
                            .substract(&A.multiply_single_precision(mat.0)),
                        A.multiply_single_precision(mat.2)
                            .substract(&B.multiply_single_precision(mat.3)),
                    )
                }
            };
            println!("updated A = {:?}, B = {:?}", A, B);
            // sign must be considered while updating coefficients of input integers
            // (U_A, U_B, V_A, V_B) = if mat.4 {
            //     (
            //         U_A * mat.0 - U_B * mat.1,
            //         U_B * mat.3 - U_A * mat.2,
            //         V_A * mat.0 - V_B * mat.1,
            //         V_B * mat.3 - V_A * mat.2,
            //     )
            // } else {
            //     (
            //         U_B * mat.1 - U_A * mat.0,
            //         U_A * mat.2 - U_B * mat.3,
            //         V_B * mat.1 - V_A * mat.0,
            //         V_A * mat.2 - V_B * mat.3,
            //     )
            // };
            (U_A, U_B, V_A, V_B) = (
                U_A * mat.0 + U_B * mat.1,
                U_A * mat.2 + U_B * mat.3,
                V_A * mat.0 + V_B * mat.1,
                V_A * mat.2 + V_B * mat.3,
            );
            signs.push(mat.4);
            trial = trial + 1;
            println!("--------------------------------- \n");
        }
        // accumulate the signs for each of coefficients
        let acc_sign = signs.into_iter().reduce(|a, b| a ^ b).unwrap();
        let (sign_U_A, sign_U_B, sign_V_A, sign_V_B) = (
            false ^ acc_sign,
            true ^ acc_sign,
            true ^ acc_sign,
            false ^ acc_sign,
        );

        //// step 2: conduct euclide extended gcd
        println!("\n\n");
        println!(
            "---- After reduced through lehmer algorithm A = {:?}, B = {:?}",
            A, B
        );
        assert!(A.data.len() == B.data.len());
        let (mut u, mut v, d, mut sign) = Self::euclid_extended_gcd(B.data[0], A.data[0]);
        // sign = true => u < 0, v > 0
        // sign = false => u > 0, v < 0
        let (sign_u, sign_v) = if sign { (false, true) } else { (true, false) };

        //// step 3: coefficients combination
        sign = !(sign_U_A ^ sign_u);
        (u, v) = (U_A * u + U_B * v, V_A * u + V_B * v);

        (u, v, d, sign)
    }
}

impl GCD for BigInteger {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_not_coprime_gcd() {
        let (u, v, d, _) = BigInteger::euclid_extended_gcd(54, 189);
        // println!("---- u = {}, v = {}", u, v);
        assert_eq!(d, 27);
    }

    #[test]
    fn test_coprime_gcd() {
        let (u, v, d, _) = BigInteger::euclid_extended_gcd(45, 127);
        // println!("---- u = {}, v = {}", u, v);
        assert_eq!(d, 1);
    }

    #[test]
    fn test_lehmer_extended_gcd() {
        let mut a_arr = vec![33, 61, 20];
        let mut b_arr = vec![45, 12, 32];
        a_arr.reverse();
        b_arr.reverse();
        let a = BigInteger {
            data: a_arr,
            basis: 1 << 6,
        };
        let b = BigInteger {
            data: b_arr,
            basis: 1 << 6,
        };
        let result = BigInteger::lehmer_extended_gcd(&a, &b);
    }
}
