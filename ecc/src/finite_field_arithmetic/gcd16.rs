use crate::finite_field_arithmetic::BigInt16;

/// lehmer matrix for lehmer extended gcd
/// referenced by Algorithm 10.46 of "Handbook of Elliptic and Hyperelliptic Curve Cryptography"
#[derive(Debug, PartialEq, Eq)]
pub struct LehmerMatrix(u16, u16, u16, u16, bool);

impl LehmerMatrix {
    pub const IDENTITY: Self = Self(1, 0, 0, 1, true);
    pub const HALF_WORD: u16 = 1 << 15;
    pub const HALF_SIZE_WORD: u16 = 1 << 8;

    fn MixedApproximation(A: &BigInt16, B: &BigInt16) -> Self {
        assert!(A.basis == B.basis);
        let b = A.basis;
        let (na, nb) = (A.size(), B.size());

        ///// step 1: get the most significant word from A and B
        let mut one_word_idx = std::cmp::max(na - 1, 0);
        let (mut hat_A, mut hat_B) = (
            A.data[one_word_idx],
            if one_word_idx + 1 > nb {
                0_u16
            } else {
                B.data[one_word_idx]
            },
        );
        // make sure the highest words are right
        assert!(hat_A >= hat_B);
        if hat_B == 0_u16 {
            return Self::IDENTITY;
        }

        // initial lehmer matrix
        let (mut alpha, mut beta, mut alpha_new, mut beta_new) = (1_u16, 0_u16, 0_u16, 1_u16);

        ///// step 2: update lehmer matrix through single precision positive integers (the most significant words)
        // first trial iteration
        let mut q = hat_A / hat_B;
        let mut T = hat_A - q * hat_B;
        let mut n_iter = 0_u16;
        // println!(
        //     "+++++++ q = {} = {} / {}, T_new = {}, half = {} ",
        //     q,
        //     hat_A,
        //     hat_B,
        //     T,
        //     Self::HALF_SIZE_WORD
        // );
        if T >= Self::HALF_SIZE_WORD {
            loop {
                // second trial iteration
                let q_new = hat_B / T;
                let T_new = hat_B - q_new * T;
                // println!(
                //     "+++++++ q_new = {} = {} / {}, T_new = {}, half = {} ",
                //     q_new,
                //     hat_B,
                //     T,
                //     T_new,
                //     Self::HALF_SIZE_WORD
                // );
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
                (q, hat_B, T) = (q_new, T, T_new);
                n_iter = n_iter + 1;
            }
        }
        // single precision failed, the loop above did not proceed more than one iteration
        if beta == 0 {
            return Self::IDENTITY;
        }

        ///// step 3: get the most two significant words from A and B
        if na < 2 {
            return Self(alpha, beta, alpha_new, beta_new, n_iter % 2 == 0);
        }

        let two_word_idx = std::cmp::max(na - 2, 0);
        let (mut hat_hat_A, mut hat_hat_B) = (
            BigInt16::new(&A.data[two_word_idx..].to_vec(), A.sign, b, None),
            if two_word_idx > nb - 1 {
                BigInt16::default()
            } else {
                BigInt16::new(&B.data[two_word_idx..].to_vec(), B.sign, b, None)
            },
        );

        // println!("++++++ update matrix through single precison: alpha = {}, beta = {}, alpha_new = {}, beta_new = {}", alpha, beta, alpha_new, beta_new);
        ///// step 4: update the two words with freshed lehmer matrix
        // println!(
        //     "+++++++ before updated hat_hat_A = {:?}, hat_hat_B = {:?}",
        //     hat_hat_A, hat_hat_B
        // );
        (hat_hat_A, hat_hat_B) = if n_iter % 2 == 0 {
            (
                &(&hat_hat_A * alpha as usize) - &(&hat_hat_B * beta as usize),
                &(&hat_hat_B * beta_new as usize) - &(&hat_hat_A * alpha_new as usize),
            )
        } else {
            (
                &(&hat_hat_B * beta as usize) - &(&hat_hat_A * alpha as usize),
                &(&hat_hat_A * alpha_new as usize) - &(&hat_hat_B * beta_new as usize),
            )
        };
        // println!(
        //     "+++++++ updated hat_hat_A = {:?}, hat_hat_B = {:?}",
        //     hat_hat_A, hat_hat_B
        // );

        ///// step 5: get the updated most significant word after refreshed two-words
        one_word_idx = std::cmp::max(hat_hat_A.size() - 1, 0);
        (hat_A, hat_B) = (
            hat_hat_A.data[one_word_idx],
            if one_word_idx > nb - 1 {
                0_u16
            } else {
                hat_hat_B.data[one_word_idx]
            },
        );
        assert!(hat_A >= hat_B);
        if hat_B == 0_u16 {
            return Self(alpha, beta, alpha_new, beta_new, n_iter % 2 == 0);
        }

        ///// step 6: update lehmer matrix through double precision positive integers
        q = hat_A / hat_B;
        T = hat_A - q * hat_B;
        // println!(
        //     "^^^^^^^^^^ q = {} = {} / {}, T_new = {}, half = {} ",
        //     q,
        //     hat_A,
        //     hat_B,
        //     T,
        //     Self::HALF_SIZE_WORD
        // );
        if T >= Self::HALF_SIZE_WORD {
            loop {
                // second trial iteration
                let q_new = hat_B / T;
                let T_new = hat_B - q_new * T;
                // println!(
                //     "^^^^^^^^ q_new = {} = {} / {}, T_new = {}, half = {} ",
                //     q_new,
                //     hat_B,
                //     T,
                //     T_new,
                //     Self::HALF_SIZE_WORD
                // );
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
    fn euclid_extended_gcd(x: u32, N: u32) -> (u32, u32, u16, bool) {
        assert!(x < N);
        let (mut A, mut B) = (N, x);
        let (mut Ua, mut Ub) = (0, 1);
        let (mut Va, mut Vb) = (1, 0);
        let mut n_iter = 0_u16;
        while B != 0 {
            let q = A / B;
            (A, B) = (B, A - q * B);
            (Ua, Ub) = (Ub, Ua + q * Ub);
            (Va, Vb) = (Vb, Va + q * Vb);
            n_iter = n_iter + 1;
        }
        let (d, u, v) = (A, Ua, Va);
        (u, v, d as u16, n_iter % 2 == 0)
    }

    /// one step forwards optimal gcd:
    /// euclid extended gcd with least remainder, which is approximately 30% faster than traditional euclid extended gcd
    fn euclid_extended_gcd_least_remainder(x: u32, N: u32) -> (u32, u32, u16, bool) {
        // make sure the smaller one is single precision, the larger one maybe double precision
        assert!((x < N) && (x < (1 << 16) as u32));

        let (mut A, mut B) = (N, x);
        let (mut Ua, mut Ub) = (0, 1);
        let (mut Va, mut Vb) = (1, 0);
        let mut n_iter = 0_u16;
        while B != 0 {
            let q = A / B;
            let r = A - q * B;
            // skip cases (iterations) specially when quotient satisfy 'q = 1'
            // therefore the overall iterations must be fewer than before
            (A, B, Ua, Ub, Va, Vb) = if r > B / 2 {
                n_iter = n_iter + 1;
                (
                    r,
                    B - r,
                    Ua + q * Ub,
                    Ub + (Ua + q * Ub),
                    Va + q * Vb,
                    Vb + (Va + q * Vb),
                )
            } else {
                (B, r, Ub, Ua + q * Ub, Vb, Va + q * Vb)
            };
            n_iter = n_iter + 1;
        }
        let (d, u, v) = (A, Ua, Va);
        (u, v, d as u16, n_iter % 2 == 0)
    }

    // specially for reduced precisions, a single precision along with a double precision at most
    fn euclid_extended_gcd_bigint(
        x: &BigInt16,
        N: &BigInt16,
    ) -> (BigInt16, BigInt16, BigInt16, bool) {
        assert!((N - x).is_positive());

        let (mut A, mut B) = (N.clone(), x.clone());
        let (mut Ua, mut Ub) = (
            BigInt16::new(vec![0_u16].as_slice(), false, 1 << 16, None),
            BigInt16::new(vec![1_u16].as_slice(), false, 1 << 16, None),
        );
        let (mut Va, mut Vb) = (
            BigInt16::new(vec![1_u16].as_slice(), false, 1 << 16, None),
            BigInt16::new(vec![0_u16].as_slice(), false, 1 << 16, None),
        );
        let mut n_iter = 0;
        while B.is_zero() == false {
            let q = &A / &B;
            let r = &A - &(&q * &B);
            (A, B, Ua, Ub, Va, Vb) = if (&r - &(&B >> 1 as usize)).is_positive() {
                n_iter = n_iter + 1;
                (
                    r.clone(),
                    &B - &r,
                    &Ua + &(&q * &Ub),
                    &Ub + &(&Ua + &(&q * &Ub)),
                    &Va + &(&q * &Vb),
                    &Vb + &(&Va + &(&q * &Vb)),
                )
            } else {
                (
                    B,
                    r,
                    Ub.clone(),
                    &Ua + &(&q * &Ub),
                    Vb.clone(),
                    &Va + &(&q * &Vb),
                )
            };
            n_iter = n_iter + 1;
        }
        let (d, u, v) = (A, Ua, Va);
        (u, v, d, n_iter % 2 == 0)
    }

    /// lehmer extended gcd for multi-precision of positive integers (big integer)
    /// referenced by Algorithm 10.45 of "Handbook of Elliptic and Hyperelliptic Curve Cryptography"
    fn lehmer_extended_gcd(x: &BigInt16, N: &BigInt16) -> (BigInt16, BigInt16, BigInt16, bool) {
        let (mut A, mut B) = (N.clone(), x.clone());
        let (mut U_A, mut U_B) = (
            BigInt16::new(vec![0_u16].as_slice(), false, 1 << 16, None),
            BigInt16::new(vec![1_u16].as_slice(), false, 1 << 16, None),
        );
        let (mut V_A, mut V_B) = (
            BigInt16::new(vec![1_u16].as_slice(), false, 1 << 16, None),
            BigInt16::new(vec![0_u16].as_slice(), false, 1 << 16, None),
        );
        assert!(!B.is_zero());

        //// step 1: reduce A/B into single precisions through lehmer mixed approximation
        // println!("\n Reducing begins ...");
        // let mut signs: Vec<bool> = vec![];
        let mut accumulated_sign = true;
        let mut trial = 0_u16;
        while B.is_zero() == false {
            // while B.size() > 1 {
            // println!(" -------- {}th trial ---------\n", trial);
            // get the optimal lehmer matrix along with its sign
            // println!("Before update lehmer matrix: A = {:?}, B = {:?}", A, B);
            let mat = LehmerMatrix::MixedApproximation(&A, &B);
            // println!("lehmer matrix {:?}", mat);
            (A, B) = if mat == LehmerMatrix::IDENTITY {
                // println!("deadlock happened, restarting outter loop is a neccessary step ...");
                // update accumulated sign
                accumulated_sign = !accumulated_sign;

                // update coefficients of gcd
                let quotient = &A / &B;
                let remainder = &A - &(&quotient * &B);
                (U_A, U_B) = (U_B.clone(), &U_A + &(&U_B * &quotient));
                (V_A, V_B) = (V_B.clone(), &V_A + &(&V_B * &quotient));

                // update source input integers
                (B, remainder)
            } else {
                // sign must be considered while updating source input integers
                // println!("lehmer matrix has already been updated ...");
                // update accumulated sign
                accumulated_sign = !(accumulated_sign ^ mat.4);

                // update coefficients of gcd
                (U_A, U_B, V_A, V_B) = (
                    &(&U_A * mat.0 as usize) + &(&U_B * mat.1 as usize),
                    &(&U_A * mat.2 as usize) + &(&U_B * mat.3 as usize),
                    &(&V_A * mat.0 as usize) + &(&V_B * mat.1 as usize),
                    &(&V_A * mat.2 as usize) + &(&V_B * mat.3 as usize),
                );

                // update source input integers
                if mat.4 {
                    (
                        &(&A * mat.0 as usize) - &(&B * mat.1 as usize),
                        &(&B * mat.3 as usize) - &(&A * mat.2 as usize),
                    )
                } else {
                    (
                        &(&B * mat.1 as usize) - &(&A * mat.0 as usize),
                        &(&A * mat.2 as usize) - &(&B * mat.3 as usize),
                    )
                }
            };
            // println!("updated input data A = {:?}, B = {:?}", A, B);
            // println!(
            //     "updated lehmer matrix alpha = {}, beta = {}, alpha_new = {}, beta_new = {}, sign = {}",
            //     mat.0, mat.1, mat.2, mat.3, mat.4
            // );
            trial = trial + 1;
            // println!("--------------------------------- \n");
        }

        // println!("\n\n");
        // println!(
        //     "---- After reduced through lehmer algorithm, A = {:?}, B = {:?}, U_A = {:?}, U_B = {:?}, V_A = {:?}, V_B = {:?},sign of lehmer matrix = {}",
        //     A, B, U_A, U_B, V_A, V_B, accumulated_sign
        // );

        // //// step 2: conduct euclide extended gcd with reduced A and B
        // let (mut u, mut v, d, reduced_sign) = Self::euclid_extended_gcd_bigint(&B, &A);
        // println!(
        //     "---- u = {:?}, v = {:?}, d = {:?}, reduced_sign = {}",
        //     u, v, d, reduced_sign
        // );
        // if reduced_sign {
        //     assert!(&(&A * &v) - &(&B * &u) == d);
        // } else {
        //     assert!(&(&B * &u) - &(&A * &v) == d);
        // }

        // //// step 3: coefficients combination
        // (u, v) = (&(&U_A * &u) + &(&U_B * &v), &(&V_A * &u) + &(&V_B * &v));
        // accumulated_sign = !(accumulated_sign ^ reduced_sign);
        // println!(
        //     "---- u = {:?}, v = {:?}, d = {:?}, accumulated_sign = {}",
        //     u, v, d, accumulated_sign
        // );

        let (u, v, d) = (U_A, V_A, A);

        (u, v, d, accumulated_sign)
    }
}

impl GCD for BigInt16 {}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_gcd_inner(a_arr: &[u16], b_arr: &[u16], c_arr: &[u16], basis: usize) {
        let (a, b, c) = (
            BigInt16::new(&a_arr, false, basis, Some(true)),
            BigInt16::new(&b_arr, false, basis, Some(true)),
            BigInt16::new(&c_arr, false, basis, Some(true)),
        );
        let (u, v, d, sign) = BigInt16::lehmer_extended_gcd(&b, &a);
        println!("u = {:?}, v = {:?}, d = {:?}, sign = {}", u, v, d, sign);
        assert_eq!(d, c);
        if sign {
            assert_eq!(&(&a * &v) - &(&b * &u), d);
        } else {
            assert_eq!(&(&b * &u) - &(&a * &v), d);
        }
    }

    #[test]
    fn test_lehmer_extended_gcd() {
        let (x, N) = (BigInt16::from("8378459450"), BigInt16::from("26498041357"));
        let (u, v, d) = (
            BigInt16::from("10055119245"),
            BigInt16::from("3179344757"),
            BigInt16::from("1"),
        );
        let (u_test, v_test, d_test, sign) = BigInt16::lehmer_extended_gcd(&x, &N);
        println!(
            "u = {:?}, v = {:?}, d = {:?}, sign = {}",
            u_test, v_test, d_test, sign
        );
        assert_eq!(u_test, u);
        assert_eq!(v_test, v);
        assert_eq!(d_test, d);
        assert_eq!(sign, false);
    }

    #[test]
    fn test_euclid_extended_gcd_1() {
        let (N, x) = (BigInt16::from("106431"), BigInt16::from("64256"));
        // let (N, x) = (BigInt16::from("197"), BigInt16::from("64"));
        let (u, v, d, sign) = BigInt16::euclid_extended_gcd_bigint(&x, &N);
        println!("u = {:?}, v = {:?}, d = {:?}, sign = {}", &u, &v, &d, sign);
        if sign {
            assert_eq!(&(&N * &v) - &(&x * &u), d);
        } else {
            assert_eq!(&(&x * &u) - &(&N * &v), d);
        }
    }

    #[test]
    fn test_euclid_extended_gcd_2() {
        // let (N, x) = (106431, 64256);
        let (N, x) = (88, 45);
        let (u, v, d, sign) = BigInt16::euclid_extended_gcd_least_remainder(x, N);
        // let (u, v, d, sign) = BigInt16::euclid_extended_gcd(x, N);
        println!("======== {}, {}, {}, {}", u, v, d, sign);
        if sign {
            assert_eq!(N as u64 * v as u64 - x as u64 * u as u64, d as u64);
        } else {
            assert_eq!(x as u64 * u as u64 - N as u64 * v as u64, d as u64);
        }
    }
}
