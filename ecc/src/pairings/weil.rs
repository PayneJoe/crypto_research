use crate::elliptic_curve_arithmetic::models::short_weierstrass_model::{AffinePoint, Curve};
use crate::finite_field_arithmetic::bigint::BigInt;
use crate::finite_field_arithmetic::pairing_friendly::field::{Field, PrimeField};

const BASE_NUM_LIMBS: usize = 6;
const SCALAR_NUM_LIMBS: usize = 4;

pub trait WeilPairing {
    type G1BaseField: PrimeField<BASE_NUM_LIMBS>;
    type ScalarField: Field<SCALAR_NUM_LIMBS>;
    type Gt: Field<BASE_NUM_LIMBS, BasePrimeField = Self::G1BaseField>;

    type G1Curve: Curve<Self::G1BaseField, Self::ScalarField, Self::Gt>;
    type G2Curve: Curve<Self::Gt, Self::ScalarField, Self::Gt>;

    fn double_line(
        T: &mut AffinePoint<Self::Gt, Self::ScalarField, Self::Gt, Self::G2Curve>,
        Q: &AffinePoint<Self::Gt, Self::ScalarField, Self::Gt, Self::G2Curve>,
    ) -> (Self::Gt, Self::Gt) {
        let constant_3 = Self::Gt::from_small_base_prime_field_str("3");
        let constant_2 = Self::Gt::from_small_base_prime_field_str("2");
        let (x_t, y_t) = (T.x, T.y);
        let (x_q, y_q) = (Q.x, Q.y);
        let alpha = (constant_3 * x_t.square()) / (constant_2 * y_t);
        T.x = alpha.square() - constant_2 * x_t;
        T.y = -y_t - alpha * (T.x - x_t);

        (y_q - y_t - alpha * (x_q - x_t), x_q - T.x)
    }

    fn add_line(
        T: &mut AffinePoint<Self::Gt, Self::ScalarField, Self::Gt, Self::G2Curve>,
        P: &AffinePoint<Self::Gt, Self::ScalarField, Self::Gt, Self::G2Curve>,
        Q: &AffinePoint<Self::Gt, Self::ScalarField, Self::Gt, Self::G2Curve>,
    ) -> (Self::Gt, Self::Gt) {
        let (x_t, y_t) = (T.x, T.y);
        let (x_p, y_p) = (P.x, P.y);
        let (x_q, y_q) = (Q.x, Q.y);
        let alpha = (y_t - y_p) / (x_t - x_p);
        T.x = alpha.square() - x_t - x_p;
        T.y = -y_t - alpha * (T.x - x_t);

        (y_q - y_t - alpha * (x_q - x_t), x_q - T.x)
    }

    fn miller_loop(
        P: &AffinePoint<Self::Gt, Self::ScalarField, Self::Gt, Self::G2Curve>,
        Q: &AffinePoint<Self::Gt, Self::ScalarField, Self::Gt, Self::G2Curve>,
        r: &BigInt<SCALAR_NUM_LIMBS>,
    ) -> Self::Gt {
        let r_bits: Vec<u8> = r.to_bits().into_iter().rev().collect();
        let bit_length = r_bits.len();
        let mut T = P.clone();
        let (mut f1, mut f2) = (Self::Gt::ONE(), Self::Gt::ONE());

        for i in 1..bit_length {
            // last bit
            if (i == bit_length - 1) && (r_bits[i] == 0_u8) {
                f1 = f1 * (Q.x - T.x);
                T = T + T;
                break;
            }
            let (e1, e2) = Self::double_line(&mut T, &Q);
            (f1, f2) = (f1.square() * e1, f2.square() * e2);

            // last bit
            if (i == bit_length - 1) && (r_bits[i] == 1_u8) {
                f1 = f1 * (Q.x - T.x);
                T = T + *P;
                break;
            }
            if r_bits[i] == 1_u8 {
                let (e1, e2) = Self::add_line(&mut T, P, Q);
                (f1, f2) = (f1 * e1, f2 * e2);
            }
        }
        assert_eq!(T.is_identity(), true);
        f1 / f2
    }

    fn weil_pairing(
        P: &AffinePoint<Self::G1BaseField, Self::ScalarField, Self::Gt, Self::G1Curve>,
        Q: &AffinePoint<Self::Gt, Self::ScalarField, Self::Gt, Self::G2Curve>,
        r: &BigInt<SCALAR_NUM_LIMBS>,
    ) -> Self::Gt {
        // convert P from Fp into Fp12
        let px = AffinePoint::<Self::Gt, Self::ScalarField, Self::Gt, Self::G2Curve>::new(
            Self::Gt::from_base_prime_field_elem(P.x),
            Self::Gt::from_base_prime_field_elem(P.y),
        );
        let f_rp_q = Self::miller_loop(&px, Q, r);
        let f_rq_p = Self::miller_loop(Q, &px, r);

        match r.0[0] % 2 == 0 {
            true => f_rp_q / f_rq_p,
            false => -f_rp_q / f_rq_p,
        }
    }
}
