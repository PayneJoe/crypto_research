/// Practical Implementation of IPA-PCS specially over pallas curve
///
use crate::ipa::{Commitment, IPAProof};
use crate::polynomial::{dense::DensePolynomial, sparse::SparsePolynomial};

use ecc::elliptic_curve_arithmetic::weierstrass_model::Curve;
use ecc::elliptic_curve_arithmetic::weierstrass_model::{pallas::curve::Pallas, AffinePoint};
use ecc::finite_field_arithmetic::bigint::BigInt;
use ecc::finite_field_arithmetic::pallas::{fq::Fq, fr::Fr};
use ecc::finite_field_arithmetic::traits::weierstrass_field::PrimeField;
use hash::shake128::transcript::Shake128Transcript;
use std::marker::PhantomData;

const WORD_SIZE: usize = 64;
const NUM_LIMBS: usize = 4;
const STATE_SIZE: usize = NUM_LIMBS * (WORD_SIZE / 8);
type Word = u64;

// initiate field/point/curve/oracle/polynomial with pallas
type PallasBaseField = Fq<NUM_LIMBS>;
type PallasScalarField = Fr<NUM_LIMBS>;
type PallasCurve = Pallas;
type PallasPoint = AffinePoint<PallasBaseField, PallasScalarField, PallasCurve>;
// bytes -> scalar field (used for challenge factor)
type RandomOracle = Shake128Transcript<PallasScalarField>;
// coefficients (need to be committed) over scalar field
type PallasPoly = SparsePolynomial<PallasScalarField>;
type PallasDensePoly = DensePolynomial<PallasScalarField>;
// type PallasProof = IPAProof<const K: usize, PallasBaseField, PallasScalarField, PallasCurve> ;

#[derive(Debug, Clone)]
pub struct IpaPcs<const D: usize, const K: usize> {
    pub K: usize,
    pub G: [PallasPoint; D],
    pub H: PallasPoint,
    pub RO: RandomOracle,
}

impl<const K: usize, const D: usize> Commitment<K, PallasBaseField, PallasScalarField, PallasCurve>
    for IpaPcs<D, K>
{
    const GENERATOR: PallasPoint = AffinePoint {
        x: Fq(BigInt([
            18294172133682577413,
            12349148269572578697,
            0,
            4611686018427387904,
        ])),
        y: Fq(BigInt([
            14970995975005405177,
            1157936496307941438,
            18446744073709551615,
            4611686018427387903,
        ])),
        _p1: PhantomData,
        _p2: PhantomData,
    };

    // generate universal SRS with prepared generator
    fn setup() -> Self {
        assert!(((1 as usize) << K) == D);
        let mut random_G = [PallasCurve::IDENTITY; D];
        for i in 0..D {
            random_G[i] = PallasCurve::GENERATOR * PallasScalarField::random();
        }
        Self {
            K,
            G: random_G,
            H: PallasCurve::GENERATOR * PallasScalarField::random(),
            RO: RandomOracle::new(b"IPA-PCS"),
        }
    }

    // commit a sparse polynomial with universal SRS
    // todo!: can be parallized
    fn commit_sparse(&self, poly: &PallasPoly) -> PallasPoint {
        assert!((poly.degreee() + 1) < self.G.len());
        poly.coefficients
            .iter()
            .fold(PallasCurve::IDENTITY, |acc, (idx, coeff)| {
                acc + self.G[*idx] * *coeff
            })
    }

    // commit a dense polynomial with universal SRS
    // todo!: can be parallized
    fn commit_dense(poly: &PallasDensePoly, bases: &Vec<PallasPoint>) -> PallasPoint {
        assert!(poly.degreee() + 1 == bases.len());
        poly.coefficients
            .iter()
            .enumerate()
            .fold(PallasCurve::IDENTITY, |acc, (idx, coeff)| {
                acc + bases[idx] * *coeff
            })
    }

    // inner product for two scalar field
    // todo!: can be parallized
    fn inner_product(
        lft: &Vec<PallasScalarField>,
        rht: &Vec<PallasScalarField>,
    ) -> PallasScalarField {
        lft.iter()
            .zip(rht.iter())
            .fold(PallasScalarField::ZERO(), |acc, (a, b)| acc + *a * *b)
    }

    // generator proof for the specified opening point x
    fn prove(
        &mut self,
        poly: &PallasPoly,
        x: PallasScalarField,
    ) -> IPAProof<K, PallasBaseField, PallasScalarField, PallasCurve> {
        assert!((poly.degreee() + 1) < self.G.len());
        let mut aX = poly.clone();
        let v = aX.evaluation(x);
        // let intercept_poly = PallasPoly::from_sparse_vec(vec![(0 as usize, v)]);
        // aX = aX - intercept_poly;

        // prepare for a, b, and G
        let mut b = vec![x.clone()];
        for _ in 1..self.G.len() {
            let e = *b.last().unwrap();
            b.push(e * e)
        }
        let mut a = DensePolynomial::from(&aX);
        let mut G = self.G.clone().to_vec();

        // random a base point U for aggregation of polynomial commitment <a, G> and evaluation <a, b>
        // <a, G> + U * <a, b>
        self.RO.absorb(b"x_for_U", x.to_bytes().as_slice());
        let r_U = self.RO.squeeze(b"r_U");
        let U = Self::GENERATOR * r_U;
        self.RO.absorb(b"r_U", r_U.to_bytes().as_slice());

        let mut L: Vec<PallasPoint> = Vec::new();
        let mut R: Vec<PallasPoint> = Vec::new();
        let mut size = self.G.len();
        for j in (0..self.K).rev() {
            size = size / 2;

            // split for a, b, G
            let (a_left, a_right) = (
                DensePolynomial::from(a.coefficients[..size].to_vec()),
                DensePolynomial::from(a.coefficients[size..].to_vec()),
            );
            let (G_left, G_right) = (G[..size].to_vec(), G[size..].to_vec());
            let (b_left, b_right) = (b[..size].to_vec(), b[size..].to_vec());

            // compute cross terms for a * G + U * (a * b)
            let l_aG = Self::commit_dense(&a_right, &G_left);
            let r_aG = Self::commit_dense(&a_left, &G_right);
            let l_ab = Self::inner_product(&a_right.coefficients, &b_left);
            let r_ab = Self::inner_product(&a_left.coefficients, &b_right);
            L.push(l_aG + U * l_ab);
            R.push(r_aG + U * r_ab);

            // challenge factor for folding a, b, and G
            // r <- hash(U, L[..i])
            for i in 0..L.len() {
                let label_x = format!("{}{}-x", j, i).clone().leak();
                self.RO
                    .absorb(label_x.as_bytes(), L[i].x.to_bytes().as_slice());
                let label_y = format!("{}{}-y", j, i).clone().leak();
                self.RO
                    .absorb(label_y.as_bytes(), L[i].y.to_bytes().as_slice());
            }
            let label_r = format!("{}-r", j).clone().leak();
            let r = self.RO.squeeze(label_r.as_bytes());
            let r_inv = r.inv();

            // fold a, b, and G
            // a <- a_L + r * a_R
            // b <- b_L + r^{-1} * b_R
            // G <- G_L + r^{-1} * G_R
            let mut a_new: Vec<PallasScalarField> = Vec::new();
            let mut b_new: Vec<PallasScalarField> = Vec::new();
            let mut G_new: Vec<PallasPoint> = Vec::new();
            for i in 0..size {
                a_new.push(a.coefficients[i] + a.coefficients[i + size] * r);
                b_new.push(b[i] + b[i + size] * r_inv);
                G_new.push(G[i] + G[i + size] * r_inv)
            }
            a = PallasDensePoly::from(a_new);
            b = b_new;
            G = G_new;
        }
        L.reverse();
        R.reverse();
        assert!(a.degreee() == 0);
        assert!(b.len() == 1);
        assert!(G.len() == 1);

        IPAProof::<K, PallasBaseField, PallasScalarField, PallasCurve>::new(
            vec2array_uncheck::<PallasPoint, K>(L),
            vec2array_uncheck::<PallasPoint, K>(R),
            a.coefficients[0],
            v,
        )
    }

    // verify the evaluation with provided proof, at the cost of non-linear time complexity
    // since the time cost of compression of universal SRS is linear, so this is not succinct
    // it will be much more efficient when multiple proof need to be verified (batching for multi-opening)
    fn verify(
        &mut self,
        cm_aG: PallasPoint,
        x: PallasScalarField,
        proof: IPAProof<K, PallasBaseField, PallasScalarField, PallasCurve>,
    ) -> bool {
        let (L, R, a, v) = (proof.L, proof.R, proof.a, proof.v);

        // restore U through random oracle
        self.RO.absorb(b"x_for_U", x.to_bytes().as_slice());
        let r_U = self.RO.squeeze(b"r_U");
        let U = Self::GENERATOR * r_U;
        self.RO.absorb(b"r_U", r_U.to_bytes().as_slice());

        // verifier need to restore b and G by himself
        let mut b = PallasScalarField::ONE();
        let mut power_of_b = b;
        let mut sX = PallasPoly::from_sparse_vec(vec![(0 as usize, PallasScalarField::ONE())]);

        // C_0 = <a, G> + U * <a, b>
        let mut Ck = cm_aG + U * v;
        // computing compressed commitment Ck, restoring b and G, with time complexity O(log d) or O(k)
        for j in 0..self.K {
            // restore challenge factor for folding a, b and G
            let mut L_rev = L[j..].to_vec();
            L_rev.reverse();
            for i in 0..L_rev.len() {
                let label_x = format!("{}{}-x", j, i).clone().leak();
                self.RO
                    .absorb(label_x.as_bytes(), L[i].x.to_bytes().as_slice());
                let label_y = format!("{}{}-y", j, i).clone().leak();
                self.RO
                    .absorb(label_y.as_bytes(), L[i].y.to_bytes().as_slice());
            }
            let label_r = format!("{}-r", j).clone().leak();
            let r = self.RO.squeeze(label_r.as_bytes());
            let r_inv = r.inv();

            // Ck <- Ck + L[j] * r + R[j] * r^{-1}
            Ck = Ck + L[j] * r + R[j] * r_inv;

            // b <- b * (1 + r^{-1} * b^{2^j})
            b = b * (PallasScalarField::ONE() + power_of_b * r_inv);
            power_of_b = power_of_b * power_of_b;

            // sX <- sX * (1 + r^{-1} X^{2^j})
            sX = sX
                * PallasPoly::from_sparse_vec(vec![
                    (0 as usize, PallasScalarField::ONE()),
                    ((1 as usize) << j, r_inv),
                ]);
        } // computing Ck and restoring for b are all done here

        // restore G
        // dense computing, with time complexity O(d)
        let G = Self::commit_dense(&PallasDensePoly::from(&sX), &self.G.to_vec());

        // C0 = <a, G> + <a, b> * U
        // Ck <- \sum_{i = 0..k} {Ci + r * Ci_left + r^{-1} * Ci_right}
        // Ck == a * G + a * b * U
        Ck == G * a + U * a * b
    }
}

fn vec2array_uncheck<T, const N: usize>(v: Vec<T>) -> [T; N] {
    v.try_into()
        .unwrap_or_else(|v: Vec<T>| panic!("Expected a Vec of length {} but it was {}", N, v.len()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_random() {
        let a = PallasScalarField::random();
        let b = PallasScalarField::random();
        assert_ne!(a, b);
    }
}
