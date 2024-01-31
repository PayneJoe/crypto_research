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
        x: Fq(BigInt([7256640077462241284, 9879318615658062958, 0, 0])),
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
            let r = PallasScalarField::random();
            random_G[i] = PallasCurve::GENERATOR * r;
            // println!(
            //     "G[{}] = ({:?}, {:?}), \n G = ({:?}, {:?}), \n r = {:?} \n\n",
            //     i,
            //     random_G[i].x.rev_reduce(),
            //     random_G[i].y.rev_reduce(),
            //     PallasCurve::GENERATOR.x.rev_reduce(),
            //     PallasCurve::GENERATOR.y.rev_reduce(),
            //     r.rev_reduce(),
            // );
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
        assert!((poly.degreee() + 1) <= self.G.len());
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
        assert!((poly.degreee() + 1) <= self.G.len());
        let mut aX = poly.clone();
        let v = aX.evaluation(x);
        // let intercept_poly = PallasPoly::from_sparse_vec(vec![(0 as usize, v)]);
        // aX = aX - intercept_poly;
        // println!("------ evaluation for done, {:?} = a({:?})", v, x);

        // prepare for a, b, and G
        let mut a = DensePolynomial::from(&aX);
        a.resize(self.G.len());
        let mut b = vec![PallasScalarField::ONE()];
        for _ in 1..self.G.len() {
            let e = *b.last().unwrap();
            b.push(e * x)
        }
        let mut G = self.G.clone().to_vec();
        // println!("------ prepare a, b and G done!");

        // random a base point U for aggregation of polynomial commitment <a, G> and evaluation <a, b>
        // <a, G> + U * <a, b>
        self.RO.absorb(b"x_for_U", x.to_bytes().as_slice());
        let r_U = self.RO.squeeze(b"r_U");
        assert!(r_U.is_zero() == false);
        let U = Self::GENERATOR * r_U;
        // println!(
        //     "#[Prover Side]: random factor r_U = {:?}, U = ({:?}, {:?})",
        //     r_U,
        //     U.x.rev_reduce(),
        //     U.y.rev_reduce()
        // );
        // for i in 0..G.len() {
        //     println!(
        //         "\n #[Prover Side]: a[{}] = {:?}, b[{}] = {:?}, G[{}] = ({:?}, {:?}) \n",
        //         i,
        //         a.coefficients[i].rev_reduce(),
        //         i,
        //         b[i].rev_reduce(),
        //         i,
        //         G[i].x.rev_reduce(),
        //         G[i].y.rev_reduce()
        //     );
        // }

        /// FOR DEBUG
        let mut Ck = Self::commit_dense(&a, &G.to_vec()) + U * v;

        let mut L: Vec<PallasPoint> = Vec::new();
        let mut R: Vec<PallasPoint> = Vec::new();
        let mut size = self.G.len();
        for j in (0..self.K).rev() {
            size = size / 2;

            // split for a, b, G
            let (a_left, a_right) = (
                DensePolynomial::from(&a.coefficients[..size].to_vec()),
                DensePolynomial::from(&a.coefficients[size..].to_vec()),
            );
            let (G_left, G_right) = (G[..size].to_vec(), G[size..].to_vec());
            let (b_left, b_right) = (b[..size].to_vec(), b[size..].to_vec());

            // compute cross terms for a * G + U * (a * b)
            let (l_aG, r_aG) = (
                Self::commit_dense(&a_right, &G_left),
                Self::commit_dense(&a_left, &G_right),
            );
            let (l_ab, r_ab) = (
                Self::inner_product(&a_right.coefficients, &b_left),
                Self::inner_product(&a_left.coefficients, &b_right),
            );
            L.push(l_aG + U * l_ab);
            R.push(r_aG + U * r_ab);
            // println!("## compute cross terms for <a, G> + U * <a, b> done!");

            // challenge factor for folding a, b, and G
            // r <- hash(U, L[..i])
            let mut tmp_bytes: Vec<u8> = vec![];
            for i in 0..L.len() {
                let (label_x, label_y) = (
                    format!("{}{}-x", j, i).clone().leak(),
                    format!("{}{}-y", j, i).clone().leak(),
                );
                self.RO
                    .absorb(label_x.as_bytes(), L[i].x.to_bytes().as_slice());
                self.RO
                    .absorb(label_y.as_bytes(), L[i].y.to_bytes().as_slice());

                tmp_bytes.extend(label_x.as_bytes());
                tmp_bytes.extend(label_y.as_bytes());
                tmp_bytes.extend(L[i].x.to_bytes().as_slice());
                tmp_bytes.extend(L[i].y.to_bytes().as_slice());
            }
            let label_r = format!("{}-r", j).clone().leak();
            let r = self.RO.squeeze(label_r.as_bytes());
            let r_inv = r.inv();
            // println!(
            //     "#[Prover Side]: r[{}] = {:?}, \n r^-1 = {:?}, \n Ck = ({:?}, {:?}), \n L = ({:?}, {:?}), \n R = ({:?}, {:?}) \n ",
            //     j,
            //     r.rev_reduce(),
            //     r_inv.rev_reduce(),
            //     Ck.x.rev_reduce(),
            //     Ck.y.rev_reduce(),
            //     L[K - 1 - j].x.rev_reduce(),
            //     L[K - 1 - j].y.rev_reduce(),
            //     R[K - 1 - j].x.rev_reduce(),
            //     R[K - 1 - j].y.rev_reduce(),
            // );

            assert!(r * r_inv == PallasScalarField::ONE());

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
                G_new.push(G[i] + G[i + size] * r_inv);
            }
            (a, b, G) = (PallasDensePoly::from(&a_new), b_new, G_new);

            // for i in 0..G.len() {
            //     println!(
            //         "\n #[Prover Side]: a[{}] = {:?}, b[{}] = {:?}, G[{}] = ({:?}, {:?}) \n",
            //         i,
            //         a.coefficients[i].rev_reduce(),
            //         i,
            //         b[i].rev_reduce(),
            //         i,
            //         G[i].x.rev_reduce(),
            //         G[i].y.rev_reduce()
            //     );
            // }

            Ck = Ck + *(L.last().unwrap()) * r + *(R.last().unwrap()) * r_inv;

            // println!(
            //     "\n #[Prover Side]: Ck[{}] = ({:?}, {:?}) \n \n \n",
            //     j,
            //     Ck.x.rev_reduce(),
            //     Ck.y.rev_reduce()
            // );
            assert!(
                Ck == Self::commit_dense(&a, &G) + U * Self::inner_product(&a.coefficients, &b)
            );
        }
        assert!(L.len() == K);
        assert!(R.len() == K);
        assert!(a.degreee() == 0);
        assert!(b.len() == 1);
        assert!(G.len() == 1);

        // reset random oracle
        self.RO = RandomOracle::new(b"IPA-PCS");

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
        aG: PallasPoint,
        x: PallasScalarField,
        proof: &IPAProof<K, PallasBaseField, PallasScalarField, PallasCurve>,
    ) -> bool {
        let (L, R, a, ab) = (proof.L, proof.R, proof.a, proof.v);

        // restore U through random oracle
        self.RO.absorb(b"x_for_U", x.to_bytes().as_slice());
        let r_U = self.RO.squeeze(b"r_U");
        let U = Self::GENERATOR * r_U;
        // println!(
        //     "#[Verifier Side]: random factor r_U = {:?}, U = ({:?}, {:?})",
        //     r_U.rev_reduce(),
        //     U.x.rev_reduce(),
        //     U.y.rev_reduce()
        // );
        // self.RO.absorb(b"r_U", r_U.to_bytes().as_slice());

        // verifier need to restore b and G by himself
        let mut b = PallasScalarField::ONE();
        let mut power_of_x = x;
        let mut sX = PallasPoly::from_sparse_vec(vec![(0 as usize, PallasScalarField::ONE())]);

        // C_0 = <a, G> + U * <a, b>
        let mut Ck = aG + U * ab;

        // computing compressed commitment Ck, with time complexity O(log d) or O(k)
        let (mut r_seq, mut r_inv_seq) = (vec![], vec![]);
        for j in 0..self.K {
            // restore challenge factor for folding a, b and G
            let r_idx = K - 1 - j;
            for i in 0..(j + 1) {
                let (label_x, label_y) = (
                    format!("{}{}-x", r_idx, i).clone().leak(),
                    format!("{}{}-y", r_idx, i).clone().leak(),
                );
                self.RO
                    .absorb(label_x.as_bytes(), L[i].x.to_bytes().as_slice());
                self.RO
                    .absorb(label_y.as_bytes(), L[i].y.to_bytes().as_slice());
            }
            let label_r = format!("{}-r", r_idx).clone().leak();
            let r = self.RO.squeeze(label_r.as_bytes());
            let r_inv = r.inv();
            r_seq.push(r);
            r_inv_seq.push(r_inv);
            // println!(
            //     "#[Verifier Side]: r[{}] = {:?}, r^-1[{}] = {:?}",
            //     r_idx,
            //     r.rev_reduce(),
            //     r_idx,
            //     r_inv.rev_reduce()
            // );

            // Ck <- Ck + L[j] * r + R[j] * r^{-1}
            Ck = Ck + L[j] * r + R[j] * r_inv;

            // println!(
            //     "\n #[Verifier Side]: Ck[{}] = ({:?}, {:?})",
            //     r_idx,
            //     Ck.x.rev_reduce(),
            //     Ck.y.rev_reduce()
            // );
        }

        // restore b and G
        for j in 0..self.K {
            // b <- b * (1 + r^{-1} * b^{2^j})
            b = b * (PallasScalarField::ONE() + power_of_x * r_inv_seq[K - 1 - j]);
            power_of_x = power_of_x * power_of_x;

            // sX <- sX * (1 + r^{-1} X^{2^j})
            sX = sX
                * PallasPoly::from_sparse_vec(vec![
                    (0 as usize, PallasScalarField::ONE()),
                    ((1 as usize) << j, r_inv_seq[K - 1 - j]),
                ]);
        }
        // dense computing, with time complexity O(d)
        let G = Self::commit_dense(&PallasDensePoly::from(&sX), &self.G.to_vec());
        // println!(
        //     "\n #[Verifier Side]: b = {:?}, G = ({:?}, {:?})",
        //     b.rev_reduce(),
        //     G.x.rev_reduce(),
        //     G.y.rev_reduce()
        // );

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
    use std::str::FromStr;

    #[test]
    fn test_random() {
        let a = PallasScalarField::random();
        let b = PallasScalarField::random();
        assert_ne!(a, b);
    }

    #[test]
    fn test_setup() {
        use std::collections::HashSet;

        const K: usize = 2;
        const D: usize = 4;
        let pcs = IpaPcs::<D, K>::setup();

        let mut data = HashSet::new();
        for i in 0..D {
            assert_eq!(data.contains(&pcs.G[i]), false);
            assert_eq!(pcs.G[i].is_on_curve(), true);
            data.insert(pcs.G[i]);
        }
        assert_eq!(data.contains(&pcs.H), false);
        // data.insert(pcs.H);
    }

    #[test]
    fn test_commit_sparse() {
        const K: usize = 6;
        const D: usize = 64;
        let pcs = IpaPcs::<D, K>::setup();
        let srs = pcs.G;

        let params_a = [(1, 10), (3, 20), (5, 30)];
        let poly_a = SparsePolynomial::from_sparse_vec(
            params_a
                .map(|(idx, v)| {
                    (
                        idx as usize,
                        PallasScalarField::from_str(v.to_string().as_str()).unwrap(),
                    )
                })
                .to_vec(),
        );
        assert_eq!(poly_a.degreee(), 5);

        let comm_a = pcs.commit_sparse(&poly_a);
        let result = srs[1] * PallasScalarField::from_str("10").unwrap()
            + srs[3] * PallasScalarField::from_str("20").unwrap()
            + srs[5] * PallasScalarField::from_str("30").unwrap();
        assert_eq!(comm_a, result);
    }

    #[test]
    fn test_commit_dense() {
        const K: usize = 6;
        const D: usize = 64;
        let pcs = IpaPcs::<D, K>::setup();
        let srs = pcs.G;

        let params_a = [(1, 10), (3, 20), (5, 30)];
        let poly_a = SparsePolynomial::from_sparse_vec(
            params_a
                .map(|(idx, v)| {
                    (
                        idx as usize,
                        PallasScalarField::from_str(v.to_string().as_str()).unwrap(),
                    )
                })
                .to_vec(),
        );
        assert_eq!(poly_a.degreee(), 5);

        let comm_a = pcs.commit_sparse(&poly_a);
        let comm_b = IpaPcs::<D, K>::commit_dense(
            &PallasDensePoly::from(&poly_a),
            &srs[..(poly_a.degreee() + 1)].to_vec(),
        );
        assert_eq!(comm_a, comm_b);
    }

    #[test]
    fn test_prove() {
        const K: usize = 6;
        const D: usize = 64;
        let mut pcs = IpaPcs::<D, K>::setup();

        // y = 10 + 20 x
        let params_a = [(0, 10), (1, 20)];
        let poly_a = SparsePolynomial::from_sparse_vec(
            params_a
                .map(|(idx, v)| {
                    (
                        idx as usize,
                        PallasScalarField::from_str(v.to_string().as_str()).unwrap(),
                    )
                })
                .to_vec(),
        );
        assert_eq!(poly_a.degreee(), 1);

        let x = PallasScalarField::from_str("3").unwrap();
        pcs.prove(&poly_a, x);
    }

    #[test]
    fn test_pcs() {
        const K: usize = 2;
        const D: usize = 4;
        let mut pcs = IpaPcs::<D, K>::setup();

        let params_a = [(0, 10), (1, 20), (3, 30)];
        let poly_a = SparsePolynomial::from_sparse_vec(
            params_a
                .map(|(idx, v)| {
                    (
                        idx as usize,
                        PallasScalarField::from_str(v.to_string().as_str()).unwrap(),
                    )
                })
                .to_vec(),
        );
        assert!(poly_a.degreee() < D);

        let x = PallasScalarField::from_str("3").unwrap();
        let y = PallasScalarField::from_str("880").unwrap();
        let comm_a = pcs.commit_sparse(&poly_a);
        let proof = pcs.prove(&poly_a, x);
        assert_eq!(proof.v, y);
        let pass = pcs.verify(comm_a, x, &proof);
        assert_eq!(pass, true);
    }
}
