//// This is a practical implementation for pallas curve over standard Weierstrass model
use crate::finite_field_arithmetic::bigint::BigInt;
use crate::finite_field_arithmetic::pairing_friendly::bls12::{fq::Fq, fr::Fr};
use crate::finite_field_arithmetic::pairing_friendly::field::Field;
use std::marker::PhantomData;

use crate::elliptic_curve_arithmetic::models::short_weierstrass_model::{AffinePoint, Curve};
use crate::utils;

// make sure WINDOW_SIZE < 8
const WINDOW_SIZE: usize = 6;
const BASE_NUM_LIMBS: usize = 6;
const SCALAR_NUM_LIMBS: usize = 4;
const COFACTOR_NUM_LIMBS: usize = 8;
const WORD_SIZE: usize = 64;
type Word = u64;

// type BigInteger = BigInt<BASE_NUM_LIMBS>;

// custom curve instance
#[derive(Debug, Clone, PartialEq, Eq, Copy, Hash)]
pub struct G1;

type ScalarField = Fr<SCALAR_NUM_LIMBS>;
type BaseField = Fq<BASE_NUM_LIMBS>;

// implementation of custom curve instance
impl Curve<BaseField, ScalarField> for G1 {
    // y^2 = x^3 + 4
    // for the purpose of research, we fill these constant parameters through mannual computation
    // actually these constant parameters need to be determined dynamicly at compile-time
    const a4: Fq<BASE_NUM_LIMBS> = Fq(BigInt([0 as Word; BASE_NUM_LIMBS]));
    const a6: Fq<BASE_NUM_LIMBS> = Fq(BigInt([
        12260768510540316659,
        6038201419376623626,
        5156596810353639551,
        12813724723179037911,
        10288881524157229871,
        708830206584151678,
    ]));
    // faked identity representing the point at infinity
    const IDENTITY: AffinePoint<BaseField, ScalarField, Self> = AffinePoint {
        x: Fq(BigInt([0 as Word; BASE_NUM_LIMBS])),
        y: Fq(BigInt([0 as Word; BASE_NUM_LIMBS])),
        _p1: PhantomData,
        _p2: PhantomData,
    };
    const COFACTOR: BigInt<COFACTOR_NUM_LIMBS> =
        BigInt([10088250816726084267, 4137836090706223446, 0, 0, 0, 0, 0, 0]);
    // generator of this elliptic curve
    const GENERATOR: AffinePoint<BaseField, ScalarField, Self> = AffinePoint {
        x: Fq(BigInt([
            14075975981840609075,
            11234243599598702937,
            7066590567892948242,
            17488881670571171319,
            11959160094107449812,
            1260506501383538803,
        ])),
        y: Fq(BigInt([
            14015390473939241526,
            1418777928756678030,
            2198044176203395939,
            9785802262355231165,
            3059124926027740001,
            543439022639199291,
        ])),
        _p1: PhantomData,
        _p2: PhantomData,
    };
    // frobenius map on G1 is trival
    const FROB_TWIST_X: BaseField = Fq(BigInt([0, 0, 0, 0, 0, 0]));
    const FROB_TWIST_Y: BaseField = Fq(BigInt([0, 0, 0, 0, 0, 0]));

    fn is_on_curve(p: &AffinePoint<BaseField, ScalarField, Self>) -> bool {
        let (x, y) = (p.x, p.y);
        // let lft = y * y + Self::a1 * x * y + Self::a3 * y;
        // let rht = x * x * x + Self::a2 * x * x + Self::a4 * x + Self::a6;
        let lft = y.square();
        let rht = x.square() * x + Self::a4 * x + Self::a6;
        lft == rht
    }

    fn to_y(x: &Fq<BASE_NUM_LIMBS>) -> Fq<BASE_NUM_LIMBS> {
        unimplemented!()
    }

    fn neg(
        p: &AffinePoint<BaseField, ScalarField, Self>,
    ) -> AffinePoint<BaseField, ScalarField, Self> {
        AffinePoint {
            x: p.x,
            y: -p.y,
            _p1: PhantomData,
            _p2: PhantomData,
        }
    }

    fn is_negate(
        p1: &AffinePoint<BaseField, ScalarField, Self>,
        p2: &AffinePoint<BaseField, ScalarField, Self>,
    ) -> bool {
        // (p1.x == p2.x) && (p1.y + p2.y == -Self::a1 * p1.x - Self::a3)
        (p1.x == p2.x) && (p1.y + p2.y == BaseField::ZERO())
    }

    // referenced from P.270 of "handbook of elliptic and hyperelliptic curve cryptography"
    fn addition(
        p1: &AffinePoint<BaseField, ScalarField, Self>,
        p2: &AffinePoint<BaseField, ScalarField, Self>,
    ) -> AffinePoint<BaseField, ScalarField, Self> {
        if *p1 == Self::IDENTITY {
            return p2.clone();
        }
        if *p2 == Self::IDENTITY {
            return p1.clone();
        }
        if Self::is_negate(p1, p2) {
            return Self::IDENTITY;
        }
        let (x1, y1, x2, y2) = (p1.x, p1.y, p2.x, p2.y);
        // chord and tangle method
        let (denominator, nominator) = if p1 == p2 {
            (
                // x1 * x1 * (3 as Word) + Self::a2 * x1 * (2 as Word) + Self::a4 - Self::a1 * y1,
                // y1 * (2 as Word) + Self::a1 * x1 + Self::a3,
                x1.square() * (3 as Word) + Self::a4,
                y1 * (2 as Word),
            )
        } else {
            (y1 - y2, x1 - x2)
        };
        let lambda = denominator * nominator.inverse().unwrap();
        // FOR DEBUG
        assert!(nominator * nominator.inverse().unwrap() == BaseField::ONE());
        // let x3 = lambda * lambda + Self::a1 * lambda - Self::a2 - x1 - x2;
        let x3 = lambda.square() - x1 - x2;
        AffinePoint {
            x: x3,
            // y: lambda * (x1 - x3) - y1 - Self::a1 * x3 - Self::a3,
            y: lambda * (x1 - x3) - y1,
            _p1: PhantomData,
            _p2: PhantomData,
        }
    }

    // referenced from Algorithm 13.6 of "handbook of elliptic and hyperelliptic curve cryptography"
    fn scalar_mul(
        base: &AffinePoint<BaseField, ScalarField, Self>,
        scalar: &Fr<SCALAR_NUM_LIMBS>,
    ) -> AffinePoint<BaseField, ScalarField, Self> {
        if scalar.is_zero() {
            return AffinePoint::<BaseField, ScalarField, Self>::IDENTITY();
        }
        // k < 8, make sure Word is big enough for store precomputated points
        // let k = 6;
        assert!(WINDOW_SIZE < 8);
        let scalar_limbs: BigInt<SCALAR_NUM_LIMBS> = (*scalar).into();
        let scalar_bits: Vec<u8> = scalar_limbs.to_bits();

        // precomputation table
        let mut table = vec![base.clone()];
        let double_base = base + base;

        for i in 1..(1 << (WINDOW_SIZE - 1)) {
            let ele = &table[i - 1] + &double_base;
            table.push(ele);
        }

        let (mut q, mut i) = (Self::IDENTITY, (scalar_bits.len() - 1) as i32);
        // println!("scalar_bits = {:?}, {}", scalar_bits, scalar_bits.len());
        // let (mut times, mut old_times) = (0 as u128, 0 as u128);

        while i >= 0 {
            // println!(
            //     "##[Progress] {}/{}",
            //     (scalar_bits.len() - 1) as i32 - i,
            //     scalar_bits.len()
            // );
            // left shift skipping zeros, doubling
            if scalar_bits[i as usize] == 0 {
                (q, i) = (&q + &q, i - 1);
                // old_times = times;
                // times = times * 2;
                // println!("##{}: {}q -> {}q", i + 1, old_times, times);
            } else {
                // left shift, doubling
                let mut s = std::cmp::max((i as i32) - (WINDOW_SIZE as i32) + 1, 0);
                let right_bound = s;
                while scalar_bits[s as usize] == 0 {
                    s = s + 1;
                }
                let left_bound = s;
                for _ in 0..(i - s + 1) {
                    q = &q + &q;
                }
                // old_times = times;
                // times = times << (i - s + 1);
                // println!("##{}: {}q -> {}q", i, old_times, times);

                // then addition with precomputated table
                let u = utils::bits_to_word(&scalar_bits[(s as usize)..(i + 1) as usize]);
                // println!("i = {}, lookup u = {} ({}-{})", i, u, s, i + 1);
                assert!(u % 2 == 1);

                q = &q + &table[((u - 1) / 2) as usize];
                // old_times = times;
                // times = times + u as u128;
                // println!("##{}: {}q -> {}q", i, old_times, times);

                i = if s >= 1 { s - 1 } else { -1 };
            }
        }
        q
    }
}
