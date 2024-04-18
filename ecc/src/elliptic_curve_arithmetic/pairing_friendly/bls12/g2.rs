//// This is a practical implementation for pallas curve over standard Weierstrass model
use crate::finite_field_arithmetic::bigint::BigInt;
use crate::finite_field_arithmetic::pairing_friendly::bls12::{
    fq::Fq, fq12::Fq12, fq2::Fq2, fq6::Fq6, fr::Fr,
};
use crate::finite_field_arithmetic::pairing_friendly::field::Field;

use crate::elliptic_curve_arithmetic::models::short_weierstrass_model::{AffinePoint, Curve};
use crate::elliptic_curve_arithmetic::pairing_friendly::bls12::g12::G12;

use crate::utils;
use std::marker::PhantomData;

// make sure WINDOW_SIZE < 8
const WINDOW_SIZE: usize = 6;
const BASE_NUM_LIMBS: usize = 6;
const SCALAR_NUM_LIMBS: usize = 4;
const WORD_SIZE: usize = 64;
type Word = u64;

type FullExtensionField = Fq12;
type FullExtensionCurve = G12;

// type BigInteger = BigInt<BASE_NUM_LIMBS>;

// custom curve instance
#[derive(Debug, Clone, PartialEq, Eq, Copy, Hash)]
pub struct G2;

type ScalarField = Fr<SCALAR_NUM_LIMBS>;
type BaseField = Fq2;

const Fq2ZERO: Fq2 = Fq2 {
    c0: Fq(BigInt([0; BASE_NUM_LIMBS])),
    c1: Fq(BigInt([0; BASE_NUM_LIMBS])),
};

const Fq6ZERO: Fq6 = Fq6 {
    c0: Fq2ZERO,
    c1: Fq2ZERO,
    c2: Fq2ZERO,
};

const Fq12ZERO: Fq12 = Fq12 {
    c0: Fq6ZERO,
    c1: Fq6ZERO,
};

impl G2 {}

// implementation of custom curve instance
impl Curve<BaseField, ScalarField, FullExtensionField> for G2 {
    // y^2 = x^3 + 4 * (u + 1)
    // for the purpose of research, we fill these constant parameters through mannual computation
    // actually these constant parameters need to be determined dynamicly at compile-time
    const a4: BaseField = BaseField {
        c0: Fq(BigInt([0; BASE_NUM_LIMBS])),
        c1: Fq(BigInt([0; BASE_NUM_LIMBS])),
    };
    const a6: BaseField = BaseField {
        c0: Fq(BigInt([
            12260768510540316659,
            6038201419376623626,
            5156596810353639551,
            12813724723179037911,
            10288881524157229871,
            708830206584151678,
        ])),
        c1: Fq(BigInt([
            12260768510540316659,
            6038201419376623626,
            5156596810353639551,
            12813724723179037911,
            10288881524157229871,
            708830206584151678,
        ])),
    };
    // faked identity representing the point at infinity
    const IDENTITY: AffinePoint<BaseField, ScalarField, FullExtensionField, Self> = AffinePoint {
        x: Fq2ZERO,
        y: Fq2ZERO,
        _p1: PhantomData,
        _p2: PhantomData,
        _p3: PhantomData,
    };
    const COFACTOR: BigInt<8> = BigInt([
        14923865813284960485,
        1591719477533150320,
        2401401741857165742,
        11973085464847352559,
        12000439733073903071,
        14812865038524111546,
        656769601966088706,
        420316534768199665,
    ]);
    // generator of this elliptic curve, g = (-1, 2)
    const GENERATOR: AffinePoint<BaseField, ScalarField, FullExtensionField, Self> = AffinePoint {
        x: BaseField {
            c0: Fq(BigInt([
                17190305259240202339,
                7296294217700534644,
                435183553914099211,
                11694783461588566454,
                14099298684380714084,
                1677353267918785999,
            ])),
            c1: Fq(BigInt([
                5685763873039570767,
                6419095246139654681,
                14293532027625969943,
                16410173736989223150,
                10929602565445287157,
                1191836563837268320,
            ])),
        },

        y: BaseField {
            c0: Fq(BigInt([
                16570605391750952557,
                10487159702080788072,
                144171039000839840,
                10179971595502506675,
                8906873232918241193,
                110972420657261769,
            ])),
            c1: Fq(BigInt([
                2470239471287016936,
                5401817021023576796,
                16653359143730736950,
                7660848434194476848,
                9728435152041756584,
                1780737560784417686,
            ])),
        },
        _p1: PhantomData,
        _p2: PhantomData,
        _p3: PhantomData,
    };

    // 1 / beta_t_x
    const TWIST_X_INV: Fq12 = Fq12 {
        c0: Fq6 {
            c0: Fq2ZERO,
            c1: Fq2ZERO,
            c2: Fq2 {
                c0: Fq(BigInt([
                    8505329371266088957,
                    17002214543764226050,
                    6865905132761471162,
                    8632934651105793861,
                    6631298214892334189,
                    1582556514881692819,
                ])),
                c1: Fq(BigInt([0; BASE_NUM_LIMBS])),
            },
        },
        c1: Fq6 {
            c0: Fq2ZERO,
            c1: Fq2ZERO,
            c2: Fq2 {
                c0: Fq(BigInt([
                    11671922859260663127,
                    11050707557586042878,
                    284884720401305268,
                    17749945728364010941,
                    8613774818643959860,
                    145621051382923523,
                ])),
                c1: Fq(BigInt([0; BASE_NUM_LIMBS])),
            },
        },
    };

    // 1 / beta_t_y
    const TWIST_Y_INV: Fq12 = Fq12 {
        c0: Fq6 {
            c0: Fq2ZERO,
            c1: Fq2 {
                c0: Fq(BigInt([0; BASE_NUM_LIMBS])),
                c1: Fq(BigInt([
                    8505329371266088957,
                    17002214543764226050,
                    6865905132761471162,
                    8632934651105793861,
                    6631298214892334189,
                    1582556514881692819,
                ])),
            },
            c2: Fq2ZERO,
        },
        c1: Fq6 {
            c0: Fq2ZERO,
            c1: Fq2 {
                c0: Fq(BigInt([0; BASE_NUM_LIMBS])),
                c1: Fq(BigInt([
                    11671922859260663127,
                    11050707557586042878,
                    284884720401305268,
                    17749945728364010941,
                    8613774818643959860,
                    145621051382923523,
                ])),
            },
            c2: Fq2ZERO,
        },
    };

    // twist_beta_t_x^{-(p - 1)}
    const FROB_TWIST_X: BaseField = BaseField {
        c0: Fq(BigInt([0, 0, 0, 0, 0, 0])),
        c1: Fq(BigInt([
            9875771541238924739,
            3094855109658912213,
            5802897354862067244,
            11677019699073781796,
            1505592401347711080,
            1505729768134575418,
        ])),
    };
    // twist_beta_t_y^{-(p - 1)}
    const FROB_TWIST_Y: BaseField = BaseField {
        c0: Fq(BigInt([
            4480897313486445265,
            4797496051193971075,
            4046559893315008306,
            10569151167044009496,
            2123814803385151673,
            852749317591686856,
        ])),
        c1: Fq(BigInt([
            8921533702591418330,
            15859389534032789116,
            3389114680249073393,
            15116930867080254631,
            3288288975085550621,
            1021049300055853010,
        ])),
    };

    // for weil pairing
    fn untwist<C: Curve<FullExtensionField, ScalarField, FullExtensionField>>(
        p: &AffinePoint<BaseField, ScalarField, FullExtensionField, Self>,
    ) -> AffinePoint<FullExtensionField, ScalarField, FullExtensionField, C> {
        AffinePoint {
            x: Fq12 {
                c0: Fq6 {
                    c0: Self::TWIST_X_INV.c0.c0 * p.x,
                    c1: Self::TWIST_X_INV.c0.c1 * p.x,
                    c2: Self::TWIST_X_INV.c0.c2 * p.x,
                },
                c1: Fq6 {
                    c0: Self::TWIST_X_INV.c1.c0 * p.x,
                    c1: Self::TWIST_X_INV.c1.c1 * p.x,
                    c2: Self::TWIST_X_INV.c1.c2 * p.x,
                },
            },
            y: Fq12 {
                c0: Fq6 {
                    c0: Self::TWIST_Y_INV.c0.c0 * p.y,
                    c1: Self::TWIST_Y_INV.c0.c1 * p.y,
                    c2: Self::TWIST_Y_INV.c0.c2 * p.y,
                },
                c1: Fq6 {
                    c0: Self::TWIST_Y_INV.c1.c0 * p.y,
                    c1: Self::TWIST_Y_INV.c1.c1 * p.y,
                    c2: Self::TWIST_Y_INV.c1.c2 * p.y,
                },
            },
            _p1: PhantomData,
            _p2: PhantomData,
            _p3: PhantomData,
        }
    }

    fn is_on_curve(p: &AffinePoint<BaseField, ScalarField, FullExtensionField, Self>) -> bool {
        let (x, y) = (p.x, p.y);
        // let lft = y * y + Self::a1 * x * y + Self::a3 * y;
        // let rht = x * x * x + Self::a2 * x * x + Self::a4 * x + Self::a6;
        let lft = y.square();
        let rht = x.square() * x + Self::a4 * x + Self::a6;
        lft == rht
    }

    fn to_y(x: &BaseField) -> BaseField {
        unimplemented!()
    }

    fn neg(
        p: &AffinePoint<BaseField, ScalarField, FullExtensionField, Self>,
    ) -> AffinePoint<BaseField, ScalarField, FullExtensionField, Self> {
        AffinePoint {
            x: p.x,
            y: -p.y,
            _p1: PhantomData,
            _p2: PhantomData,
            _p3: PhantomData,
        }
    }

    fn is_negate(
        p1: &AffinePoint<BaseField, ScalarField, FullExtensionField, Self>,
        p2: &AffinePoint<BaseField, ScalarField, FullExtensionField, Self>,
    ) -> bool {
        // (p1.x == p2.x) && (p1.y + p2.y == -Self::a1 * p1.x - Self::a3)
        (p1.x == p2.x) && (p1.y + p2.y == BaseField::ZERO())
    }

    // referenced from P.270 of "handbook of elliptic and hyperelliptic curve cryptography"
    fn addition(
        p1: &AffinePoint<BaseField, ScalarField, FullExtensionField, Self>,
        p2: &AffinePoint<BaseField, ScalarField, FullExtensionField, Self>,
    ) -> AffinePoint<BaseField, ScalarField, FullExtensionField, Self> {
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
                x1.square() * BaseField::from_small_base_prime_field_str("3") + Self::a4,
                y1 * BaseField::from_small_base_prime_field_str("2"),
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
            _p3: PhantomData,
        }
    }

    // referenced from Algorithm 13.6 of "handbook of elliptic and hyperelliptic curve cryptography"
    fn scalar_mul(
        base: &AffinePoint<BaseField, ScalarField, FullExtensionField, Self>,
        scalar: &Fr<SCALAR_NUM_LIMBS>,
    ) -> AffinePoint<BaseField, ScalarField, FullExtensionField, Self> {
        if scalar.is_zero() {
            return AffinePoint::<BaseField, ScalarField, FullExtensionField, Self>::IDENTITY();
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
