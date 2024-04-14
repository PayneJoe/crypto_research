/// This is a practical implementation for pallas curve over standard Weierstrass model
///
use crate::finite_field_arithmetic::bigint::BigInt;
use crate::finite_field_arithmetic::cycle_friendly::pallas::{fq::Fq, fr::Fr};
use crate::finite_field_arithmetic::traits::weierstrass_field::PrimeField;
use std::marker::PhantomData;

use crate::elliptic_curve_arithmetic::models::weierstrass_model::{AffinePoint, Curve};
use crate::utils;

// make sure WINDOW_SIZE < 8
const WINDOW_SIZE: usize = 6;
const NUM_LIMBS: usize = 4;
const WORD_SIZE: usize = 64;
type Word = u64;

type BigInteger = BigInt<NUM_LIMBS>;

// custom curve instance
#[derive(Debug, Clone, PartialEq, Eq, Copy, Hash)]
pub struct Pallas;

type ScalarField = Fr<NUM_LIMBS>;
type BaseField = Fq<NUM_LIMBS>;

// implementation of custom curve instance
impl Curve<BaseField, ScalarField> for Pallas {
    // y^2 = x^3 + 5
    // for the purpose of research, we fill these constant parameters through mannual computation
    // actually these constant parameters need to be determined dynamicly at compile-time
    const a1: Fq<NUM_LIMBS> = Fq(BigInt([0 as Word; NUM_LIMBS]));
    const a3: Fq<NUM_LIMBS> = Fq(BigInt([0 as Word; NUM_LIMBS]));
    const a2: Fq<NUM_LIMBS> = Fq(BigInt([0 as Word; NUM_LIMBS]));
    const a4: Fq<NUM_LIMBS> = Fq(BigInt([0 as Word; NUM_LIMBS]));
    const a6: Fq<NUM_LIMBS> = Fq(BigInt([
        11647819816328232941,
        8413468796752855795,
        18446744073709551613,
        4611686018427387903,
    ]));
    // faked identity representing the point at infinity
    const IDENTITY: AffinePoint<BaseField, ScalarField, Self> = AffinePoint {
        x: Fq(BigInt([0 as Word; NUM_LIMBS])),
        y: Fq(BigInt([0 as Word; NUM_LIMBS])),
        _p1: PhantomData,
        _p2: PhantomData,
    };
    // generator of this elliptic curve, g = (-1, 2)
    const GENERATOR: AffinePoint<BaseField, ScalarField, Self> = AffinePoint {
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
    // order of pallas curve is 28948022309329048855892746252171976963363056481941647379679742748393362948097
    // which is the modulus of vesta curve
    const ORDER: BigInteger = BigInt([
        10108024940646105089,
        2469829653919213789,
        0,
        4611686018427387904,
    ]);

    fn is_on_curve(p: &AffinePoint<BaseField, ScalarField, Self>) -> bool {
        let (x, y) = (p.x, p.y);
        let lft = y * y + Self::a1 * x * y + Self::a3 * y;
        let rht = x * x * x + Self::a2 * x * x + Self::a4 * x + Self::a6;
        lft == rht
    }

    fn to_y(x: &Fq<NUM_LIMBS>) -> Fq<NUM_LIMBS> {
        unimplemented!()
    }

    fn neg(
        p: &AffinePoint<BaseField, ScalarField, Self>,
    ) -> AffinePoint<BaseField, ScalarField, Self> {
        AffinePoint {
            x: p.x,
            y: -p.y - Self::a1 * p.x - Self::a3,
            _p1: PhantomData,
            _p2: PhantomData,
        }
    }

    fn is_negate(
        p1: &AffinePoint<BaseField, ScalarField, Self>,
        p2: &AffinePoint<BaseField, ScalarField, Self>,
    ) -> bool {
        (p1.x == p2.x) && (p1.y + p2.y == -Self::a1 * p1.x - Self::a3)
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
                x1 * x1 * (3 as Word) + Self::a2 * x1 * (2 as Word) + Self::a4 - Self::a1 * y1,
                y1 * (2 as Word) + Self::a1 * x1 + Self::a3,
            )
        } else {
            (y1 - y2, x1 - x2)
        };
        // println!(
        //     "denominator = {:?}, nominator = {:?}",
        //     denominator.rev_reduce(),
        //     nominator.rev_reduce()
        // );
        let lambda = denominator * nominator.inv();
        // FOR DEBUG
        assert!(nominator * nominator.inv() == BaseField::ONE());
        let x3 = lambda * lambda + Self::a1 * lambda - Self::a2 - x1 - x2;
        AffinePoint {
            x: x3,
            y: lambda * (x1 - x3) - y1 - Self::a1 * x3 - Self::a3,
            _p1: PhantomData,
            _p2: PhantomData,
        }
    }

    // referenced from Algorithm 13.6 of "handbook of elliptic and hyperelliptic curve cryptography"
    fn scalar_mul(
        base: &AffinePoint<BaseField, ScalarField, Self>,
        scalar: &Fr<NUM_LIMBS>,
    ) -> AffinePoint<BaseField, ScalarField, Self> {
        if scalar.is_zero() {
            return AffinePoint::<BaseField, ScalarField, Self>::IDENTITY();
        }
        // k < 8, make sure Word is big enough for store precomputated points
        // let k = 6;
        assert!(WINDOW_SIZE < 8);
        let scalar_limbs: BigInteger = (*scalar).into();
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

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::*;

    #[test]
    fn test_is_nonsingular() {
        assert_eq!(Pallas::is_nonsingular(), true)
    }

    #[test]
    fn test_generator() {
        let generator = Pallas::GENERATOR;
        let (x, y) = (
            "8846324870586583739111697172863888851520659981074829056268361825029048600336",
            "26279573376377669217484965006105184765343200454744742625721925061231198562202",
        );
        let generator_for_urs = AffinePoint::<BaseField, ScalarField, Pallas>::from_str(x, y);
        assert_eq!(generator.is_on_curve(), true);
        assert_eq!(generator_for_urs.is_on_curve(), true);
    }

    #[test]
    fn test_addition() {
        let ((ax, ay), (bx, by), (cx, cy)) = (
            (
                "27832571899082851396112614260331381199907587448938948538656515152087941036303",
                "24016138054826863494902451010045394661513815511955280984179509803134200597424",
            ),
            (
                "27832571899082851396112614260331381199907587448938948538656515152087941036303",
                "24016138054826863494902451010045394661513815511955280984179509803134200597424",
            ),
            (
                "1616285261656811539606192746606747756620678341484325945708180334022554109686",
                "25742719387172571165925478012183202284853392851062803962279894474527460895168",
            ), // (
               //     "9908836592418296524525327456052162783984680435925575327364933045222333063285",
               //     "8141182330133757207695661486101802544682001213708043141844125426486190597129",
               // ),
               // (
               //     "24705484178801360521248032223897679619079690363058394526611267317500158565519",
               //     "20116172318493854331373158539065837890266930911788035342907408169251780891917",
               // ),
               // (
               //     "28209655690986056210243374189450018186603428622824425839925213070114995880237",
               //     "12981965780070749014963957098443976813651938384094761603182376058552077049181",
               // ),
        );
        let (lft, rht, result) = (
            AffinePoint::<BaseField, ScalarField, Pallas>::from_str(ax, ay),
            AffinePoint::<BaseField, ScalarField, Pallas>::from_str(bx, by),
            AffinePoint::<BaseField, ScalarField, Pallas>::from_str(cx, cy),
        );
        let ret = lft + rht;
        assert_eq!(lft.is_on_curve(), true);
        assert_eq!(rht.is_on_curve(), true);
        assert_eq!(result.is_on_curve(), true);
        assert_eq!(ret.is_on_curve(), true);
        // println!("lft = {:?}, rht = {:?}, result = {:?}", lft, rht, result);
        assert_eq!(lft + rht, result);
    }

    #[test]
    fn test_scalar_mul() {
        let ((ax, ay), b, (cx, cy)) = (
            (
                "28948022309329048855892746252171976963363056481941560715954676764349967630336",
                "2",
            ),
            // "21712623616103011952848252560147696209404897876251782782637736724240743453506",
            // (
            //     "21897919324007031737251702070883943178876384330244155651599442226564320960121",
            //     "19992762491079946800453099425190032771033651151145773319238789454742728844037",
            // ),
            // "6918041010455432824",
            "19011490841973700698332973919617636205534468959166330882313323025205883376085",
            (
                "8796806641042381430351535528292925080765933408623233116962084588733109519666",
                "24000051983946361489643647773470857482924029485753766325300523193930768017055",
            ),
        );
        let (lft, rht, result) = (
            AffinePoint::<BaseField, ScalarField, Pallas>::from_str(ax, ay),
            ScalarField::from_str(b).unwrap(),
            AffinePoint::<BaseField, ScalarField, Pallas>::from_str(cx, cy),
        );
        let ret = lft * rht;

        assert_eq!(lft.is_on_curve(), true);
        assert_eq!(result.is_on_curve(), true);
        assert_eq!(ret.is_on_curve(), true);
        assert_eq!(ret, result);
    }

    #[test]
    fn test_negation() {
        let ((ax, ay), (cx, cy)) = (
            (
                "248525666866812560507066916710240758017531002361014375814245558108119597731",
                "3362287877906691235374814057740343136641788570196024186086103136382541086240",
            ),
            (
                "248525666866812560507066916710240758017531002361014375814245558108119597731",
                "25585734431422357620517932194431633826721267911745536529868573627967426544097",
            ),
        );
        let (lft, result) = (
            AffinePoint::<BaseField, ScalarField, Pallas>::from_str(ax, ay),
            AffinePoint::<BaseField, ScalarField, Pallas>::from_str(cx, cy),
        );
        assert_eq!(-lft, result);
    }

    #[test]
    fn test_substraction() {
        let ((ax, ay), (bx, by), (cx, cy)) = (
            (
                "19695523566563000590934567785358500667230186155941079306438240655728104147429",
                "9288424016392432538980327549656364256373834966733680910665980909883736936946",
            ),
            (
                "15966562881909561331802305473338881241257244978114549904895326206091026310095",
                "25573892053698188490698169163048014526514100808900586223876505459414980177669",
            ),
            (
                "20965631353460842626848471227656821850716571331093248540601546860419307242271",
                "21120711471727820202962953462315232896557399839414464009488072562085621651983",
            ),
        );
        let (lft, rht, result) = (
            AffinePoint::<BaseField, ScalarField, Pallas>::from_str(ax, ay),
            AffinePoint::<BaseField, ScalarField, Pallas>::from_str(bx, by),
            AffinePoint::<BaseField, ScalarField, Pallas>::from_str(cx, cy),
        );
        assert_eq!(lft - rht, result);
    }
}
