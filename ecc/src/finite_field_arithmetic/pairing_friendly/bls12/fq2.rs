use crate::finite_field_arithmetic::bigint::BigInt;
use crate::finite_field_arithmetic::pairing_friendly::bls12::fq::Fq;
use crate::finite_field_arithmetic::pairing_friendly::field::{Field, PrimeField};
use crate::finite_field_arithmetic::pairing_friendly::quadratic_extension::{
    QuadraticExtension, QuadraticExtensionConfig,
};

const NUM_LIMBS: usize = 6;

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct Fq2Config;

impl QuadraticExtensionConfig<NUM_LIMBS> for Fq2Config {
    type BasePrimeField = Fq<NUM_LIMBS>;
    type BaseField = Fq<NUM_LIMBS>;
    type FrobCoeff = Fq<NUM_LIMBS>;

    // absolute degree of current (extension) field
    const DEGREE_OVER_BASE_PRIME_FIELD: usize = 2;

    // Fq[X] / X^2 - alpha, where alpha = -1 in Fq2/Fq over bls12
    const NON_QUADRATIC_RESIDUAL: Fq<NUM_LIMBS> = Fq(BigInt([
        4897101644811774638,
        3654671041462534141,
        569769440802610537,
        17053147383018470266,
        17227549637287919721,
        291242102765847046,
    ]));

    // coefficients of frobenius map over Fp2
    // alpha^{(p^d - 1)/2}
    // [1, -1] for bls12-381
    const FROBENIUS_COEFF_C1: &'static [Fq<NUM_LIMBS>] = &[
        Fq(BigInt([
            8505329371266088957,
            17002214543764226050,
            6865905132761471162,
            8632934651105793861,
            6631298214892334189,
            1582556514881692819,
        ])),
        Fq(BigInt([
            4897101644811774638,
            3654671041462534141,
            569769440802610537,
            17053147383018470266,
            17227549637287919721,
            291242102765847046,
        ])),
    ];

    fn multiply_frobenius_coeff(c: &mut Self::BaseField, power: usize) {
        *c = *c * Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
    }
}

// fq2 is an abstract type of quadratic extension over base prime field Fq
// F_q^2 = F_q[X] / X^2 - alpha
pub type Fq2 = QuadraticExtension<NUM_LIMBS, Fq2Config>;
impl Fq2 {}

mod tests {
    use std::str::FromStr;

    use super::*;
    #[test]
    fn test_addition() {
        let a_raw = vec![
            Fq::from_str("2600519326539830558472509080271623368201766255277844531820768498222883817097811802988058217645752353302064093108850").unwrap(), 
            Fq::from_str("3773918811240353028132708301698978168292380996566825694015677018101572110945786257480550114727558899711214737389297").unwrap()
        ];
        let b_raw = vec![
            Fq::from_str("1072443203470471287375910660069807232650870352399015460462646395708207704323743652585202041826373844903210482208294").unwrap(),
            Fq::from_str("1937629177128637437019849210056487608119787191192930938568869938274614678466303024924457365106197793056188608524635").unwrap(),
        ];
        let c_raw = vec![
            Fq::from_str("3672962530010301845848419740341430600852636607676859992283414893931091521421555455573260259472126198205274575317144").unwrap(),
            Fq::from_str("1709138433147323071734767686019561619855285367820748747252488820252155138921251417962319850704741028729509073354145").unwrap(),
        ];
        let (a, b, c) = (
            Fq2::from_base_prime_field_elems(a_raw).unwrap(),
            Fq2::from_base_prime_field_elems(b_raw).unwrap(),
            Fq2::from_base_prime_field_elems(c_raw).unwrap(),
        );
        assert_eq!(a + b, c);
    }

    #[test]
    fn test_substraction() {
        let a_raw = vec![
            Fq::from_str("2600519326539830558472509080271623368201766255277844531820768498222883817097811802988058217645752353302064093108850").unwrap(), 
            Fq::from_str("3773918811240353028132708301698978168292380996566825694015677018101572110945786257480550114727558899711214737389297").unwrap()
        ];
        let b_raw = vec![
            Fq::from_str("1072443203470471287375910660069807232650870352399015460462646395708207704323743652585202041826373844903210482208294").unwrap(),
            Fq::from_str("1937629177128637437019849210056487608119787191192930938568869938274614678466303024924457365106197793056188608524635").unwrap(),
        ];
        let c_raw = vec![
            Fq::from_str("1528076123069359271096598420201816135550895902878829071358122102514676112774068150402856175819378508398853610900556").unwrap(),
            Fq::from_str("1836289634111715591112859091642490560172593805373894755446807079826957432479483232556092749621361106655026128864662").unwrap(),
        ];
        let (a, b, c) = (
            Fq2::from_base_prime_field_elems(a_raw).unwrap(),
            Fq2::from_base_prime_field_elems(b_raw).unwrap(),
            Fq2::from_base_prime_field_elems(c_raw).unwrap(),
        );
        assert_eq!(a - b, c);
    }

    #[test]
    fn test_multiplication() {
        let a_raw = vec![
            Fq::from_str("2600519326539830558472509080271623368201766255277844531820768498222883817097811802988058217645752353302064093108850").unwrap(), 
            Fq::from_str("3773918811240353028132708301698978168292380996566825694015677018101572110945786257480550114727558899711214737389297").unwrap()
        ];
        let b_raw = vec![
            Fq::from_str("1072443203470471287375910660069807232650870352399015460462646395708207704323743652585202041826373844903210482208294").unwrap(),
            Fq::from_str("1937629177128637437019849210056487608119787191192930938568869938274614678466303024924457365106197793056188608524635").unwrap(),
        ];
        let c_raw = vec![
            Fq::from_str("556208749584661639957942654044599295690469944532141304145531378490586797153983612243811963241573572295298421905755").unwrap(),
            Fq::from_str("866544074811559545842722076162038564993428573250080593420943059316839286192636948251281627438375968628346956814291").unwrap(),
        ];
        let (a, b, c) = (
            Fq2::from_base_prime_field_elems(a_raw).unwrap(),
            Fq2::from_base_prime_field_elems(b_raw).unwrap(),
            Fq2::from_base_prime_field_elems(c_raw).unwrap(),
        );
        assert_eq!(a * b, c);
    }

    #[test]
    fn test_square() {
        let a_raw = vec![
            Fq::from_str("2600519326539830558472509080271623368201766255277844531820768498222883817097811802988058217645752353302064093108850").unwrap(), 
            Fq::from_str("3773918811240353028132708301698978168292380996566825694015677018101572110945786257480550114727558899711214737389297").unwrap()
        ];
        let c_raw = vec![
            Fq::from_str("2414268780359673411536291021229035896167802922632347653049330443285718833634818934749623933015380980921129376733478").unwrap(),
            Fq::from_str("2990207912459564482841439440112516037418304406499886569718075035418182988839528336665105475099950172019326796987650").unwrap(),
        ];
        let (a, c) = (
            Fq2::from_base_prime_field_elems(a_raw).unwrap(),
            Fq2::from_base_prime_field_elems(c_raw).unwrap(),
        );
        assert_eq!(a.square(), c);
    }

    #[test]
    fn test_square_root() {
        let a_raw = vec![
            Fq::from_str("1693922075044290799586349725248644656585209853597348488898111566020612993050193363830654421807029150461916786662716").unwrap(), 
            Fq::from_str("1307028274328876631654407493870275350503167406226236850937444676077287355413536344735451832959901439122302764659852").unwrap()
        ];
        let c_raw = vec![
            Fq::from_str("580410647848592055830372521352618176370581346071501298766544056530108483378644042998188696526419569327799015666059").unwrap(),
            Fq::from_str("3723377844579724351991614428907507999026362056113827237743761756407371509030296871660599888607067526098598631645156").unwrap(),
        ];
        let (a, c) = (
            Fq2::from_base_prime_field_elems(a_raw).unwrap(),
            Fq2::from_base_prime_field_elems(c_raw).unwrap(),
        );
        assert_eq!(a.sqrt().unwrap(), c);
    }

    #[test]
    fn test_division() {
        let a_raw = vec![
            Fq::from_str("1693922075044290799586349725248644656585209853597348488898111566020612993050193363830654421807029150461916786662716").unwrap(), 
            Fq::from_str("1307028274328876631654407493870275350503167406226236850937444676077287355413536344735451832959901439122302764659852").unwrap()
        ];
        let b_raw = vec![
            Fq::from_str("1072443203470471287375910660069807232650870352399015460462646395708207704323743652585202041826373844903210482208294").unwrap(),
            Fq::from_str("1937629177128637437019849210056487608119787191192930938568869938274614678466303024924457365106197793056188608524635").unwrap(),
        ];
        let c_raw = vec![
            Fq::from_str("578682175539542899404344330670147622954437015586306023176264550146069209468914364927419163097250145798034531111025").unwrap(),
            Fq::from_str("3224287095935887440332365375898050098092543654413624768616755955039132080286676720145658097465025463630610208936124").unwrap(),
        ];
        let (a, b, c) = (
            Fq2::from_base_prime_field_elems(a_raw).unwrap(),
            Fq2::from_base_prime_field_elems(b_raw).unwrap(),
            Fq2::from_base_prime_field_elems(c_raw).unwrap(),
        );
        assert_eq!(a / b, c);
    }
}
