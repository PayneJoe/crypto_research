pub mod addition;
pub mod division;
pub mod euclid_gcd;
pub mod modulo;
pub mod multiplication;
pub mod substraction;

use core::fmt::Debug;
#[derive(Debug, Clone, PartialEq)]
pub struct BigInteger {
    data: Vec<i32>,
    basis: i32,
}
