pub mod addition;
pub mod division;
pub mod euclid_gcd;
pub mod modular;
pub mod multiplication;

use core::fmt::Debug;
#[derive(Debug, Clone, PartialEq)]
pub struct BigInteger {
    data: Vec<i32>,
    basis: i32,
}
