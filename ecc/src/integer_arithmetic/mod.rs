pub mod addition;
pub mod division;
pub mod gcd;
pub mod inversion;
pub mod modulo;
pub mod multiplication;
pub mod substraction;

use core::fmt::Debug;
#[derive(Debug, Clone, PartialEq)]
pub struct BigInteger {
    data: Vec<u8>,
    basis: u8,
}
