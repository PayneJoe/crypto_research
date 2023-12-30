pub mod addition;
pub mod division;
pub mod gcd;
pub mod modulo;
pub mod multiplication;
pub mod squre_root;
pub mod substraction;

use core::fmt::Debug;
#[derive(Debug, Clone, PartialEq, Default)]
pub struct BigInteger {
    data: Vec<u8>,
    basis: u8,
}

impl BigInteger {
    fn is_zero(&self) -> bool {
        assert!(self.data.len() >= 1);
        self.data[0] == 0_u8
    }
}
