pub mod basic_ops;

use std::ops::{Add, AddAssign, Mul, Sub, SubAssign};

pub trait BasicOps<Rhs = Self, Output = Self>:
    Add<Rhs, Output = Output>
    + Sub<Rhs, Output = Output>
    + AddAssign<Rhs>
    + SubAssign<Rhs>
    + Mul<Rhs, Output = Output>
    + Default
    + PartialEq
{
}

impl<T, Rhs, Output> BasicOps<Rhs, Output> for T where
    T: Add<Rhs, Output = Output>
        + Sub<Rhs, Output = Output>
        + AddAssign<Rhs>
        + SubAssign<Rhs>
        + Mul<Rhs, Output = Output>
        + Default
        + PartialEq
{
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct BigInteger<T> {
    data: Vec<T>,
    sign: bool,
    basis: usize,
}

impl<T> Default for BigInteger<T>
where
    T: Default,
{
    fn default() -> Self {
        BigInteger {
            data: vec![T::default()],
            sign: false,
            basis: 1 as usize,
        }
    }
}

impl<T> BigInteger<T>
where
    T: Default + PartialEq + Clone,
{
    #[inline(always)]
    fn new(v: &[T], sign: bool, basis: usize) -> Self {
        let mut data = v.to_vec();
        data.reverse();
        BigInteger {
            data: data,
            sign: sign,
            basis: basis,
        }
    }

    #[inline(always)]
    fn size(&self) -> usize {
        self.data.len()
    }

    // highest digit is zero except the default ones
    #[inline(always)]
    fn is_leading_zero(&self) -> bool {
        (self.size() > 1) && (self.data[self.size() - 1] == T::default())
    }

    #[inline(always)]
    fn is_zero(&self) -> bool {
        let n = self.data.len();
        assert!(n != 0);
        (n == 1) && (self.data[0] == T::default())
    }

    #[inline(always)]
    fn is_positive(&self) -> bool {
        (!self.is_zero()) && (self.sign == false)
    }

    #[inline(always)]
    fn is_negative(&self) -> bool {
        (!self.is_zero()) && (self.sign == true)
    }

    fn strip_leading_zeros(&mut self) {
        if (!self.is_zero()) && (self.is_leading_zero()) {
            let mut end = self.size() - 1;
            while (end > 0) && (self.data[end] == T::default()) {
                end = end - 1;
            }
            self.data = self.data[..(end + 1)].to_vec();
        }
    }
}
