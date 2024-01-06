use std::ops::{Add, Div, Mul, Sub};
use std::str::FromStr;

//////////////////////////////////// Implementation of BigInteger Specially for Finite Field
#[derive(Debug, PartialEq)]
pub struct BI<const N: usize>(pub [u8; N]);

impl FromStr for BI<2> {
    type Err = ParseStrErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        unimplemented!()
    }
}

impl Add for BI<2> {
    type Output = (BI<2>, bool);
    fn add(self, other: BI<2>) -> Self::Output {
        unimplemented!()
    }
}

impl<'a, 'b> Add<&'b BI<2>> for &'a BI<2> {
    type Output = (BI<2>, bool);
    fn add(self, other: &'b BI<2>) -> Self::Output {
        unimplemented!()
    }
}

impl Sub for BI<2> {
    type Output = (BI<2>, bool);
    fn sub(self, other: BI<2>) -> Self::Output {
        unimplemented!()
    }
}

impl<'a, 'b> Sub<&'b BI<2>> for &'a BI<2> {
    type Output = (BI<2>, bool);
    fn sub(self, other: &'b BI<2>) -> Self::Output {
        unimplemented!()
    }
}

impl Mul for BI<2> {
    type Output = (BI<2>, bool);
    fn mul(self, other: BI<2>) -> Self::Output {
        unimplemented!()
    }
}

impl<'a, 'b> Mul<&'b BI<2>> for &'a BI<2> {
    type Output = (BI<2>, bool);
    fn mul(self, other: &'b BI<2>) -> Self::Output {
        unimplemented!()
    }
}

impl Div for BI<2> {
    type Output = (BI<2>, bool);
    fn div(self, other: BI<2>) -> Self::Output {
        unimplemented!()
    }
}

impl<'a, 'b> Div<&'b BI<2>> for &'a BI<2> {
    type Output = (BI<2>, bool);
    fn div(self, other: &'b BI<2>) -> Self::Output {
        unimplemented!()
    }
}

// field configurations including basic behavious on fields
pub trait Field<const N: usize>: FromStr + From<BI<N>> + Into<BI<N>> {
    // finite field modulus
    const MODULUS: BI<N>;
    // next power of ONE-WORD behind modulus, and then mod modulus, for Montgomery reduce
    const R: BI<N>;
    // inversion of least significant word of modulus, also convenient to Montgomery reduce
    const M0: u8;

    // reduce a long-long integer into the specified modulus
    fn reduce(u: &[u8]) -> Self;

    // basic arithmetic of finite field
    fn add_reduce(&self, other: &Self) -> Self;
    fn sub_reduce(&self, other: &Self) -> Self;
    fn mul_reduce(&self, other: &Self) -> Self;
    fn div_reduce(&self, other: &Self) -> Self;
    fn inv_reduce(&self) -> Self;
}

//////////////////////////////////////// custom finite field
#[derive(Debug)]
pub struct ParseStrErr;

#[derive(Debug, PartialEq)]
pub struct Foo<const T: usize>(BI<T>);

impl FromStr for Foo<2> {
    type Err = ParseStrErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        unimplemented!()
    }
}

impl From<BI<2>> for Foo<2> {
    fn from(value: BI<2>) -> Self {
        unimplemented!()
    }
}

impl Into<BI<2>> for Foo<2> {
    fn into(self) -> BI<2> {
        unimplemented!()
    }
}

impl Add for Foo<2> {
    type Output = Foo<2>;

    fn add(self, other: Self) -> Foo<2> {
        self.add_reduce(&other)
    }
}

impl<'a, 'b> Add<&'b Foo<2>> for &'a Foo<2> {
    type Output = Foo<2>;
    fn add(self, other: &'b Foo<2>) -> Foo<2> {
        self.add_reduce(other)
    }
}

impl Sub for Foo<2> {
    type Output = Foo<2>;

    fn sub(self, other: Self) -> Foo<2> {
        self.sub_reduce(&other)
    }
}

impl<'a, 'b> Sub<&'b Foo<2>> for &'a Foo<2> {
    type Output = Foo<2>;
    fn sub(self, other: &'b Foo<2>) -> Foo<2> {
        self.sub_reduce(other)
    }
}

impl Mul for Foo<2> {
    type Output = Foo<2>;

    fn mul(self, other: Self) -> Foo<2> {
        self.mul_reduce(&other)
    }
}

impl<'a, 'b> Mul<&'b Foo<2>> for &'a Foo<2> {
    type Output = Foo<2>;
    fn mul(self, other: &'b Foo<2>) -> Foo<2> {
        self.mul_reduce(other)
    }
}

impl Div for Foo<2> {
    type Output = Foo<2>;

    fn div(self, other: Self) -> Foo<2> {
        self.div_reduce(&other)
    }
}

impl<'a, 'b> Div<&'b Foo<2>> for &'a Foo<2> {
    type Output = Foo<2>;
    fn div(self, other: &'b Foo<2>) -> Foo<2> {
        self.div_reduce(other)
    }
}

// define custom finite field
impl Field<2> for Foo<2> {
    // fabricated configure parameters of custom finite field
    // these parameters need to determined at compile time
    const MODULUS: BI<2> = BI([5, 2]);
    const R: BI<2> = BI([138, 1]);
    const M0: u8 = 51_u8;

    fn reduce(u: &[u8]) -> Self {
        unimplemented!()
    }

    fn add_reduce(&self, other: &Self) -> Self {
        unimplemented!();
    }
    fn sub_reduce(&self, other: &Self) -> Self {
        unimplemented!();
    }
    fn mul_reduce(&self, other: &Self) -> Self {
        let s = self.0 .0.len();
        let t = vec![0_u8; s + 2];
        let carrier = 0_u8;

        unimplemented!();
    }
    fn div_reduce(&self, other: &Self) -> Self {
        unimplemented!();
    }
    fn inv_reduce(&self) -> Self {
        unimplemented!();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_config() {
        let result: BI<2> = Foo::<2>::from_str("259").unwrap().into();
        let actual = BI([3, 1]);
        assert_eq!(result, actual);
    }

    #[test]
    fn test_basic_ops() {
        let a = Foo::<2>::from_str("259").unwrap();
        let b = Foo::<2>::from_str("258").unwrap();
        let c = Foo::<2>::from_str("517").unwrap();
        assert_eq!(a + b, c);
    }
}
