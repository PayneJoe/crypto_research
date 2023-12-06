/// standard lib
use core::fmt::Debug;
use core::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

// third-party lib
use ff::{PrimeField, PrimeFieldBits};
use serde::{Deserialize, Serialize};

// custom defined lib
use crate::error::MyError;

///////////////////////////////////////////
/// basic marker trait for Add/AddAssign/Sub/SubAssign
/// 1. a + b
/// 2. a += b
/// 3. a - b
/// 4. a -= b
pub trait GroupOps<Rhs = Self, Output = Self>:
    Add<Rhs, Output = Output> + Sub<Rhs, Output = Output> + AddAssign<Rhs> + SubAssign<Rhs>
{
}

/// marker trait for Mul/MulAssign
/// 1. a * b
/// 2. a *= b
pub trait ScalarMul<Rhs, Output = Self>: Mul<Rhs, Output = Output> + MulAssign<Rhs> {}

////////////////////////////////////////////  impl ops trait for T
// impl<T, Rhs, Output> GroupOps<Rhs, Output> for T where
//   T: Add<Rhs, Output = Output> + Sub<Rhs, Output = Output> + AddAssign<Rhs> + SubAssign<Rhs>
// {
// }
// impl<T, Rhs, Output> ScalarMul<Rhs, Output> for T where T: Mul<Rhs, Output = Output> + MulAssign<Rhs>
// {}

//////////////////////////////////////////  owned ops trait
/// owned marker trait
// pub trait GroupOpsOwned<Rhs = Self, Output = Self>: for<'r> GroupOps<&'r Rhs, Output> {}
// pub trait ScalarMulOwned<Rhs, Output = Self>: for<'r> ScalarMul<&'r Rhs, Output> {}

///////////////////////////////////////////  impl owned ops trait for T
// impl<T, Rhs, Output> GroupOpsOwned<Rhs, Output> for T where T: for<'r> GroupOps<&'r Rhs, Output> {}
// impl<T, Rhs, Output> ScalarMulOwned<Rhs, Output> for T where T: for<'r> ScalarMul<&'r Rhs, Output> {}

////////////////////////////////////////// transcript trait
pub trait TranscriptReprTrait<G: Group> {
    /// returns a byte representation of self to be added to the transcript
    fn to_transcript_bytes(&self) -> Vec<u8>;
}

pub trait TranscriptEngineTrait<G: Group> {
    /// initializes the transcript
    fn new(label: &'static [u8]) -> Self;

    /// returns a scalar element of the group as a challenge
    fn squeeze(&mut self, label: &'static [u8]) -> Result<G::Scalar, MyError>;

    /// absorbs any type that implements TranscriptReprTrait under a label
    fn absorb<T: TranscriptReprTrait<G>>(&mut self, label: &'static [u8], o: &T);

    /// adds a domain separator
    fn dom_sep(&mut self, bytes: &'static [u8]);
}

/////////////////////////////////////////// compressed point
/// output x-coordinate along with the sign y-coordinate provided a AffinePoint (x, y)
pub trait CompressedGroup:
    Copy + Clone + Debug + Eq + Sized + Serialize + for<'de> Deserialize<'de>
{
    /// A type that holds the decompressed version of the compressed group element
    type GroupElement: Group + Serialize + for<'de> Deserialize<'de>;

    /// Decompresses the compressed group element
    fn decompress(&self) -> Option<Self::GroupElement>;
}

////////////////////////////////////////// uncompressed point
pub trait Group:
    Copy
    + Clone
    + Debug
    + Eq
    + Sized
    + GroupOps
    + ScalarMul<Self::Scalar>
    + Serialize
    + for<'de> Deserialize<'de>
{
    // field type: inner type representing scalar field of group element (#E(Fp))
    type Scalar: PrimeField + PrimeFieldBits + Serialize + for<'de> Deserialize<'de>;

    // field type: inner type representing base field of group element (Fp)
    type Base: PrimeField + PrimeFieldBits + Serialize + for<'de> Deserialize<'de>;

    // point type: inner type representing compressed group element (E(Fp))
    type CompressedGroupElement: CompressedGroup<GroupElement = Self>;

    // point type: inner type representing preprocessed group element (E(Fp))
    type PreprocessedGroupElement: Clone + Debug + Serialize + for<'de> Deserialize<'de>;

    // hash type: inner type representing hasher for the purpose of folding factor (#E(Fp))
    // type RO: ROTrait<Self::Base, Self::Scalar> + Serialize + for<'de> Deserialize<'de>;

    // hash type: inner type representing hasher for Fiat-Shamir transcript
    type TE: TranscriptEngineTrait<Self>;

    // pcs type: inner type representing commitment scheme
    // type CE: CommitmentEngineTrait<Self> + Serialize + for<'de> Deserialize<'de>;
}
