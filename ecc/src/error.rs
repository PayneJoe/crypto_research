//! This module defines errors returned by the library.
use core::fmt::Debug;
use thiserror::Error;

/// Errors returned by Nova
#[derive(Clone, Debug, Eq, PartialEq, Error)]
pub enum MyError {
    /// group error
    #[error("group error")]
    GroupError,
    #[error("keccak error")]
    KeccakError,
}
