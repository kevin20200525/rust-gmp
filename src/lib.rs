#![crate_name = "gmp"]

#![warn(deprecated)]
#![allow(non_camel_case_types)]

extern crate libc;
extern crate num_traits;

#[cfg(feature="serde_support")]
extern crate serde;

#[cfg(feature="serde_support")]
extern crate serde_derive;

#[cfg(feature="serde_support")]
extern crate serde_json;

mod ffi;
pub mod mpz;
pub mod mpq;
pub mod mpf;
pub mod rand;
pub mod sign;

#[cfg(test)]
mod test;
