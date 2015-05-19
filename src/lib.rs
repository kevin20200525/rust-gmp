#![crate_name = "gmp"]

#![warn(deprecated)]
#![allow(non_camel_case_types)]

macro_rules! gen_overloads_inner {
    ($tr:ident, $meth:ident, $T:ident) => {
        impl<'a> $tr <&'a $T> for $T {
            type Output = $T;
            fn $meth(self, other: &'a $T) -> $T {
                (&self).$meth(other)
            }
        }
        impl<'a> $tr <$T> for &'a $T {
            type Output = $T;
            fn $meth(self, other: $T) -> $T {
                self.$meth(&other)
            }
        }
        impl $tr<$T> for $T {
            type Output = $T;
            fn $meth(self, other: $T) -> $T {
                (&self).$meth(&other)
            }
        }
    }
}

macro_rules! gen_overloads {
    ($T:ident) => {
        gen_overloads_inner!(Add, add, $T);
        gen_overloads_inner!(Sub, sub, $T);
        gen_overloads_inner!(Mul, mul, $T);
        gen_overloads_inner!(Div, div, $T);
        impl Neg for $T {
            type Output = $T;
            fn neg(self) -> $T {
                -(&self)
            }
        }
    }
}

pub mod mpz;
pub mod mpq;
pub mod mpf;
pub mod rand;

#[cfg(test)]
mod test;

#[test]
fn test_limb_size() {
    // We are assuming that the limb size is the same as the pointer size.
    assert_eq!(std::mem::size_of::<mp_limb_t>() * 8, __gmp_bits_per_limb as usize);
}
