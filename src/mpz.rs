use libc::{c_char, c_int, c_long, c_ulong, c_void, c_double, size_t};
use super::rand::gmp_randstate_t;
use std::convert::{From, Into};
use std::mem::{uninitialized,size_of};
use std::{fmt, hash};
use std::cmp::Ordering::{self, Greater, Less, Equal};
use std::str::FromStr;
use std::ops::{Div, Mul, Add, Sub, Neg, Shl, Shr, BitXor, BitAnd, BitOr, Rem};
use std::ffi::CString;

#[repr(C)]
pub struct mpz_struct {
    _mp_alloc: c_int,
    _mp_size: c_int,
    _mp_d: *mut c_void
}

pub type mp_limb_t = usize; // TODO: Find a way to use __gmp_bits_per_limb instead.
pub type mp_bitcnt_t = c_ulong;
pub type mpz_srcptr = *const mpz_struct;
pub type mpz_ptr = *mut mpz_struct;

#[link(name = "gmp")]
extern "C" {
    static __gmp_bits_per_limb: c_int;
    fn __gmpz_init(x: mpz_ptr);
    fn __gmpz_init2(x: mpz_ptr, n: mp_bitcnt_t);
    fn __gmpz_init_set(rop: mpz_ptr, op: mpz_srcptr);
    fn __gmpz_init_set_ui(rop: mpz_ptr, op: c_ulong);
    fn __gmpz_init_set_str(rop: mpz_ptr, s: *const c_char, base: c_int) -> c_int;
    fn __gmpz_clear(x: mpz_ptr);
    fn __gmpz_realloc2(x: mpz_ptr, n: mp_bitcnt_t);
    fn __gmpz_set(rop: mpz_ptr, op: mpz_srcptr);
    fn __gmpz_set_str(rop: mpz_ptr, s: *const c_char, base: c_int) -> c_int;
    fn __gmpz_get_str(s: *mut c_char, base: c_int, op: mpz_srcptr) -> *mut c_char;
    fn __gmpz_get_ui(op: mpz_srcptr) -> c_ulong;
    fn __gmpz_fits_ulong_p(op: mpz_srcptr) -> c_int;
    fn __gmpz_get_si(op: mpz_srcptr) -> c_ulong;
    fn __gmpz_get_d(op: mpz_srcptr) -> c_double;
    fn __gmpz_fits_slong_p(op: mpz_srcptr) -> c_long;
    fn __gmpz_sizeinbase(op: mpz_srcptr, base: c_int) -> size_t;
    fn __gmpz_cmp(op1: mpz_srcptr, op2: mpz_srcptr) -> c_int;
    fn __gmpz_cmp_ui(op1: mpz_srcptr, op2: c_ulong) -> c_int;
    fn __gmpz_add(rop: mpz_ptr, op1: mpz_srcptr, op2: mpz_srcptr);
    fn __gmpz_add_ui(rop: mpz_ptr, op1: mpz_srcptr, op2: c_ulong);
    fn __gmpz_sub(rop: mpz_ptr, op1: mpz_srcptr, op2: mpz_srcptr);
    fn __gmpz_sub_ui(rop: mpz_ptr, op1: mpz_srcptr, op2: c_ulong);
    fn __gmpz_ui_sub(rop: mpz_ptr, op1: c_ulong, op2: mpz_srcptr);
    fn __gmpz_mul(rop: mpz_ptr, op1: mpz_srcptr, op2: mpz_srcptr);
    fn __gmpz_mul_si(rop: mpz_ptr, op1: mpz_srcptr, op2: c_long);
    fn __gmpz_mul_2exp(rop: mpz_ptr, op1: mpz_srcptr, op2: mp_bitcnt_t);
    fn __gmpz_neg(rop: mpz_ptr, op: mpz_srcptr);
    fn __gmpz_abs(rop: mpz_ptr, op: mpz_srcptr);
    fn __gmpz_tdiv_q(q: mpz_ptr, n: mpz_srcptr, d: mpz_srcptr);
    fn __gmpz_tdiv_r(r: mpz_ptr, n: mpz_srcptr, d: mpz_srcptr);
    fn __gmpz_tdiv_q_ui(q: mpz_ptr, n: mpz_srcptr, d: c_ulong);
    fn __gmpz_tdiv_r_ui(r: mpz_ptr, n: mpz_srcptr, d: c_ulong);
    fn __gmpz_fdiv_q(q: mpz_ptr, n: mpz_srcptr, d: mpz_srcptr);
    fn __gmpz_fdiv_r(r: mpz_ptr, n: mpz_srcptr, d: mpz_srcptr);
    fn __gmpz_fdiv_q_2exp(q: mpz_ptr, n: mpz_srcptr, b: mp_bitcnt_t);
    fn __gmpz_mod(r: mpz_ptr, n: mpz_srcptr, d: mpz_srcptr);
    fn __gmpz_divisible_p(n: mpz_srcptr, d: mpz_srcptr) -> c_int;
    fn __gmpz_and(rop: mpz_ptr, op1: mpz_srcptr, op2: mpz_srcptr);
    fn __gmpz_ior(rop: mpz_ptr, op1: mpz_srcptr, op2: mpz_srcptr);
    fn __gmpz_xor(rop: mpz_ptr, op1: mpz_srcptr, op2: mpz_srcptr);
    fn __gmpz_com(rop: mpz_ptr, op: mpz_srcptr);
    fn __gmpz_popcount(op: mpz_srcptr) -> mp_bitcnt_t;
    fn __gmpz_pow_ui(rop: mpz_ptr, base: mpz_srcptr, exp: c_ulong);
    fn __gmpz_ui_pow_ui(rop: mpz_ptr, base: c_ulong, exp: c_ulong);
    fn __gmpz_powm(rop: mpz_ptr, base: mpz_srcptr, exp: mpz_srcptr, modulo: mpz_srcptr);
    fn __gmpz_hamdist(op1: mpz_srcptr, op2: mpz_srcptr) -> mp_bitcnt_t;
    fn __gmpz_setbit(rop: mpz_ptr, bit_index: mp_bitcnt_t);
    fn __gmpz_clrbit(rop: mpz_ptr, bit_index: mp_bitcnt_t);
    fn __gmpz_combit(rop: mpz_ptr, bit_index: mp_bitcnt_t);
    fn __gmpz_tstbit(rop: mpz_srcptr, bit_index: mp_bitcnt_t) -> c_int;
    fn __gmpz_nextprime(rop: mpz_ptr, op: mpz_srcptr);
    fn __gmpz_gcd(rop: mpz_ptr, op1: mpz_srcptr, op2: mpz_srcptr);
    fn __gmpz_gcdext(g: mpz_ptr, s: mpz_ptr, t: mpz_ptr, a: mpz_srcptr, b: mpz_srcptr);
    fn __gmpz_lcm(rop: mpz_ptr, op1: mpz_srcptr, op2: mpz_srcptr);
    fn __gmpz_invert(rop: mpz_ptr, op1: mpz_srcptr, op2: mpz_srcptr) -> c_int;
    fn __gmpz_import(rop: mpz_ptr, count: size_t, order: c_int, size: size_t,
                     endian: c_int, nails: size_t, op: *const c_void);
    fn __gmpz_export(rop: *mut c_void, countp: *mut size_t, order: c_int, size: size_t, 
                     endian: c_int, nails: size_t, op: mpz_srcptr);
    fn __gmpz_root(rop: mpz_ptr, op: mpz_srcptr, n: c_ulong) -> c_int;
    fn __gmpz_sqrt(rop: mpz_ptr, op: mpz_srcptr);
    fn __gmpz_millerrabin(n: mpz_srcptr, reps: c_int) -> c_int;
    fn __gmpz_urandomb(rop: mpz_ptr, state: gmp_randstate_t, n: mp_bitcnt_t);
    fn __gmpz_urandomm(rop: mpz_ptr, state: gmp_randstate_t, n: mpz_srcptr);
}

pub struct Mpz {
    pub mpz: mpz_struct,
}

unsafe impl Send for Mpz { }

impl Drop for Mpz {
    fn drop(&mut self) { unsafe { __gmpz_clear(&mut self.mpz) } }
}

impl Mpz {
    pub fn new() -> Mpz {
        unsafe {
            let mut mpz = uninitialized();
            __gmpz_init(&mut mpz);
            Mpz { mpz: mpz }
        }
    }

    pub fn new_reserve(n: usize) -> Mpz {
        unsafe {
            let mut mpz = uninitialized();
            __gmpz_init2(&mut mpz, n as c_ulong);
            Mpz { mpz: mpz }
        }
    }

    pub fn reserve(&mut self, n: usize) {
        if self.bit_length() < n {
            unsafe { __gmpz_realloc2(&mut self.mpz, n as c_ulong) }
        }
    }

    // TODO: fail on an invalid base
    // FIXME: Unfortunately it isn't currently possible to use the fmt::RadixFmt
    //        machinery for a custom type.
    pub fn to_str_radix(&self, base: u8) -> String {
        unsafe {
            // Extra two bytes are for possible minus sign and null terminator
            let len = __gmpz_sizeinbase(&self.mpz, base as c_int) as usize + 2;

            // Allocate and write into a raw *c_char of the correct length
            let mut vector: Vec<u8> = Vec::with_capacity(len);
            vector.set_len(len);

            __gmpz_get_str(vector.as_mut_ptr() as *mut _, base as c_int, &self.mpz);

            let mut first_nul = None;
            let mut index : usize = 0;
            for elem in &vector {
                if *elem == 0 {
                    first_nul = Some(index);
                    break;
                }
                index += 1;
            }
            let first_nul = first_nul.unwrap_or(len);

            vector.truncate(first_nul);
            match String::from_utf8(vector) {
                Ok(s)  => s,
                Err(_) => panic!("GMP returned invalid UTF-8!")
            }
        }
    }

    pub fn from_str_radix(s: &str, base: u8) -> Result<Mpz, ()> {
        unsafe {
            assert!(base == 0 || (base >= 2 && base <= 62));
            let mut mpz = uninitialized();
            let s = CString::new(s.to_string()).unwrap();
            let r = __gmpz_init_set_str(&mut mpz, s.as_ptr(), base as c_int);
            if r == 0 {
                Ok(Mpz { mpz: mpz })
            } else {
                __gmpz_clear(&mut mpz);
                Err(())
            }
        }
    }

    pub fn set(&mut self, other: &Mpz) {
        unsafe { __gmpz_set(&mut self.mpz, &other.mpz) }
    }

    // TODO: too easy to forget to check this return value - rename?
    pub fn set_from_str_radix(&mut self, s: &str, base: u8) -> bool {
        assert!(base == 0 || (base >= 2 && base <= 62));
        let s = CString::new(s.to_string()).unwrap();
        unsafe { __gmpz_set_str(&mut self.mpz, s.as_ptr(), base as c_int) == 0 }
    }

    pub fn bit_length(&self) -> usize {
        unsafe { __gmpz_sizeinbase(&self.mpz, 2) as usize }
    }

    pub fn compl(&self) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_com(&mut res.mpz, &self.mpz);
            res
        }
    }

    pub fn abs(&self) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_abs(&mut res.mpz, &self.mpz);
            res
        }
    }

    pub fn div_floor(&self, other: &Mpz) -> Mpz {
        unsafe {
            if other.is_zero() {
                panic!("divide by zero")
            }

            let mut res = Mpz::new();
            __gmpz_fdiv_q(&mut res.mpz, &self.mpz, &other.mpz);
            res
        }
    }

    pub fn mod_floor(&self, other: &Mpz) -> Mpz {
        unsafe {
            if other.is_zero() {
                panic!("divide by zero")
            }

            let mut res = Mpz::new();
            __gmpz_fdiv_r(&mut res.mpz, &self.mpz, &other.mpz);
            res
        }
    }

    pub fn nextprime(&self) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_nextprime(&mut res.mpz, &self.mpz);
            res
        }
    }

    pub fn gcd(&self, other: &Mpz) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_gcd(&mut res.mpz, &self.mpz, &other.mpz);
            res
        }
    }

    /// Given (a, b), return (g, s, t) such that g = gcd(a, b) = s*a + t*b.
    pub fn gcdext(&self, other: &Mpz) -> (Mpz, Mpz, Mpz) {
        unsafe {
            let mut g = Mpz::new();
            let mut s = Mpz::new();
            let mut t = Mpz::new();
            __gmpz_gcdext(&mut g.mpz, &mut s.mpz, &mut t.mpz,
                          &self.mpz, &other.mpz);
            (g, s, t)
        }
    }

    pub fn lcm(&self, other: &Mpz) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_lcm(&mut res.mpz, &self.mpz, &other.mpz);
            res
        }
    }

    pub fn is_multiple_of(&self, other: &Mpz) -> bool {
        unsafe {
            __gmpz_divisible_p(&self.mpz, &other.mpz) != 0
        }
    }

    #[inline]
    pub fn divides(&self, other: &Mpz) -> bool {
        other.is_multiple_of(self)
    }

    pub fn modulus(&self, modulo: &Mpz) -> Mpz {
        unsafe {
            if modulo.is_zero() {
                panic!("divide by zero")
            }

            let mut res = Mpz::new();
            __gmpz_mod(&mut res.mpz, &self.mpz, &modulo.mpz);
            res
        }
    }

    // TODO: handle a zero modulo
    pub fn invert(&self, modulo: &Mpz) -> Option<Mpz> {
        unsafe {
            let mut res = Mpz::new();
            if __gmpz_invert(&mut res.mpz, &self.mpz, &modulo.mpz) == 0 {
                None
            } else {
                Some(res)
            }
        }
    }

    pub fn popcount(&self) -> usize {
        unsafe { __gmpz_popcount(&self.mpz) as usize }
    }

    pub fn pow(&self, exp: u32) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_pow_ui(&mut res.mpz, &self.mpz, exp as c_ulong);
            res
        }
    }

    pub fn powm(&self, exp: &Mpz, modulus: &Mpz) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_powm(&mut res.mpz, &self.mpz, &exp.mpz, &modulus.mpz);
            res
        }
    }
    
    pub fn ui_pow_ui(x: u64, y: u64) -> Mpz {
    	unsafe {
    		let mut res = Mpz::new();
    		__gmpz_ui_pow_ui(&mut res.mpz, x as c_ulong, y as c_ulong);
    		res
    	}
    }

    pub fn hamdist(&self, other: &Mpz) -> usize {
        unsafe { __gmpz_hamdist(&self.mpz, &other.mpz) as usize }
    }

    pub fn setbit(&mut self, bit_index: usize) {
        unsafe { __gmpz_setbit(&mut self.mpz, bit_index as c_ulong) }
    }

    pub fn clrbit(&mut self, bit_index: usize) {
        unsafe { __gmpz_clrbit(&mut self.mpz, bit_index as c_ulong) }
    }

    pub fn combit(&mut self, bit_index: usize) {
        unsafe { __gmpz_combit(&mut self.mpz, bit_index as c_ulong) }
    }

    pub fn tstbit(&self, bit_index: usize) -> bool {
        unsafe { __gmpz_tstbit(&self.mpz, bit_index as c_ulong) == 1 }
    }

    pub fn root(&self, n: u64) -> Mpz {
        assert!(*self >= Mpz::zero());
        unsafe {
            let mut res = Mpz::new();
            let _perfect_root
                = match __gmpz_root(&mut res.mpz, &self.mpz, n as c_ulong) {
                    0 => false,
                    _ => true,
            };
            // TODO: consider returning `_perfect_root`
            res
        }
    }

    pub fn sqrt(&self) -> Mpz {
        assert!(*self >= Mpz::zero());
        unsafe {
            let mut res = Mpz::new();
            __gmpz_sqrt(&mut res.mpz, &self.mpz);
            res
        }
    }

    pub fn millerrabin(&self, reps: i32) -> i32 {
        unsafe {
            __gmpz_millerrabin(&self.mpz, reps as c_int)
        }
    }

    pub fn one() -> Mpz {
        unsafe {
            let mut mpz = uninitialized();
            __gmpz_init_set_ui(&mut mpz, 1);
            Mpz { mpz: mpz }
        }
    }

    pub fn zero() -> Mpz { Mpz::new() }

    pub fn is_zero(&self) -> bool {
        unsafe { __gmpz_cmp_ui(&self.mpz, 0) == 0 }
    }
}

impl Clone for Mpz {
    fn clone(&self) -> Mpz {
        unsafe {
            let mut mpz = uninitialized();
            __gmpz_init_set(&mut mpz, &self.mpz);
            Mpz { mpz: mpz }
        }
    }
}

impl Eq for Mpz { }

impl PartialEq for Mpz {
    fn eq(&self, other: &Mpz) -> bool {
        unsafe { __gmpz_cmp(&self.mpz, &other.mpz) == 0 }
    }
}

impl Ord for Mpz {
    fn cmp(&self, other: &Mpz) -> Ordering {
        let cmp = unsafe { __gmpz_cmp(&self.mpz, &other.mpz) };
        if cmp == 0 {
            Equal
        } else if cmp < 0 {
            Less
        } else {
            Greater
        }
    }
}

impl PartialOrd for Mpz {
    fn partial_cmp(&self, other: &Mpz) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<'a, 'b> Add<&'a Mpz> for &'b Mpz {
    type Output = Mpz;
    fn add(self, other: &Mpz) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_add(&mut res.mpz, &self.mpz, &other.mpz);
            res
        }
    }
}

impl Add<u64> for Mpz {
	type Output = Mpz;
	fn add(self, other: u64) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_add_ui(&mut res.mpz, &self.mpz, other as c_ulong);
            res
        }
	}
}

impl<'a> Add<u64> for &'a Mpz {
	type Output = Mpz;
	fn add(self, other: u64) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_add_ui(&mut res.mpz, &self.mpz, other as c_ulong);
            res
        }
	}
}

impl Add<Mpz> for u64 {
	type Output = Mpz;
	fn add(self, other: Mpz) -> Mpz {
		unsafe {
            let mut res = Mpz::new();
            __gmpz_add_ui(&mut res.mpz, &other.mpz, self as c_ulong);
            res
		}
	}
}

impl<'a> Add<&'a Mpz> for u64 {
	type Output = Mpz;
	fn add(self, other: &'a Mpz) -> Mpz {
		unsafe {
            let mut res = Mpz::new();
            __gmpz_add_ui(&mut res.mpz, &other.mpz, self as c_ulong);
            res
		}
	}
}

impl<'a, 'b> Sub<&'a Mpz> for &'b Mpz {
    type Output = Mpz;
    fn sub(self, other: &Mpz) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_sub(&mut res.mpz, &self.mpz, &other.mpz);
            res
        }
    }
}

impl Sub<u64> for Mpz {
	type Output = Mpz;
	fn sub(self, other: u64) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_sub_ui(&mut res.mpz, &self.mpz, other as c_ulong);
            res
        }
	}
}

impl<'a> Sub<u64> for &'a Mpz {
	type Output = Mpz;
	fn sub(self, other: u64) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_sub_ui(&mut res.mpz, &self.mpz, other as c_ulong);
            res
        }
	}
}

impl Sub<Mpz> for u64 {
	type Output = Mpz;
	fn sub(self, other: Mpz) -> Mpz {
		unsafe {
            let mut res = Mpz::new();
            __gmpz_ui_sub(&mut res.mpz, self as c_ulong, &other.mpz);
            res
		}
	}
}

impl<'a> Sub<&'a Mpz> for u64 {
	type Output = Mpz;
	fn sub(self, other: &'a Mpz) -> Mpz {
		unsafe {
            let mut res = Mpz::new();
            __gmpz_ui_sub(&mut res.mpz, self as c_ulong, &other.mpz);
            res
		}
	}
}

impl<'a, 'b> Mul<&'a Mpz> for &'b Mpz {
    type Output = Mpz;
    fn mul(self, other: &Mpz) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_mul(&mut res.mpz, &self.mpz, &other.mpz);
            res
        }
    }
}

impl Mul<i64> for Mpz {
	type Output = Mpz;
	fn mul(self, other: i64) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_mul_si(&mut res.mpz, &self.mpz, other as c_long);
            res
        }
	}
}

impl<'a> Mul<i64> for &'a Mpz {
	type Output = Mpz;
	fn mul(self, other: i64) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_mul_si(&mut res.mpz, &self.mpz, other as c_long);
            res
        }
	}
}

impl Mul<Mpz> for i64 {
	type Output = Mpz;
	fn mul(self, other: Mpz) -> Mpz {
		unsafe {
            let mut res = Mpz::new();
            __gmpz_mul_si(&mut res.mpz, &other.mpz, self as c_long);
            res
		}
	}
}

impl<'a> Mul<&'a Mpz> for i64 {
	type Output = Mpz;
	fn mul(self, other: &'a Mpz) -> Mpz {
		unsafe {
            let mut res = Mpz::new();
            __gmpz_mul_si(&mut res.mpz, &other.mpz, self as c_long);
            res
		}
	}
}

impl<'a, 'b> Div<&'a Mpz> for &'b Mpz {
    type Output = Mpz;
    fn div(self, other: &Mpz) -> Mpz {
        unsafe {
            if other.is_zero() {
                panic!("divide by zero")
            }

            let mut res = Mpz::new();
            __gmpz_tdiv_q(&mut res.mpz, &self.mpz, &other.mpz);
            res
        }
    }
}

impl Div<u64> for Mpz {
    type Output = Mpz;
    fn div(self, other: u64) -> Mpz {
        unsafe {
            if other == 0 {
                panic!("divide by zero")
            }

            let mut res = Mpz::new();
            __gmpz_tdiv_q_ui(&mut res.mpz, &self.mpz, other as c_ulong);
            res
        }
    }
}

impl<'a> Div<u64> for &'a Mpz {
    type Output = Mpz;
    fn div(self, other: u64) -> Mpz {
        unsafe {
            if other == 0 {
                panic!("divide by zero")
            }

            let mut res = Mpz::new();
            __gmpz_tdiv_q_ui(&mut res.mpz, &self.mpz, other as c_ulong);
            res
        }
    }
}

impl<'a, 'b> Rem<&'a Mpz> for &'b Mpz {
    type Output = Mpz;
    fn rem(self, other: &Mpz) -> Mpz {
        unsafe {
            if other.is_zero() {
                panic!("divide by zero")
            }

            let mut res = Mpz::new();
            __gmpz_tdiv_r(&mut res.mpz, &self.mpz, &other.mpz);
            res
        }
    }
}

impl Rem<u64> for Mpz {
    type Output = Mpz;
    fn rem(self, other: u64) -> Mpz {
        unsafe {
            if other == 0 {
                panic!("divide by zero")
            }

            let mut res = Mpz::new();
            __gmpz_tdiv_r_ui(&mut res.mpz, &self.mpz, other as c_ulong);
            res
        }
    }
}

impl<'a> Rem<u64> for &'a Mpz {
    type Output = Mpz;
    fn rem(self, other: u64) -> Mpz {
        unsafe {
            if other == 0 {
                panic!("divide by zero")
            }

            let mut res = Mpz::new();
            __gmpz_tdiv_r_ui(&mut res.mpz, &self.mpz, other as c_ulong);
            res
        }
    }
}

impl<'b> Neg for &'b Mpz {
    type Output = Mpz;
    fn neg(self) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_neg(&mut res.mpz, &self.mpz);
            res
        }
    }
}

impl<'b> Into<Option<i64>> for &'b Mpz {
    fn into(self) -> Option<i64> {
        unsafe {
            let negative = __gmpz_cmp_ui(&self.mpz, 0) < 0;
            let mut to_export = Mpz::new();

            if negative {
                __gmpz_com(&mut to_export.mpz, &self.mpz);
            } else {
                __gmpz_set(&mut to_export.mpz, &self.mpz);
            }

            if __gmpz_sizeinbase(&to_export.mpz, 2) <= 63 {
                let mut result : i64 = 0;
                __gmpz_export(&mut result as *mut i64 as *mut c_void, 0 as *mut size_t, -1, size_of::<i64>() as size_t, 0, 0, &to_export.mpz);
                if negative {
                    Some(result ^ -1i64)
                } else {
                    Some(result)
                }
            } else {
                return None;
            }
        }
    }
}

impl<'b> Into<Option<u64>> for &'b Mpz {
    fn into(self) -> Option<u64> {
        unsafe {
            if __gmpz_sizeinbase(&self.mpz, 2) <= 64 && __gmpz_cmp_ui(&self.mpz, 0) >= 0 {
                let mut result : u64 = 0;
                __gmpz_export(&mut result as *mut u64 as *mut c_void, 0 as *mut size_t, -1, size_of::<u64>() as size_t, 0, 0, &self.mpz);
                Some(result)
            } else {
                None
            }
        }
    }
}

impl<'a> Into<f64> for &'a Mpz {
    fn into(self) -> f64 {
        unsafe {
            __gmpz_get_d(&self.mpz) as f64
        }
    }
}

impl From<u64> for Mpz {
    fn from(other: u64) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_import(&mut res.mpz, 1, -1, size_of::<u64>() as size_t, 0, 0,
                          &other as *const u64 as *const c_void);
            res
        }
    }
}

impl From<i64> for Mpz {
    fn from(other: i64) -> Mpz {
        unsafe {
            let mut res = Mpz::new();

            if other.is_negative() {
                __gmpz_import(&mut res.mpz, 1, -1, size_of::<i64>() as size_t, 0, 0,
                              &(other ^ -1i64) as *const i64 as *const c_void);
                __gmpz_com(&mut res.mpz, &res.mpz);
            } else {
                __gmpz_import(&mut res.mpz, 1, -1, size_of::<i64>() as size_t, 0, 0,
                              &other as *const i64 as *const c_void);
            }
            res
        }
    }
}


impl<'a, 'b> BitAnd<&'a Mpz> for &'b Mpz {
    type Output = Mpz;
    fn bitand(self, other: &Mpz) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_and(&mut res.mpz, &self.mpz, &other.mpz);
            res
        }
    }
}

impl<'a, 'b> BitOr<&'a Mpz> for &'b Mpz {
    type Output = Mpz;
    fn bitor(self, other: &Mpz) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_ior(&mut res.mpz, &self.mpz, &other.mpz);
            res
        }
    }
}

impl<'a, 'b> BitXor<&'a Mpz> for &'b Mpz {
    type Output = Mpz;
    fn bitxor(self, other: &Mpz) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_xor(&mut res.mpz, &self.mpz, &other.mpz);
            res
        }
    }
}

impl<'b> Shl<usize> for &'b Mpz {
    type Output = Mpz;
    fn shl(self, other: usize) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_mul_2exp(&mut res.mpz, &self.mpz, other as c_ulong);
            res
        }
    }
}

impl<'b> Shr<usize> for &'b Mpz {
    type Output = Mpz;
    fn shr(self, other: usize) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_fdiv_q_2exp(&mut res.mpz, &self.mpz, other as c_ulong);
            res
        }
    }
}

impl Shl<usize> for Mpz {
    type Output = Mpz;
    fn shl(self, other: usize) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_mul_2exp(&mut res.mpz, &self.mpz, other as c_ulong);
            res
        }
    }
}

impl Shr<usize> for Mpz {
    type Output = Mpz;
    fn shr(self, other: usize) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_fdiv_q_2exp(&mut res.mpz, &self.mpz, other as c_ulong);
            res
        }
    }
}

impl FromStr for Mpz {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Mpz::from_str_radix(s, 10)
    }
}

impl fmt::Display for Mpz {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.to_str_radix(10))
    }
}

impl fmt::Debug for Mpz {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.to_str_radix(10))
    }
}

impl hash::Hash for Mpz {
    fn hash<S: hash::Hasher>(&self, state: &mut S) {
        unsafe {
            for i in 0..self.mpz._mp_size {
                let limb = self.mpz._mp_d as *const mp_limb_t;
                let limb = *(limb.offset(i as isize));
                limb.hash(state);
            }
        }
    }
}

gen_overloads!(Mpz);

gen_overloads_inner!(BitXor, bitxor, Mpz);
gen_overloads_inner!(BitAnd, bitand, Mpz);
gen_overloads_inner!(BitOr, bitor, Mpz);
gen_overloads_inner!(Rem, rem, Mpz);