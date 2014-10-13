#![crate_name = "gmp"]

#![comment = "gmp bindings"]
#![license = "MIT"]
#![crate_type = "lib"]

#![feature(globs)]
#![allow(non_camel_case_types)]

extern crate libc;

use libc::{c_char, c_double, c_int, c_long, c_ulong, c_void, size_t};
use std::num::{One, Zero, ToStrRadix};
use std::mem::{uninitialized,size_of};
use std::{cmp, fmt};
use std::from_str::FromStr;

#[cfg(test)]
mod test;

#[repr(C)]
struct mpz_struct {
    _mp_alloc: c_int,
    _mp_size: c_int,
    _mp_d: *mut c_void
}

#[repr(C)]
struct mpq_struct {
    _mp_num: mpz_struct,
    _mp_den: mpz_struct
}

type mp_exp_t = c_long;

#[repr(C)]
struct mpf_struct {
    _mp_prec: c_int,
    _mp_size: c_int,
    _mp_exp: mp_exp_t,
    _mp_d: *mut c_void
}

#[repr(C)]
struct gmp_randstate_struct {
    _mp_seed: mpz_struct,
    _mp_alg: c_int,
    _mp_algdata: *const c_void
}

type mp_bitcnt_t = c_ulong;
type mpz_srcptr = *const mpz_struct;
type mpz_ptr = *mut mpz_struct;
type mpq_srcptr = *const mpq_struct;
type mpq_ptr = *mut mpq_struct;
type mpf_srcptr = *const mpf_struct;
type mpf_ptr = *mut mpf_struct;
type gmp_randstate_t = *mut gmp_randstate_struct;

#[link(name = "gmp")]
extern "C" {
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
    fn __gmpz_sizeinbase(op: mpz_srcptr, base: c_int) -> size_t;
    fn __gmpz_cmp(op1: mpz_srcptr, op2: mpz_srcptr) -> c_int;
    fn __gmpz_cmp_ui(op1: mpz_srcptr, op2: c_ulong) -> c_int;
    fn __gmpz_add(rop: mpz_ptr, op1: mpz_srcptr, op2: mpz_srcptr);
    fn __gmpz_sub(rop: mpz_ptr, op1: mpz_srcptr, op2: mpz_srcptr);
    fn __gmpz_mul(rop: mpz_ptr, op1: mpz_srcptr, op2: mpz_srcptr);
    fn __gmpz_mul_2exp(rop: mpz_ptr, op1: mpz_srcptr, op2: mp_bitcnt_t);
    fn __gmpz_neg(rop: mpz_ptr, op: mpz_srcptr);
    fn __gmpz_abs(rop: mpz_ptr, op: mpz_srcptr);
    fn __gmpz_tdiv_q(q: mpz_ptr, n: mpz_srcptr, d: mpz_srcptr);
    fn __gmpz_tdiv_r(r: mpz_ptr, n: mpz_srcptr, d: mpz_srcptr);
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
    fn __gmpz_hamdist(op1: mpz_srcptr, op2: mpz_srcptr) -> mp_bitcnt_t;
    fn __gmpz_setbit(rop: mpz_ptr, bit_index: mp_bitcnt_t);
    fn __gmpz_clrbit(rop: mpz_ptr, bit_index: mp_bitcnt_t);
    fn __gmpz_combit(rop: mpz_ptr, bit_index: mp_bitcnt_t);
    fn __gmpz_tstbit(rop: mpz_srcptr, bit_index: mp_bitcnt_t) -> c_int;
    fn __gmpz_gcd(rop: mpz_ptr, op1: mpz_srcptr, op2: mpz_srcptr);
    fn __gmpz_lcm(rop: mpz_ptr, op1: mpz_srcptr, op2: mpz_srcptr);
    fn __gmpz_invert(rop: mpz_ptr, op1: mpz_srcptr, op2: mpz_srcptr) -> c_int;
    fn __gmpz_import(rop: mpz_ptr, count: size_t, order: c_int, size: size_t,
                     endian: c_int, nails: size_t, op: *const c_void);
    fn __gmp_randinit_default(state: gmp_randstate_t);
    fn __gmp_randinit_mt(state: gmp_randstate_t);
    fn __gmp_randinit_lc_2exp(state: gmp_randstate_t, a: mpz_srcptr, c: c_ulong, m2exp: mp_bitcnt_t);
    fn __gmp_randinit_lc_2exp_size(state: gmp_randstate_t, size: mp_bitcnt_t);
    fn __gmp_randinit_set(state: gmp_randstate_t, op: *const gmp_randstate_struct);
    fn __gmp_randclear(state: gmp_randstate_t);
    fn __gmp_randseed(state: gmp_randstate_t, seed: mpz_srcptr);
    fn __gmp_randseed_ui(state: gmp_randstate_t, seed: c_ulong);
    fn __gmpz_urandomb(rop: mpz_ptr, state: gmp_randstate_t, n: mp_bitcnt_t);
    fn __gmpz_urandomm(rop: mpz_ptr, state: gmp_randstate_t, n: mpz_srcptr);
    fn __gmpq_init(x: mpq_ptr);
    fn __gmpq_clear(x: mpq_ptr);
    fn __gmpq_set(rop: mpq_ptr, op: mpq_srcptr);
    fn __gmpq_set_z(rop: mpq_ptr, op: mpz_srcptr);
    fn __gmpq_set_ui(rop: mpq_ptr, op1: c_ulong, op2: c_ulong);
    fn __gmpq_set_d(rop: mpq_ptr, op: c_double);
    fn __gmpq_set_f(rop: mpq_ptr, op: mpf_srcptr);
    fn __gmpq_cmp(op1: mpq_srcptr, op2: mpq_srcptr) -> c_int;
    fn __gmpq_cmp_ui(op1: mpq_srcptr, num2: c_ulong, den2: c_ulong) -> c_int;
    fn __gmpq_equal(op1: mpq_srcptr, op2: mpq_srcptr) -> c_int;
    fn __gmpq_add(sum: mpq_ptr, addend1: mpq_srcptr, addend2: mpq_srcptr);
    fn __gmpq_sub(difference: mpq_ptr, minuend: mpq_srcptr, subtrahend: mpq_srcptr);
    fn __gmpq_mul(product: mpq_ptr, multiplier: mpq_srcptr, multiplicand: mpq_srcptr);
    fn __gmpq_div(product: mpq_ptr, multiplier: mpq_srcptr, multiplicand: mpq_srcptr);
    fn __gmpq_neg(negated_operand: mpq_ptr, operand: mpq_srcptr);
    fn __gmpq_abs(rop: mpq_ptr, op: mpq_srcptr);
    fn __gmpq_inv(inverted_number: mpq_ptr, number: mpq_srcptr);
    fn __gmpq_get_num(numerator: mpz_ptr, rational: mpq_srcptr);
    fn __gmpq_get_den(denominator: mpz_ptr, rational: mpq_srcptr);
    fn __gmpf_init2(x: mpf_ptr, prec: mp_bitcnt_t);
    fn __gmpf_init_set(rop: mpf_ptr, op: mpf_srcptr);
    fn __gmpf_clear(x: mpf_ptr);
    fn __gmpf_get_prec(op: mpf_srcptr) -> mp_bitcnt_t;
    fn __gmpf_set_prec(rop: mpf_ptr, prec: mp_bitcnt_t);
    fn __gmpf_set(rop: mpf_ptr, op: mpf_srcptr);
    fn __gmpf_cmp(op1: mpf_srcptr, op2: mpf_srcptr) -> c_int;
    fn __gmpf_cmp_d(op1: mpf_srcptr, op2: c_double) -> c_int;
    fn __gmpf_cmp_ui(op1: mpf_srcptr, op2: c_ulong) -> c_int;
    fn __gmpf_reldiff(rop: mpf_ptr, op1: mpf_srcptr, op2: mpf_srcptr);
    fn __gmpf_add(rop: mpf_ptr, op1: mpf_srcptr, op2: mpf_srcptr);
    fn __gmpf_sub(rop: mpf_ptr, op1: mpf_srcptr, op2: mpf_srcptr);
    fn __gmpf_mul(rop: mpf_ptr, op1: mpf_srcptr, op2: mpf_srcptr);
    fn __gmpf_div(rop: mpf_ptr, op1: mpf_srcptr, op2: mpf_srcptr);
    fn __gmpf_neg(rop: mpf_ptr, op: mpf_srcptr);
    fn __gmpf_abs(rop: mpf_ptr, op: mpf_srcptr);
    fn __gmpf_ceil(rop: mpf_ptr, op: mpf_srcptr);
    fn __gmpf_floor(rop: mpf_ptr, op: mpf_srcptr);
    fn __gmpf_trunc(rop: mpf_ptr, op: mpf_srcptr);
}

pub struct Mpz {
    mpz: mpz_struct,
}

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

    pub fn new_reserve(n: c_ulong) -> Mpz {
        unsafe {
            let mut mpz = uninitialized();
            __gmpz_init2(&mut mpz, n);
            Mpz { mpz: mpz }
        }
    }

    pub fn reserve(&mut self, n: c_ulong) {
        if (self.bit_length() as c_ulong) < n {
            unsafe { __gmpz_realloc2(&mut self.mpz, n) }
        }
    }

    pub fn from_str_radix(s: &str, base: uint) -> Option<Mpz> {
        unsafe {
            assert!(base == 0 || (base >= 2 && base <= 62));
            let mut mpz = uninitialized();
            let r = s.with_c_str(|s| {
                __gmpz_init_set_str(&mut mpz, s, base as c_int)
            });
            if r == 0 {
                Some(Mpz { mpz: mpz })
            } else {
                __gmpz_clear(&mut mpz);
                None
            }
        }
    }

    pub fn set(&mut self, other: &Mpz) {
        unsafe { __gmpz_set(&mut self.mpz, &other.mpz) }
    }

    // TODO: too easy to forget to check this return value - rename?
    pub fn set_from_str_radix(&mut self, s: &str, base: uint) -> bool {
        assert!(base == 0 || (base >= 2 && base <= 62));
        s.with_c_str(|s| {
            unsafe { __gmpz_set_str(&mut self.mpz, s, base as c_int) == 0 }
        })
    }

    pub fn bit_length(&self) -> uint {
        unsafe { __gmpz_sizeinbase(&self.mpz, 2) as uint }
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
                fail!("divide by zero")
            }

            let mut res = Mpz::new();
            __gmpz_fdiv_q(&mut res.mpz, &self.mpz, &other.mpz);
            res
        }
    }

    pub fn mod_floor(&self, other: &Mpz) -> Mpz {
        unsafe {
            if other.is_zero() {
                fail!("divide by zero")
            }

            let mut res = Mpz::new();
            __gmpz_fdiv_r(&mut res.mpz, &self.mpz, &other.mpz);
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
                fail!("divide by zero")
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

    pub fn popcount(&self) -> uint {
        unsafe { __gmpz_popcount(&self.mpz) as uint }
    }

    pub fn hamdist(&self, other: &Mpz) -> uint {
        unsafe { __gmpz_hamdist(&self.mpz, &other.mpz) as uint }
    }

    pub fn setbit(&mut self, bit_index: c_ulong) {
        unsafe { __gmpz_setbit(&mut self.mpz, bit_index) }
    }

    pub fn clrbit(&mut self, bit_index: c_ulong) {
        unsafe { __gmpz_clrbit(&mut self.mpz, bit_index) }
    }

    pub fn combit(&mut self, bit_index: c_ulong) {
        unsafe { __gmpz_combit(&mut self.mpz, bit_index) }
    }

    pub fn tstbit(&self, bit_index: c_ulong) -> bool {
        unsafe { __gmpz_tstbit(&self.mpz, bit_index) == 1 }
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

impl Eq for Mpz {

}

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

impl Add<Mpz, Mpz> for Mpz {
    fn add(&self, other: &Mpz) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_add(&mut res.mpz, &self.mpz, &other.mpz);
            res
        }
    }
}

impl Sub<Mpz, Mpz> for Mpz {
    fn sub(&self, other: &Mpz) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_sub(&mut res.mpz, &self.mpz, &other.mpz);
            res
        }
    }
}

impl Mul<Mpz, Mpz> for Mpz {
    fn mul(&self, other: &Mpz) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_mul(&mut res.mpz, &self.mpz, &other.mpz);
            res
        }
    }
}

impl Div<Mpz, Mpz> for Mpz {
    fn div(&self, other: &Mpz) -> Mpz {
        unsafe {
            if other.is_zero() {
                fail!("divide by zero")
            }

            let mut res = Mpz::new();
            __gmpz_tdiv_q(&mut res.mpz, &self.mpz, &other.mpz);
            res
        }
    }
}

impl Rem<Mpz, Mpz> for Mpz {
    fn rem(&self, other: &Mpz) -> Mpz {
        unsafe {
            if self.is_zero() {
                fail!("divide by zero")
            }

            let mut res = Mpz::new();
            __gmpz_tdiv_r(&mut res.mpz, &self.mpz, &other.mpz);
            res
        }
    }
}

impl Neg<Mpz> for Mpz {
    fn neg(&self) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_neg(&mut res.mpz, &self.mpz);
            res
        }
    }
}

impl ToPrimitive for Mpz {
    fn to_i64(&self) -> Option<i64> {
        fail!("not implemented")
    }
    fn to_u64(&self) -> Option<u64> {
        fail!("not implemented")
    }
}

impl FromPrimitive for Mpz {
    fn from_u64(other: u64) -> Option<Mpz> {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_import(&mut res.mpz, 1, 1, size_of::<u64>() as size_t, 0, 0,
                          &other as *const u64 as *const c_void);
            Some(res)
        }
    }
    fn from_i64(other: i64) -> Option<Mpz> {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_import(&mut res.mpz, 1, 1, size_of::<i64>() as size_t, 0, 0,
                          &other.abs() as *const i64 as *const c_void);
            if other.is_negative() {
                __gmpz_neg(&mut res.mpz, &res.mpz)
            }
            Some(res)
        }
    }
}

impl One for Mpz {
    fn one() -> Mpz {
        unsafe {
            let mut mpz = uninitialized();
            __gmpz_init_set_ui(&mut mpz, 1);
            Mpz { mpz: mpz }
        }
    }
}

impl Zero for Mpz {
    fn zero() -> Mpz { Mpz::new() }
    fn is_zero(&self) -> bool {
        unsafe { __gmpz_cmp_ui(&self.mpz, 0) == 0 }
    }
}

impl BitAnd<Mpz, Mpz> for Mpz {
    fn bitand(&self, other: &Mpz) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_and(&mut res.mpz, &self.mpz, &other.mpz);
            res
        }
    }
}

impl BitOr<Mpz, Mpz> for Mpz {
    fn bitor(&self, other: &Mpz) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_ior(&mut res.mpz, &self.mpz, &other.mpz);
            res
        }
    }
}

impl BitXor<Mpz, Mpz> for Mpz {
    fn bitxor(&self, other: &Mpz) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_xor(&mut res.mpz, &self.mpz, &other.mpz);
            res
        }
    }
}

impl Shl<c_ulong, Mpz> for Mpz {
    fn shl(&self, other: &c_ulong) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_mul_2exp(&mut res.mpz, &self.mpz, *other);
            res
        }
    }
}

impl Shr<c_ulong, Mpz> for Mpz {
    fn shr(&self, other: &c_ulong) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_fdiv_q_2exp(&mut res.mpz, &self.mpz, *other);
            res
        }
    }
}

impl FromStr for Mpz {
    fn from_str(s: &str) -> Option<Mpz> {
        Mpz::from_str_radix(s, 10)
    }
}

impl ToStrRadix for Mpz {
    // TODO: fail on an invalid base
    fn to_str_radix(&self, base: uint) -> String {
        unsafe {
            // Extra two bytes are for possible minus sign and null terminator
            let len = __gmpz_sizeinbase(&self.mpz, base as c_int) as uint + 2;

            // Allocate and write into a raw *c_char of the correct length
            let mut vector: Vec<u8> = Vec::with_capacity(len);
            vector.set_len(len);

            let mut cstr = vector.as_slice().to_c_str_unchecked();

            __gmpz_get_str(cstr.as_mut_ptr(), base as c_int, &self.mpz);

            match cstr.as_str() {
                Some(slice) => slice.to_string(),
                None        => fail!("GMP returned invalid UTF-8!")
            }
        }
    }
}

impl fmt::Show for Mpz {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.to_str_radix(10))
    }
}


pub struct RandState {
    state: gmp_randstate_struct,
}

impl Drop for RandState {
    fn drop(&mut self) {
        unsafe { __gmp_randclear(&mut self.state) }
    }
}

impl RandState {
    pub fn new() -> RandState {
        unsafe {
            let mut state: gmp_randstate_struct = uninitialized();
            __gmp_randinit_default(&mut state);
            RandState { state: state }
        }
    }

    pub fn new_mt() -> RandState {
        unsafe {
            let mut state: gmp_randstate_struct = uninitialized();
            __gmp_randinit_mt(&mut state);
            RandState { state: state }
        }
    }

    pub fn new_lc_2exp(a: Mpz, c: c_ulong, m2exp: c_ulong) -> RandState {
        unsafe {
            let mut state: gmp_randstate_struct = uninitialized();
            __gmp_randinit_lc_2exp(&mut state, &a.mpz, c, m2exp);
            RandState { state: state }
        }
    }

    pub fn new_lc_2exp_size(size: c_ulong) -> RandState {
        unsafe {
            let mut state: gmp_randstate_struct = uninitialized();
            __gmp_randinit_lc_2exp_size(&mut state, size);
            RandState { state: state }
        }
    }

    pub fn seed(&mut self, seed: Mpz) {
        unsafe { __gmp_randseed(&mut self.state, &seed.mpz) }
    }

    pub fn seed_ui(&mut self, seed: c_ulong) {
        unsafe { __gmp_randseed_ui(&mut self.state, seed) }
    }

    /// Generate a uniform random integer in the range 0 to n-1, inclusive
    pub fn urandom(&mut self, n: &Mpz) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_urandomm(&mut res.mpz, &mut self.state, &n.mpz);
            res
        }
    }

    /// Generate a uniformly distributed random integer in the range 0 to 2^n−1, inclusive.
    pub fn urandom_2exp(&mut self, n: c_ulong) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpz_urandomb(&mut res.mpz, &mut self.state, n);
            res
        }
    }
}

impl Clone for RandState {
    fn clone(&self) -> RandState {
        unsafe {
            let mut state: gmp_randstate_struct = uninitialized();
            __gmp_randinit_set(&mut state, &self.state);
            RandState { state: state }
        }
    }
}

pub struct Mpq {
    mpq: mpq_struct,
}

impl Drop for Mpq {
    fn drop(&mut self) { unsafe { __gmpq_clear(&mut self.mpq) } }
}

impl Mpq {
    pub fn new() -> Mpq {
        unsafe {
            let mut mpq = uninitialized();
            __gmpq_init(&mut mpq);
            Mpq { mpq: mpq }
        }
    }

    pub fn set(&mut self, other: &Mpq) {
        unsafe { __gmpq_set(&mut self.mpq, &other.mpq) }
    }

    pub fn set_z(&mut self, other: &Mpz) {
        unsafe { __gmpq_set_z(&mut self.mpq, &other.mpz) }
    }

    pub fn set_d(&mut self, other: f64) {
        unsafe { __gmpq_set_d(&mut self.mpq, other) }
    }

    pub fn set_f(&mut self, other: &Mpf) {
        unsafe { __gmpq_set_f(&mut self.mpq, &other.mpf) }
    }

    pub fn get_num(&self) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpq_get_num(&mut res.mpz, &self.mpq);
            res
        }
    }

    pub fn get_den(&self) -> Mpz {
        unsafe {
            let mut res = Mpz::new();
            __gmpq_get_den(&mut res.mpz, &self.mpq);
            res
        }
    }

    pub fn abs(&self) -> Mpq {
        unsafe {
            let mut res = Mpq::new();
            __gmpq_abs(&mut res.mpq, &self.mpq);
            res
        }
    }

    pub fn invert(&self) -> Mpq {
        unsafe {
            if self.is_zero() {
                fail!("divide by zero")
            }

            let mut res = Mpq::new();
            __gmpq_inv(&mut res.mpq, &self.mpq);
            res
        }
    }
}

impl Clone for Mpq {
    fn clone(&self) -> Mpq {
        let mut res = Mpq::new();
        res.set(self);
        res
    }
}

impl Eq for Mpq { }
impl PartialEq for Mpq {
    fn eq(&self, other: &Mpq) -> bool {
        unsafe { __gmpq_equal(&self.mpq, &other.mpq) != 0 }
    }
}

impl Ord for Mpq {
    fn cmp(&self, other: &Mpq) -> Ordering {
        let cmp = unsafe { __gmpq_cmp(&self.mpq, &other.mpq) };
        if cmp == 0 {
            Equal
        } else if cmp < 0 {
            Less
        } else {
            Greater
        }
    }
}
impl PartialOrd for Mpq {
    fn partial_cmp(&self, other: &Mpq) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Add<Mpq, Mpq> for Mpq {
    fn add(&self, other: &Mpq) -> Mpq {
        unsafe {
            let mut res = Mpq::new();
            __gmpq_add(&mut res.mpq, &self.mpq, &other.mpq);
            res
        }
    }
}

impl Sub<Mpq, Mpq> for Mpq {
    fn sub(&self, other: &Mpq) -> Mpq {
        unsafe {
            let mut res = Mpq::new();
            __gmpq_sub(&mut res.mpq, &self.mpq, &other.mpq);
            res
        }
    }
}

impl Mul<Mpq, Mpq> for Mpq {
    fn mul(&self, other: &Mpq) -> Mpq {
        unsafe {
            let mut res = Mpq::new();
            __gmpq_mul(&mut res.mpq, &self.mpq, &other.mpq);
            res
        }
    }
}

impl Div<Mpq, Mpq> for Mpq {
    fn div(&self, other: &Mpq) -> Mpq {
        unsafe {
            if other.is_zero() {
                fail!("divide by zero")
            }

            let mut res = Mpq::new();
            __gmpq_div(&mut res.mpq, &self.mpq, &other.mpq);
            res
        }
    }
}

impl Neg<Mpq> for Mpq {
    fn neg(&self) -> Mpq {
        unsafe {
            let mut res = Mpq::new();
            __gmpq_neg(&mut res.mpq, &self.mpq);
            res
        }
    }
}

impl ToPrimitive for Mpq {
    fn to_i64(&self) -> Option<i64> {
        fail!("not implemented")
    }
    fn to_u64(&self) -> Option<u64> {
        fail!("not implemented")
    }
}

impl FromPrimitive for Mpq {
    fn from_i64(other: i64) -> Option<Mpq> {
        let mut res = Mpq::new();
        res.set_z(&FromPrimitive::from_i64(other).unwrap());
        Some(res)
    }
    fn from_u64(other: u64) -> Option<Mpq> {
        let mut res = Mpq::new();
        res.set_z(&FromPrimitive::from_u64(other).unwrap());
        Some(res)
    }
}

impl One for Mpq {
    fn one() -> Mpq {
        let mut res = Mpq::new();
        unsafe { __gmpq_set_ui(&mut res.mpq, 1, 1) }
        res
    }
}

impl Zero for Mpq {
    fn zero() -> Mpq { Mpq::new() }
    fn is_zero(&self) -> bool {
        unsafe { __gmpq_cmp_ui(&self.mpq, 0, 1) == 0 }
    }
}

pub struct Mpf {
    mpf: mpf_struct,
}

impl Drop for Mpf {
    fn drop(&mut self) { unsafe { __gmpf_clear(&mut self.mpf) } }
}

impl Mpf {
    pub fn new(precision: c_ulong) -> Mpf {
        unsafe {
            let mut mpf = uninitialized();
            __gmpf_init2(&mut mpf, precision);
            Mpf { mpf: mpf }
        }
    }

    pub fn set(&mut self, other: &Mpf) {
        unsafe { __gmpf_set(&mut self.mpf, &other.mpf) }
    }

    pub fn get_prec(&self) -> c_ulong {
        unsafe { __gmpf_get_prec(&self.mpf) }
    }

    pub fn set_prec(&mut self, precision: c_ulong) {
        unsafe { __gmpf_set_prec(&mut self.mpf, precision) }
    }

    pub fn abs(&self) -> Mpf {
        unsafe {
            let mut res = Mpf::new(self.get_prec());
            __gmpf_abs(&mut res.mpf, &self.mpf);
            res
        }
    }

    pub fn ceil(&self) -> Mpf {
        unsafe {
            let mut res = Mpf::new(self.get_prec());
            __gmpf_ceil(&mut res.mpf, &self.mpf);
            res
        }
    }

    pub fn floor(&self) -> Mpf {
        unsafe {
            let mut res = Mpf::new(self.get_prec());
            __gmpf_floor(&mut res.mpf, &self.mpf);
            res
        }
    }

    pub fn trunc(&self) -> Mpf {
        unsafe {
            let mut res = Mpf::new(self.get_prec());
            __gmpf_trunc(&mut res.mpf, &self.mpf);
            res
        }
    }

    pub fn reldiff(&self, other: &Mpf) -> Mpf {
        unsafe {
            let mut res = Mpf::new(cmp::max(self.get_prec() as uint,
                                         other.get_prec() as uint) as c_ulong);
            __gmpf_reldiff(&mut res.mpf, &self.mpf, &other.mpf);
            res
        }
    }
}

impl Clone for Mpf {
    fn clone(&self) -> Mpf {
        unsafe {
            let mut mpf = uninitialized();
            __gmpf_init_set(&mut mpf, &self.mpf);
            Mpf { mpf: mpf }
        }
    }
}

impl Eq for Mpf { }
impl PartialEq for Mpf {
    fn eq(&self, other: &Mpf) -> bool {
        unsafe { __gmpf_cmp(&self.mpf, &other.mpf) == 0 }
    }
}

impl Ord for Mpf {
    fn cmp(&self, other: &Mpf) -> Ordering {
        let cmp = unsafe { __gmpf_cmp(&self.mpf, &other.mpf) };
        if cmp == 0 {
            Equal
        } else if cmp > 0 {
            Greater
        } else {
            Less
        }
    }
}
impl PartialOrd for Mpf {
    fn partial_cmp(&self, other: &Mpf) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Add<Mpf, Mpf> for Mpf {
    fn add(&self, other: &Mpf) -> Mpf {
        unsafe {
            let mut res = Mpf::new(cmp::max(self.get_prec() as uint,
                                             other.get_prec() as uint) as c_ulong);
            __gmpf_add(&mut res.mpf, &self.mpf, &other.mpf);
            res
        }
    }
}

impl Sub<Mpf, Mpf> for Mpf {
    fn sub(&self, other: &Mpf) -> Mpf {
        unsafe {
            let mut res = Mpf::new(cmp::max(self.get_prec() as uint,
                                             other.get_prec() as uint) as c_ulong);
            __gmpf_sub(&mut res.mpf, &self.mpf, &other.mpf);
            res
        }
    }
}

impl Mul<Mpf, Mpf> for Mpf {
    fn mul(&self, other: &Mpf) -> Mpf {
        unsafe {
            let mut res = Mpf::new(cmp::max(self.get_prec() as uint,
                                             other.get_prec() as uint) as c_ulong);
            __gmpf_mul(&mut res.mpf, &self.mpf, &other.mpf);
            res
        }
    }
}

impl Div<Mpf, Mpf> for Mpf {
    fn div(&self, other: &Mpf) -> Mpf {
        unsafe {
            if __gmpf_cmp_ui(&self.mpf, 0) == 0 {
                fail!("divide by zero")
            }

            let mut res = Mpf::new(cmp::max(self.get_prec() as uint,
                                             other.get_prec() as uint) as c_ulong);
            __gmpf_div(&mut res.mpf, &self.mpf, &other.mpf);
            res
        }
    }
}

impl Neg<Mpf> for Mpf {
    fn neg(&self) -> Mpf {
        unsafe {
            let mut res = Mpf::new(self.get_prec());
            __gmpf_neg(&mut res.mpf, &self.mpf);
            res
        }
    }
}
