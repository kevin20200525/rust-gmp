extern mod std;

#[abi = "rust-intrinsic"]
extern mod rusti {
  fn init<T>() -> T;
}

use clone::Clone;
use from_str::FromStr;
use libc::{c_char,c_int,c_ulong,c_void,size_t};
use num::{Num, One, Zero};
use ptr::{addr_of,mut_addr_of,to_mut_unsafe_ptr};
use str::as_c_str;

struct mpz_struct {
  _mp_alloc: c_int,
  _mp_size: c_int,
  _mp_d: *c_void
}

struct mpq_struct {
  _mp_num: mpz_struct,
  _mp_den: mpz_struct
}

struct gmp_randstate_struct {
  _mp_seed: mpz_struct,
  _mp_alg: c_int,
  _mp_algdata: *c_void
}

type mp_bitcnt_t = c_ulong;
type mpz_srcptr = *const mpz_struct;
type mpz_ptr = *mut mpz_struct;
type mpq_srcptr = *const mpq_struct;
type mpq_ptr = *mut mpq_struct;
type gmp_randstate_t = *mut gmp_randstate_struct;

extern mod gmp {
  fn __gmpz_init(x: mpz_ptr);
  fn __gmpz_init_set(rop: mpz_ptr, op: mpz_srcptr);
  fn __gmpz_init_set_ui(rop: mpz_ptr, op: c_ulong);
  fn __gmpz_init_set_str(rop: mpz_ptr, str: *c_char, base: c_int) -> c_int;
  fn __gmpz_clear(x: mpz_ptr);
  fn __gmpz_set(rop: mpz_ptr, op: mpz_srcptr);
  fn __gmpz_set_str(rop: mpz_ptr, str: *c_char, base: c_int) -> c_int;
  fn __gmpz_get_str(str: *c_char, base: c_int, op: mpz_srcptr) -> *c_char;
  pure fn __gmpz_sizeinbase(op: mpz_srcptr, base: c_int) -> size_t;
  pure fn __gmpz_cmp(op1: mpz_srcptr, op2: mpz_srcptr) -> c_int;
  pure fn __gmpz_cmp_ui(op1: mpz_srcptr, op2: c_ulong) -> c_int;
  fn __gmpz_add(rop: mpz_ptr, op1: mpz_srcptr, op2: mpz_srcptr);
  fn __gmpz_sub(rop: mpz_ptr, op1: mpz_srcptr, op2: mpz_srcptr);
  fn __gmpz_mul(rop: mpz_ptr, op1: mpz_srcptr, op2: mpz_srcptr);
  fn __gmpz_mul_2exp(rop: mpz_ptr, op1: mpz_srcptr, op2: mp_bitcnt_t);
  fn __gmpz_neg(rop: mpz_ptr, op: mpz_srcptr);
  fn __gmpz_abs(rop: mpz_ptr, op: mpz_srcptr);
  fn __gmpz_tdiv_q(q: mpz_ptr, n: mpz_srcptr, d: mpz_srcptr);
  fn __gmpz_fdiv_q_2exp(q: mpz_ptr, n: mpz_srcptr, b: mp_bitcnt_t);
  fn __gmpz_mod(r: mpz_ptr, n: mpz_srcptr, d: mpz_srcptr);
  fn __gmpz_and(rop: mpz_ptr, op1: mpz_srcptr, op2: mpz_srcptr);
  fn __gmpz_ior(rop: mpz_ptr, op1: mpz_srcptr, op2: mpz_srcptr);
  fn __gmpz_xor(rop: mpz_ptr, op1: mpz_srcptr, op2: mpz_srcptr);
  fn __gmpz_com(rop: mpz_ptr, op: mpz_srcptr);
  pure fn __gmpz_popcount(op: mpz_srcptr) -> mp_bitcnt_t;
  pure fn __gmpz_hamdist(op1: mpz_srcptr, op2: mpz_srcptr) -> mp_bitcnt_t;
  fn __gmpz_gcd(rop: mpz_ptr, op1: mpz_srcptr, op2: mpz_srcptr);
  fn __gmpz_lcm(rop: mpz_ptr, op1: mpz_srcptr, op2: mpz_srcptr);
  fn __gmpz_invert(rop: mpz_ptr, op1: mpz_srcptr, op2: mpz_srcptr) -> c_int;
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
  pure fn __gmpq_cmp(op1: mpq_srcptr, op2: mpq_srcptr) -> c_int;
  pure fn __gmpq_equal(op1: mpq_srcptr, op2: mpq_srcptr) -> c_int;
  fn __gmpq_add(sum: mpq_ptr, addend1: mpq_srcptr, addend2: mpq_srcptr);
  fn __gmpq_sub(difference: mpq_ptr, minuend: mpq_srcptr, subtrahend: mpq_srcptr);
  fn __gmpq_mul(product: mpq_ptr, multiplier: mpq_srcptr, multiplicand: mpq_srcptr);
  fn __gmpq_div(product: mpq_ptr, multiplier: mpq_srcptr, multiplicand: mpq_srcptr);
  fn __gmpq_neg(negated_operand: mpq_ptr, operand: mpq_srcptr);
  fn __gmpq_abs(rop: mpq_ptr, op: mpq_srcptr);
}

use gmp::*;

pub struct Mpz {
  priv mpz: mpz_struct,

  drop {
    __gmpz_clear(mut_addr_of(&self.mpz));
  }
}

impl Mpz {
  fn set(&mut self, other: &Mpz) {
    __gmpz_set(mut_addr_of(&self.mpz), addr_of(&other.mpz));
  }

  fn set_from_str_radix(&mut self, s: &str, base: uint) -> bool {
    assert base == 0 || base >= 2 || base <= 62;
    let mpz = to_mut_unsafe_ptr(&mut self.mpz);
    as_c_str(s, { |s| __gmpz_set_str(mpz, s, base as c_int) }) == 0
  }

  pure fn to_str_radix(&self, base: int) -> ~str unsafe {
    let len = __gmpz_sizeinbase(addr_of(&self.mpz), base as c_int) as uint + 2;
    let dst = vec::to_mut(vec::from_elem(len, '0'));
    let pdst = vec::raw::to_ptr(dst);

    str::raw::from_c_str(__gmpz_get_str(pdst as *c_char, base as c_int,
                         addr_of(&self.mpz)))
  }

  pure fn bit_length(&self) -> uint {
    return __gmpz_sizeinbase(addr_of(&self.mpz), 2) as uint
  }

  pure fn compl(&self) -> Mpz unsafe {
    let res = Mpz::new();
    __gmpz_com(mut_addr_of(&res.mpz), addr_of(&self.mpz));
    res
  }

  pure fn abs(&self) -> Mpz unsafe {
    let res = Mpz::new();
    __gmpz_abs(mut_addr_of(&res.mpz), addr_of(&self.mpz));
    res
  }

  pure fn gcd(&self, other: &Mpz) -> Mpz unsafe {
    let res = Mpz::new();
    __gmpz_gcd(mut_addr_of(&res.mpz), addr_of(&self.mpz), addr_of(&other.mpz));
    res
  }

  pure fn lcm(&self, other: &Mpz) -> Mpz unsafe {
    let res = Mpz::new();
    __gmpz_lcm(mut_addr_of(&res.mpz), addr_of(&self.mpz), addr_of(&other.mpz));
    res
  }

  pure fn invert(&self, modulo: &Mpz) -> Option<Mpz> unsafe {
    let res = Mpz::new();
    if __gmpz_invert(mut_addr_of(&res.mpz), addr_of(&self.mpz),
                     addr_of(&modulo.mpz)) == 0 {
      None
    } else {
      Some(res)
    }
  }

  pure fn popcount(&self) -> uint {
    __gmpz_popcount(addr_of(&self.mpz)) as uint
  }

  pure fn hamdist(&self, other: &Mpz) -> uint {
    __gmpz_hamdist(addr_of(&self.mpz), addr_of(&other.mpz)) as uint
  }

  static pure fn new() -> Mpz unsafe {
    let mpz = rusti::init(); // TODO: switch to rusti::uninit when implemented
    __gmpz_init(mut_addr_of(&mpz));
    Mpz { mpz: mpz }
  }

  static pure fn from_str_radix(s: &str, base: uint) -> Option<Mpz> unsafe {
    assert base == 0 || base >= 2 || base <= 62;
    let mpz = rusti::init(); // TODO: switch to rusti::uninit when implemented
    let mpz_ptr = mut_addr_of(&mpz);
    let r = as_c_str(s, { |s| __gmpz_init_set_str(mpz_ptr, s, base as c_int) });
    if r == 0 {
      Some(Mpz { mpz: mpz })
    } else {
      __gmpz_clear(mpz_ptr);
      None
    }
  }
}

impl Mpz: Clone {
  pure fn clone(&self) -> Mpz unsafe {
    let mpz = rusti::init(); // TODO: switch to rusti::uninit when implemented
    __gmpz_init_set(mut_addr_of(&mpz), addr_of(&self.mpz));
    Mpz { mpz: mpz }
  }
}

impl Mpz: cmp::Eq {
  pure fn eq(&self, other: &Mpz) -> bool {
    __gmpz_cmp(addr_of(&self.mpz), addr_of(&other.mpz)) == 0
  }
  pure fn ne(&self, other: &Mpz) -> bool {
    __gmpz_cmp(addr_of(&self.mpz), addr_of(&other.mpz)) != 0
  }
}

impl Mpz: cmp::Ord {
  pure fn lt(&self, other: &Mpz) -> bool {
    __gmpz_cmp(addr_of(&self.mpz), addr_of(&other.mpz)) < 0
  }
  pure fn le(&self, other: &Mpz) -> bool {
    __gmpz_cmp(addr_of(&self.mpz), addr_of(&other.mpz)) <= 0
  }
  pure fn gt(&self, other: &Mpz) -> bool {
    __gmpz_cmp(addr_of(&self.mpz), addr_of(&other.mpz)) > 0
  }
  pure fn ge(&self, other: &Mpz) -> bool {
    __gmpz_cmp(addr_of(&self.mpz), addr_of(&other.mpz)) >= 0
  }
}

impl Mpz: Num {
  pure fn add(&self, other: &Mpz) -> Mpz unsafe {
    let res = Mpz::new();
    __gmpz_add(mut_addr_of(&res.mpz), addr_of(&self.mpz), addr_of(&other.mpz));
    res
  }
  pure fn sub(&self, other: &Mpz) -> Mpz unsafe {
    let res = Mpz::new();
    __gmpz_sub(mut_addr_of(&res.mpz), addr_of(&self.mpz), addr_of(&other.mpz));
    res
  }
  pure fn mul(&self, other: &Mpz) -> Mpz unsafe {
    let res = Mpz::new();
    __gmpz_mul(mut_addr_of(&res.mpz), addr_of(&self.mpz), addr_of(&other.mpz));
    res
  }
  pure fn div(&self, other: &Mpz) -> Mpz unsafe {
    if __gmpz_cmp_ui(addr_of(&self.mpz), 0) == 0 {
      fail ~"divide by zero";
    }

    let res = Mpz::new();
    __gmpz_tdiv_q(mut_addr_of(&res.mpz), addr_of(&self.mpz), addr_of(&other.mpz));
    res
  }
  pure fn modulo(&self, other: &Mpz) -> Mpz unsafe {
    if __gmpz_cmp_ui(addr_of(&self.mpz), 0) == 0 {
      fail ~"divide by zero";
    }

    let res = Mpz::new();
    __gmpz_mod(mut_addr_of(&res.mpz), addr_of(&self.mpz), addr_of(&other.mpz));
    res
  }
  pure fn neg(&self) -> Mpz unsafe {
    let res = Mpz::new();
    __gmpz_neg(mut_addr_of(&res.mpz), addr_of(&self.mpz));
    res
  }
  pure fn to_int(&self) -> int {
    fail ~"not implemented";
  }
  static pure fn from_int(&self, other: int) -> Mpz {
    // the gmp functions dealing with longs aren't usable here - long is only
    // guaranteed to be at least 32-bit
    FromStr::from_str(other.to_str()).unwrap()
  }
}

impl Mpz: One {
  static pure fn one() -> Mpz unsafe {
    let mpz = rusti::init(); // TODO: switch to rusti::uninit when implemented
    __gmpz_init_set_ui(mut_addr_of(&mpz), 1);
    Mpz { mpz: mpz }
  }
}

impl Mpz: Zero {
  static pure fn zero() -> Mpz {
    Mpz::new()
  }
}

impl Mpz: BitAnd<Mpz, Mpz> {
  pure fn bitand(&self, other: &Mpz) -> Mpz unsafe {
    let res = Mpz::new();
    __gmpz_and(mut_addr_of(&res.mpz), addr_of(&self.mpz), addr_of(&other.mpz));
    res
  }
}

impl Mpz: BitOr<Mpz, Mpz> {
  pure fn bitor(&self, other: &Mpz) -> Mpz unsafe {
    let res = Mpz::new();
    __gmpz_ior(mut_addr_of(&res.mpz), addr_of(&self.mpz), addr_of(&other.mpz));
    res
  }
}

impl Mpz: BitXor<Mpz, Mpz> {
  pure fn bitxor(&self, other: &Mpz) -> Mpz unsafe {
    let res = Mpz::new();
    __gmpz_xor(mut_addr_of(&res.mpz), addr_of(&self.mpz), addr_of(&other.mpz));
    res
  }
}

impl Mpz: Shl<c_ulong, Mpz> {
  pure fn shl(&self, other: &c_ulong) -> Mpz unsafe {
    let res = Mpz::new();
    __gmpz_mul_2exp(mut_addr_of(&res.mpz), addr_of(&self.mpz), *other);
    res
  }
}

impl Mpz: Shr<c_ulong, Mpz> {
  pure fn shr(&self, other: &c_ulong) -> Mpz unsafe {
    let res = Mpz::new();
    __gmpz_fdiv_q_2exp(mut_addr_of(&res.mpz), addr_of(&self.mpz), *other);
    res
  }
}

impl Mpz : FromStr {
  static pure fn from_str(s: &str) -> Option<Mpz> {
    Mpz::from_str_radix(s, 10)
  }
}

impl Mpz : to_str::ToStr {
  pure fn to_str() -> ~str {
    self.to_str_radix(10)
  }
}

pub struct RandState {
  priv state: gmp_randstate_struct,

  drop {
    __gmp_randclear(mut_addr_of(&self.state));
  }
}

impl RandState {
  static pure fn new() -> RandState unsafe {
    let state: gmp_randstate_struct = rusti::init();
    __gmp_randinit_default(mut_addr_of(&state));
    RandState { state: state }
  }

  static pure fn new_mt() -> RandState unsafe {
    let state: gmp_randstate_struct = rusti::init();
    __gmp_randinit_mt(mut_addr_of(&state));
    RandState { state: state }
  }

  static pure fn new_lc_2exp(a: Mpz, c: c_ulong, m2exp: c_ulong) -> RandState unsafe {
    let state: gmp_randstate_struct = rusti::init();
    __gmp_randinit_lc_2exp(mut_addr_of(&state), addr_of(&a.mpz), c, m2exp);
    RandState { state: state }
  }

  static pure fn new_lc_2exp_size(size: c_ulong) -> RandState unsafe {
    let state: gmp_randstate_struct = rusti::init();
    __gmp_randinit_lc_2exp_size(mut_addr_of(&state), size);
    RandState { state: state }
  }

  fn seed(&mut self, seed: Mpz) {
    __gmp_randseed(mut_addr_of(&self.state), addr_of(&seed.mpz));
  }

  fn seed_ui(&mut self, seed: c_ulong) {
    __gmp_randseed_ui(mut_addr_of(&self.state), seed);
  }

  /// Generate a uniform random integer in the range 0 to n-1, inclusive.
  fn urandom(n: &Mpz) -> Mpz {
    let res = Mpz::new();
    __gmpz_urandomm(mut_addr_of(&res.mpz), mut_addr_of(&self.state), addr_of(&n.mpz));
    res
  }

  /// Generate a uniformly distributed random integer in the range 0 to 2^n−1, inclusive.
  fn urandom_2exp(n: c_ulong) -> Mpz {
    let res = Mpz::new();
    __gmpz_urandomb(mut_addr_of(&res.mpz), mut_addr_of(&self.state), n);
    res
  }
}

impl RandState: Clone {
  pure fn clone(&self) -> RandState unsafe {
    let state: gmp_randstate_struct = rusti::init();
    __gmp_randinit_set(mut_addr_of(&state), addr_of(&self.state));
    RandState { state: state }
  }
}

pub struct Mpq {
  priv mpq: mpq_struct,

  drop {
    __gmpq_clear(mut_addr_of(&self.mpq));
  }
}

impl Mpq {
  static pure fn new() -> Mpq unsafe {
    let mpq = rusti::init(); // TODO: switch to rusti::uninit when implemented
    __gmpq_init(mut_addr_of(&mpq));
    Mpq { mpq: mpq }
  }

  fn set(&mut self, other: &Mpq) {
    __gmpq_set(mut_addr_of(&self.mpq), addr_of(&other.mpq));
  }

  fn set_z(&mut self, other: &Mpz) {
    __gmpq_set_z(mut_addr_of(&self.mpq), addr_of(&other.mpz));
  }

  pure fn abs(&self) -> Mpq unsafe {
    let res = Mpq::new();
    __gmpq_abs(mut_addr_of(&res.mpq), addr_of(&self.mpq));
    res
  }
}

impl Mpq: Clone {
  pure fn clone(&self) -> Mpq unsafe {
    let mut res = Mpq::new();
    res.set(self);
    res
  }
}

impl Mpq: cmp::Eq {
  pure fn eq(&self, other: &Mpq) -> bool {
    __gmpq_equal(addr_of(&self.mpq), addr_of(&other.mpq)) != 0
  }
  pure fn ne(&self, other: &Mpq) -> bool {
    __gmpq_equal(addr_of(&self.mpq), addr_of(&other.mpq)) == 0
  }
}

impl Mpq: cmp::Ord {
  pure fn lt(&self, other: &Mpq) -> bool {
    __gmpq_cmp(addr_of(&self.mpq), addr_of(&other.mpq)) < 0
  }
  pure fn le(&self, other: &Mpq) -> bool {
    __gmpq_cmp(addr_of(&self.mpq), addr_of(&other.mpq)) <= 0
  }
  pure fn gt(&self, other: &Mpq) -> bool {
    __gmpq_cmp(addr_of(&self.mpq), addr_of(&other.mpq)) > 0
  }
  pure fn ge(&self, other: &Mpq) -> bool {
    __gmpq_cmp(addr_of(&self.mpq), addr_of(&other.mpq)) >= 0
  }
}

impl Mpq: Num {
  pure fn add(&self, other: &Mpq) -> Mpq unsafe {
    let res = Mpq::new();
    __gmpq_add(mut_addr_of(&res.mpq), addr_of(&self.mpq), addr_of(&other.mpq));
    res
  }
  pure fn sub(&self, other: &Mpq) -> Mpq unsafe {
    let res = Mpq::new();
    __gmpq_sub(mut_addr_of(&res.mpq), addr_of(&self.mpq), addr_of(&other.mpq));
    res
  }
  pure fn mul(&self, other: &Mpq) -> Mpq unsafe {
    let res = Mpq::new();
    __gmpq_mul(mut_addr_of(&res.mpq), addr_of(&self.mpq), addr_of(&other.mpq));
    res
  }
  // TODO: handle division by zero
  pure fn div(&self, other: &Mpq) -> Mpq unsafe {
    let res = Mpq::new();
    __gmpq_div(mut_addr_of(&res.mpq), addr_of(&self.mpq), addr_of(&other.mpq));
    res
  }
  pure fn modulo(&self, _other: &Mpq) -> Mpq unsafe {
    fail ~"not implemented";
  }
  pure fn neg(&self) -> Mpq unsafe {
    let res = Mpq::new();
    __gmpq_neg(mut_addr_of(&res.mpq), addr_of(&self.mpq));
    res
  }
  pure fn to_int(&self) -> int {
    fail ~"not implemented";
  }
  static pure fn from_int(&self, _other: int) -> Mpq {
    fail ~"not implemented";
  }
}

#[cfg(test)]
mod test_mpz {
  use Num::from_int;
  use FromStr::from_str;

  #[test]
  fn test_set() {
    let mut x: Mpz = from_int(1000);
    let y: Mpz = from_int(5000);
    assert x != y;
    x.set(&y);
    assert x == y;
  }

  #[test]
  fn test_set_from_str_radix() {
    let mut x: Mpz = from_int(1000);
    let y: Mpz = from_int(5000);
    assert x != y;
    assert x.set_from_str_radix("5000", 10);
    assert x == y;
    assert !x.set_from_str_radix("aaaa", 2);
  }

  #[test]
  fn test_eq() {
    let x: Mpz = from_int(4242142195);
    let y: Mpz = from_int(4242142195);
    let z: Mpz = from_int(4242142196);

    assert x == y;
    assert x != z;
    assert y != z;
  }

  #[test]
  fn test_ord() {
    let x: Mpz = from_str("40000000000000000000000").unwrap();
    let y: Mpz = from_str("45000000000000000000000").unwrap();
    let z: Mpz = from_str("50000000000000000000000").unwrap();

    assert x < y && x < z && y < z;
    assert x <= x && x <= y && x <= z && y <= z;
    assert z > y && z > x && y > x;
    assert z >= z && z >= y && z >= x && y >= x;
  }

  #[test]
  #[should_fail]
  fn test_div_zero() {
    let x = Mpz::new();
    x / x;
  }

  #[test]
  #[should_fail]
  fn test_modulo_zero() {
    let x = Mpz::new();
    x % x;
  }

  #[test]
  fn test_div_round() {
    let x: Mpz = from_int(2);
    let y: Mpz = from_int(3);
    assert (x / y).to_str() == (2 / 3).to_str();
    assert (x / -y).to_str() == (2 / -3).to_str();
  }

  #[test]
  fn test_to_str_radix() {
    let x: Mpz = from_int(255);
    assert x.to_str_radix(16) == ~"ff";
  }

  #[test]
  fn test_to_str() {
    let x: Mpz = from_str("1234567890").unwrap();
    assert x.to_str() == ~"1234567890";
  }

  #[test]
  fn test_invalid_str() {
    assert from_str::<Mpz>("foobar").is_none();
  }

  #[test]
  fn test_clone() {
    let a = from_int::<Mpz>(100);
    let b = a.clone();
    assert b == a;
    assert a + b == from_int::<Mpz>(200);
  }

  #[test]
  fn test_from_int() {
    let x: Mpz = from_int(150);
    assert x.to_str() == ~"150";
    assert x == from_str("150").unwrap();
  }

  #[test]
  fn test_abs() {
    let x: Mpz = from_int(1000);
    let y: Mpz = from_int(-1000);
    assert -x == y;
    assert x == -y;
    assert x == y.abs();
    assert x.abs() == y.abs();
  }

  #[test]
  fn test_bitand() {
    let a = 0b1001_0111;
    let b = 0b1100_0100;
    assert from_int::<Mpz>(a) & from_int::<Mpz>(b) == from_int::<Mpz>(a & b);
  }

  #[test]
  fn test_bitor() {
    let a = 0b1001_0111;
    let b = 0b1100_0100;
    assert from_int::<Mpz>(a) | from_int::<Mpz>(b) == from_int::<Mpz>(a | b);
  }

  #[test]
  fn test_bitxor() {
    let a = 0b1001_0111;
    let b = 0b1100_0100;
    assert from_int::<Mpz>(a) ^ from_int::<Mpz>(b) == from_int::<Mpz>(a ^ b);
  }

  #[test]
  fn test_shifts() {
    let i = 227;
    let j: Mpz = from_int(i);
    assert (i << 4).to_str() == (j << 4).to_str();
    assert (-i << 4).to_str() == (-j << 4).to_str();
    assert (i >> 4).to_str() == (j >> 4).to_str();
    assert (-i >> 4).to_str() == (-j >> 4).to_str();
  }

  #[test]
  fn test_compl() {
    assert from_int::<Mpz>(13).compl().to_str() == (!13).to_str();
    assert from_int::<Mpz>(-442).compl().to_str() == (!-442).to_str();
  }

  #[test]
  fn test_popcount() {
    Mpz::from_str_radix("1010010011", 2).unwrap().popcount() == 5;
  }

  #[test]
  fn test_hamdist() {
    assert from_int::<Mpz>(0b1011_0001).hamdist(&from_int(0b0010_1011)) == 4;
  }

  #[test]
  fn test_bit_length() {
    assert from_int::<Mpz>(0b1011_0000_0001_0000).bit_length() == 16;
    assert from_int::<Mpz>(0b101).bit_length() == 3;
  }

  #[test]
  fn test_gcd() {
    assert from_int::<Mpz>(0).gcd(&from_int(0)) == from_int(0);
    assert from_int::<Mpz>(3).gcd(&from_int(6)) == from_int(3);
    assert from_int::<Mpz>(18).gcd(&from_int(24)) == from_int(6);
  }

  #[test]
  fn test_lcm() {
    assert from_int::<Mpz>(0).lcm(&from_int(5)) == from_int(0);
    assert from_int::<Mpz>(5).lcm(&from_int(0)) == from_int(0);
    assert from_int::<Mpz>(3).lcm(&from_int(6)) == from_int(6);
    assert from_int::<Mpz>(18).lcm(&from_int(24)) == from_int(72);
  }

  #[test]
  fn test_invert() {
    assert from_int::<Mpz>(3).invert(&from_int(11)) == Some(from_int(4));
    assert from_int::<Mpz>(4).invert(&from_int(11)) == Some(from_int(3));
    assert from_int::<Mpz>(2).invert(&from_int(5)) == Some(from_int(3));
    assert from_int::<Mpz>(3).invert(&from_int(5)) == Some(from_int(2));
    assert from_int::<Mpz>(2).invert(&from_int(4)).is_none();
  }

  #[test]
  fn test_one() {
    assert One::one::<Mpz>() == Num::from_int(1);
  }
}

#[cfg(test)]
mod test_rand {
  #[test]
  fn test_randstate() {
    let mut state = RandState::new();
    state.seed_ui(42);
    for uint::range(1, 1000) |_| {
      for int::range(1, 10) |x| {
        let upper = Num::from_int(x);
        assert state.urandom(&upper) < upper;
      }
    }
  }
}

#[cfg(test)]
mod test_mpq {
  #[test]
  fn test_mpq() {
    assert Mpq::new() == Mpq::new();
  }
}
