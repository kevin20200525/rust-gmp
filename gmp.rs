extern mod std;

use libc::{c_char,c_int,c_ulong,c_void,size_t};
use num::{from_int,Num};
use ptr::{addr_of,mut_addr_of,null,to_mut_unsafe_ptr};
use str::as_c_str;

struct mpz_struct {
  _mp_alloc: c_int,
  _mp_size: c_int,
  _mp_d: *c_void
}

type mp_bitcnt_t = c_ulong;
type mpz_srcptr = *const mpz_struct;
type mpz_ptr = *mut mpz_struct;

extern mod gmp {
  fn __gmpz_init(x: mpz_ptr);
  fn __gmpz_init_set(rop: mpz_ptr, op: mpz_srcptr);
  fn __gmpz_init_set_str(rop: mpz_ptr, str: *c_char, base: c_int) -> c_int;
  fn __gmpz_clear(x: mpz_ptr);
  fn __gmpz_set(rop: mpz_ptr, op: mpz_srcptr);
  fn __gmpz_set_str(rop: mpz_ptr, str: *c_char, base: c_int) -> c_int;
  fn __gmpz_get_str(str: *c_char, base: c_int, op: mpz_srcptr) -> *c_char;
  pure fn __gmpz_sizeinbase(op: mpz_srcptr, base: c_int) -> size_t;
  pure fn __gmpz_cmp(op: mpz_srcptr, op2: mpz_srcptr) -> c_int;
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
    assert(base == 0 || base >= 2 || base <= 62);
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

  pure fn bit_length() -> uint {
    return __gmpz_sizeinbase(addr_of(&self.mpz), 2) as uint
  }

  // TODO: implement the clone::Clone trait when 0.5 is released
  pure fn clone() -> Mpz unsafe {
    let mpz = mpz_struct { _mp_alloc: 0, _mp_size: 0, _mp_d: null() };
    __gmpz_init_set(mut_addr_of(&mpz), addr_of(&self.mpz));
    Mpz { mpz: mpz }
  }

  pure fn compl() -> Mpz unsafe {
    let res = init();
    __gmpz_com(mut_addr_of(&res.mpz), addr_of(&self.mpz));
    res
  }

  pure fn abs() -> Mpz unsafe {
    let res = init();
    __gmpz_abs(mut_addr_of(&res.mpz), addr_of(&self.mpz));
    res
  }

  pure fn popcount() -> uint {
    __gmpz_popcount(addr_of(&self.mpz)) as uint
  }

  pure fn hamdist(other: &Mpz) -> uint {
    __gmpz_hamdist(addr_of(&self.mpz), addr_of(&other.mpz)) as uint
  }
}

impl Mpz: cmp::Eq {
  pure fn eq(other: &Mpz) -> bool {
    __gmpz_cmp(addr_of(&self.mpz), addr_of(&other.mpz)) == 0
  }
  pure fn ne(other: &Mpz) -> bool {
    __gmpz_cmp(addr_of(&self.mpz), addr_of(&other.mpz)) != 0
  }
}

impl Mpz: cmp::Ord {
  pure fn lt(other: &Mpz) -> bool {
    __gmpz_cmp(addr_of(&self.mpz), addr_of(&other.mpz)) < 0
  }
  pure fn le(other: &Mpz) -> bool {
    __gmpz_cmp(addr_of(&self.mpz), addr_of(&other.mpz)) <= 0
  }
  pure fn gt(other: &Mpz) -> bool {
    __gmpz_cmp(addr_of(&self.mpz), addr_of(&other.mpz)) > 0
  }
  pure fn ge(other: &Mpz) -> bool {
    __gmpz_cmp(addr_of(&self.mpz), addr_of(&other.mpz)) >= 0
  }
}

impl Mpz: Num {
  pure fn add(other: &Mpz) -> Mpz unsafe {
    let res = init();
    __gmpz_add(mut_addr_of(&res.mpz), addr_of(&self.mpz), addr_of(&other.mpz));
    res
  }
  pure fn sub(other: &Mpz) -> Mpz unsafe {
    let res = init();
    __gmpz_sub(mut_addr_of(&res.mpz), addr_of(&self.mpz), addr_of(&other.mpz));
    res
  }
  pure fn mul(other: &Mpz) -> Mpz unsafe {
    let res = init();
    __gmpz_mul(mut_addr_of(&res.mpz), addr_of(&self.mpz), addr_of(&other.mpz));
    res
  }
  pure fn div(other: &Mpz) -> Mpz unsafe {
    if __gmpz_cmp_ui(addr_of(&self.mpz), 0) == 0 {
      fail ~"divide by zero";
    }

    let res = init();
    __gmpz_tdiv_q(mut_addr_of(&res.mpz), addr_of(&self.mpz), addr_of(&other.mpz));
    res
  }
  pure fn modulo(other: &Mpz) -> Mpz unsafe {
    if __gmpz_cmp_ui(addr_of(&self.mpz), 0) == 0 {
      fail ~"divide by zero";
    }

    let res = init();
    __gmpz_mod(mut_addr_of(&res.mpz), addr_of(&self.mpz), addr_of(&other.mpz));
    res
  }
  pure fn neg() -> Mpz unsafe {
    let res = init();
    __gmpz_neg(mut_addr_of(&res.mpz), addr_of(&self.mpz));
    res
  }
  pure fn to_int() -> int {
    fail ~"not implemented";
  }
  static pure fn from_int(other: int) -> Mpz {
    // the gmp functions dealing with longs aren't usable here - long is only
    // guaranteed to be at least 32-bit
    option::unwrap(from_str(other.to_str()))
  }
}

impl Mpz: BitAnd<Mpz, Mpz> {
  pure fn bitand(other: &Mpz) -> Mpz unsafe {
    let res = init();
    __gmpz_and(mut_addr_of(&res.mpz), addr_of(&self.mpz), addr_of(&other.mpz));
    res
  }
}

impl Mpz: BitOr<Mpz, Mpz> {
  pure fn bitor(other: &Mpz) -> Mpz unsafe {
    let res = init();
    __gmpz_ior(mut_addr_of(&res.mpz), addr_of(&self.mpz), addr_of(&other.mpz));
    res
  }
}

impl Mpz: BitXor<Mpz, Mpz> {
  pure fn bitxor(other: &Mpz) -> Mpz unsafe {
    let res = init();
    __gmpz_xor(mut_addr_of(&res.mpz), addr_of(&self.mpz), addr_of(&other.mpz));
    res
  }
}

impl Mpz: Shl<c_ulong, Mpz> {
  pure fn shl(other: &c_ulong) -> Mpz unsafe {
    let res = init();
    __gmpz_mul_2exp(mut_addr_of(&res.mpz), addr_of(&self.mpz), *other);
    res
  }
}

impl Mpz: Shr<c_ulong, Mpz> {
  pure fn shr(other: &c_ulong) -> Mpz unsafe {
    let res = init();
    __gmpz_fdiv_q_2exp(mut_addr_of(&res.mpz), addr_of(&self.mpz), *other);
    res
  }
}

pub pure fn from_str_radix(s: &str, base: uint) -> Option<Mpz> unsafe {
  assert(base == 0 || base >= 2 || base <= 62);
  let mpz = mpz_struct { _mp_alloc: 0, _mp_size: 0, _mp_d: null() };
  let mpz_ptr = mut_addr_of(&mpz);
  let r = as_c_str(s, { |s| __gmpz_init_set_str(mpz_ptr, s, base as c_int) });
  if r == 0 {
    Some(Mpz { mpz: mpz })
  } else {
    __gmpz_clear(mpz_ptr);
    None
  }
}

pub pure fn from_str(s: &str) -> Option<Mpz> {
  from_str_radix(s, 10)
}

impl Mpz : from_str::FromStr {
  static fn from_str(s: &str) -> Option<Mpz> {
    from_str(s)
  }
}

impl Mpz : to_str::ToStr {
  pure fn to_str() -> ~str {
    self.to_str_radix(10)
  }
}

pub pure fn init() -> Mpz unsafe {
  let mpz = mpz_struct { _mp_alloc: 0, _mp_size: 0, _mp_d: null() };
  __gmpz_init(mut_addr_of(&mpz));
  Mpz { mpz: mpz }
}

#[cfg(test)]
mod test_mpz {
  #[test]
  fn test_set() {
    let mut x: Mpz = from_int(1000);
    let y: Mpz = from_int(5000);
    assert(x != y);
    x.set(&y);
    assert(x == y);
  }

  #[test]
  fn test_set_from_str_radix() {
    let mut x: Mpz = from_int(1000);
    let y: Mpz = from_int(5000);
    assert(x != y);
    assert(x.set_from_str_radix("5000", 10));
    assert(x == y);
    assert(!x.set_from_str_radix("aaaa", 2));
  }

  #[test]
  fn test_eq() {
    let x = option::unwrap(from_str("4242142195"));
    let y = option::unwrap(from_str("4242142195"));
    let z = option::unwrap(from_str("4242142196"));

    assert(x == y);
    assert(x != z);
    assert(y != z);
  }

  #[test]
  fn test_ord() {
    let x = option::unwrap(from_str("40000000000000000000000"));
    let y = option::unwrap(from_str("45000000000000000000000"));
    let z = option::unwrap(from_str("50000000000000000000000"));

    assert(x < y && x < z && y < z);
    assert(x <= x && x <= y && x <= z && y <= z);
    assert(z > y && z > x && y > x);
    assert(z >= z && z >= y && z >= x && y >= x);
  }

  #[test]
  #[should_fail]
  fn test_div_zero() {
    let x = init();
    x / x;
  }

  #[test]
  #[should_fail]
  fn test_modulo_zero() {
    let x = init();
    x % x;
  }

  #[test]
  fn test_div_round() {
    let x: Mpz = from_int(2);
    let y: Mpz = from_int(3);
    assert((x / y).to_str() == (2 / 3).to_str());
    assert((x / -y).to_str() == (2 / -3).to_str());
  }

  #[test]
  fn test_to_str_radix() {
    let x = option::unwrap(from_str("255"));
    assert(x.to_str_radix(16) == ~"ff");
  }

  #[test]
  fn test_to_str() {
    let x = option::unwrap(from_str("1234567890"));
    assert(x.to_str() == ~"1234567890");
  }

  #[test]
  fn test_invalid_str() {
    assert(from_str("foobar").is_none());
  }

  #[test]
  fn test_clone() {
    let a = option::unwrap(from_str("100"));
    let b = a.clone();
    assert(b == a);
    assert(a + b == option::unwrap(from_str("200")));
  }

  #[test]
  fn test_from_int() {
    let x: Mpz = from_int(150);
    assert(x.to_str() == ~"150");
    assert(x == option::unwrap(from_str("150")));
  }

  #[test]
  fn test_abs() {
    let x: Mpz = from_int(1000);
    let y: Mpz = from_int(-1000);
    assert(-x == y);
    assert(x == -y);
    assert(x == y.abs());
    assert(x.abs() == y.abs());
  }

  #[test]
  fn test_bitand() {
    let a = 0b1001_0111;
    let b = 0b1100_0100;
    assert(from_int::<Mpz>(a) & from_int::<Mpz>(b) == from_int::<Mpz>(a & b));
  }

  #[test]
  fn test_bitor() {
    let a = 0b1001_0111;
    let b = 0b1100_0100;
    assert(from_int::<Mpz>(a) | from_int::<Mpz>(b) == from_int::<Mpz>(a | b));
  }

  #[test]
  fn test_bitxor() {
    let a = 0b1001_0111;
    let b = 0b1100_0100;
    assert(from_int::<Mpz>(a) ^ from_int::<Mpz>(b) == from_int::<Mpz>(a ^ b));
  }

  #[test]
  fn test_shifts() {
    let i = 227;
    let j: Mpz = from_int(i);
    assert((i << 4).to_str() == (j << 4).to_str());
    assert((-i << 4).to_str() == (-j << 4).to_str());
    assert((i >> 4).to_str() == (j >> 4).to_str());
    assert((-i >> 4).to_str() == (-j >> 4).to_str());
  }

  #[test]
  fn test_compl() {
    assert(from_int::<Mpz>(13).compl().to_str() == (!13).to_str());
    assert(from_int::<Mpz>(-442).compl().to_str() == (!-442).to_str());
  }

  #[test]
  fn test_popcount() {
    let a = option::unwrap(from_str_radix("1010010011", 2));
    assert(a.popcount() == 5);
  }

  #[test]
  fn test_hamdist() {
    let a: Mpz = from_int(0b1011_0001);
    let b: Mpz = from_int(0b0010_1011);
    assert(a.hamdist(&b) == 4);
  }

  #[test]
  fn test_bit_length() {
    assert(from_int::<Mpz>(0b1011_0000_0001_0000).bit_length() == 16);
    assert(from_int::<Mpz>(0b101).bit_length() == 3);
  }
}
