extern mod std;

use libc::c_char;
use libc::c_int;
use libc::size_t;
use libc::c_void;
use ptr::null;
use str::as_c_str;
use ptr::addr_of;

struct mpz_struct {
  _mp_alloc: c_int,
  _mp_size: c_int,
  _mp_d: *c_void
}

type mpz_t = *mpz_struct;

extern mod gmp {
  fn __gmpz_init(x: mpz_t);
  fn __gmpz_clear(x: mpz_t);
  fn __gmpz_set_str(rop: mpz_t, str: *c_char, base: c_int) -> c_int;
  fn __gmpz_get_str(str: *c_char, base: c_int, op: mpz_t) -> *c_char;
  fn __gmpz_sizeinbase(op: mpz_t, base: c_int) -> size_t;
  fn __gmpz_cmp(op: mpz_t, op2: mpz_t) -> c_int;
}

use gmp::*;

pub struct Mpz {
  priv mpz: mpz_struct,

  drop {
    __gmpz_clear(addr_of(&self.mpz));
  }
}

impl Mpz {
  fn set_str(&self, s: &str, base: int) -> bool {
    let mpz = addr_of(&self.mpz);
    let r = as_c_str(s, { |s| __gmpz_set_str(mpz, s, base as c_int) });
    r == 0
  }

  fn size_in_base(&self, base: int) -> uint {
    __gmpz_sizeinbase(addr_of(&self.mpz), base as c_int) as uint
  }
}

impl Mpz: cmp::Eq {
  pure fn eq(other: &Mpz) -> bool unsafe {
    __gmpz_cmp(addr_of(&self.mpz), addr_of(&other.mpz)) == 0
  }
  pure fn ne(other: &Mpz) -> bool unsafe {
    __gmpz_cmp(addr_of(&self.mpz), addr_of(&other.mpz)) != 0
  }
}

impl Mpz: cmp::Ord {
  pure fn lt(other: &Mpz) -> bool unsafe {
    __gmpz_cmp(addr_of(&self.mpz), addr_of(&other.mpz)) < 0
  }
  pure fn le(other: &Mpz) -> bool unsafe {
    __gmpz_cmp(addr_of(&self.mpz), addr_of(&other.mpz)) <= 0
  }
  pure fn gt(other: &Mpz) -> bool unsafe {
    __gmpz_cmp(addr_of(&self.mpz), addr_of(&other.mpz)) > 0
  }
  pure fn ge(other: &Mpz) -> bool unsafe {
    __gmpz_cmp(addr_of(&self.mpz), addr_of(&other.mpz)) >= 0
  }
}

fn init() -> Mpz {
  let mpz = mpz_struct { _mp_alloc: 0, _mp_size: 0, _mp_d: null() };
  __gmpz_init(addr_of(&mpz));
  Mpz { mpz: mpz }
}

#[cfg(test)]
mod tests {
  #[test]
  fn size_in_base() {
    let x = init();
    x.set_str("150000", 10);
    assert(x.size_in_base(10) == 6);
  }
  #[test]
  fn eq() {
    let x = init();
    x.set_str("4242142195", 10);
    let y = init();
    y.set_str("4242142195", 10);
    let z = init();
    z.set_str("4242142196", 10);

    assert(x == y);
    assert(x != z);
    assert(y != z);
  }
  #[test]
  fn ord() {
    let x = init();
    x.set_str("40000000000000000000000", 10);
    let y = init();
    y.set_str("45000000000000000000000", 10);
    let z = init();
    z.set_str("50000000000000000000000", 10);

    assert(x < y && x < z && y < z);
    assert(x <= x && x <= y && x <= z && y <= z);
    assert(z > y && z > x && y > x);
    assert(z >= z && z >= y && z >= x && y >= x);
  }
}
