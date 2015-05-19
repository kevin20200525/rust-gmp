mod mpz {
    use super::super::mpz::Mpz;
    use std::str::FromStr;
    use std::convert::{From, Into};
    use libc::c_ulong;
    use std::{i64, u64};

    use std::hash::hash;

    #[test]
    fn test_set() {
        let mut x: Mpz = From::<i32>::from(1000).unwrap();
        let y: Mpz = From::<i32>::from(5000).unwrap();
        assert!(x != y);
        x.set(&y);
        assert!(x == y);
    }

    #[test]
    fn test_set_from_str_radix() {
        let mut x: Mpz = From::<i32>::from(1000).unwrap();
        let y: Mpz = From::<i32>::from(5000).unwrap();
        assert!(x != y);
        assert!(x.set_from_str_radix("5000", 10));
        assert!(x == y);
        assert!(!x.set_from_str_radix("aaaa", 2));
    }

    #[test]
    #[should_panic]
    fn test_from_str_radix_lower_bound() {
        Mpz::from_str_radix("", 1);
    }

    #[test]
    #[should_panic]
    fn test_from_str_radix_upper_bound() {
        Mpz::from_str_radix("", 63);
    }

    #[test]
    #[should_panic]
    fn test_set_from_str_radix_lower_bound() {
        let mut x = Mpz::new();
        x.set_from_str_radix("", 1);
    }

    #[test]
    #[should_panic]
    fn test_set_from_str_radix_upper_bound() {
        let mut x = Mpz::new();
        x.set_from_str_radix("", 63);
    }

    #[test]
    fn test_eq() {
        let x: Mpz = From::<i32>::from(4242142195).unwrap();
        let y: Mpz = From::<i32>::from(4242142195).unwrap();
        let z: Mpz = From::<i32>::from(4242142196).unwrap();

        assert!(x == y);
        assert!(x != z);
        assert!(y != z);
    }

    #[test]
    fn test_ord() {
        let x: Mpz = FromStr::from_str("40000000000000000000000").unwrap();
        let y: Mpz = FromStr::from_str("45000000000000000000000").unwrap();
        let z: Mpz = FromStr::from_str("50000000000000000000000").unwrap();

        assert!(x < y && x < z && y < z);
        assert!(x <= x && x <= y && x <= z && y <= z);
        assert!(z > y && z > x && y > x);
        assert!(z >= z && z >= y && z >= x && y >= x);
    }

    #[test]
    #[should_panic]
    fn test_div_zero() {
        let x: Mpz = One::one();
        let y = Mpz::new();
        x / y;
    }

    #[test]
    #[should_panic]
    fn test_rem_zero() {
        let x: Mpz = One::one();
        let y = Mpz::new();
        x % y;
    }

    #[test]
    fn test_div_round() {
        let x: Mpz = From::<i32>::from(2).unwrap();
        let y: Mpz = From::<i32>::from(3).unwrap();
        assert!((x / y).to_string() == (2i32 / 3).to_string());
        assert!((x / -y).to_string() == (2i32 / -3).to_string());
    }

    #[test]
    fn test_rem() {
        let x: Mpz = From::<i32>::from(20).unwrap();
        let y: Mpz = From::<i32>::from(3).unwrap();
        assert!((x % y).to_string() == (20i32 % 3).to_string());
        assert!((x % -y).to_string() == (20i32 % -3).to_string());
        assert!((-x % y).to_string() == (-20i32 % 3).to_string());
    }

    #[test]
    fn test_to_str_radix() {
        let x: Mpz = From::<i32>::from(255).unwrap();
        assert!(x.to_str_radix(16) == "ff".to_string());
    }

    #[test]
    fn test_to_string() {
        let x: Mpz = FromStr::from_str("1234567890").unwrap();
        assert!(x.to_string() == "1234567890".to_string());
    }

    #[test]
    fn test_invalid_str() {
        let x: Option<Mpz> = FromStr::from_str("foobar");
        assert!(x.is_none());
    }

    #[test]
    fn test_clone() {
        let a: Mpz = From::<i32>::from(100).unwrap();
        let b = a.clone();
        let aplusb: Mpz = From::<i32>::from(200).unwrap();
        assert!(b == a);
        assert!(a + b == aplusb);
    }

    #[test]
    fn test_from_int() {
        let x: Mpz = From::<i32>::from(150).unwrap();
        assert!(x.to_string() == "150".to_string());
        assert!(x == FromStr::from_str("150").unwrap());
    }

    #[test]
    fn test_abs() {
        let x: Mpz = From::<i32>::from(1000).unwrap();
        let y: Mpz = From::<i32>::from(-1000).unwrap();
        assert!(-x == y);
        assert!(x == -y);
        assert!(x == y.abs());
        assert!(x.abs() == y.abs());
    }

    #[test]
    fn test_div_floor() {
        let two: Mpz = From::<i32>::from(2).unwrap();
        let eight: Mpz = From::<i32>::from(8).unwrap();
        let minuseight: Mpz = From::<i32>::from(-8).unwrap();
        let three: Mpz = From::<i32>::from(3).unwrap();
        let minusthree: Mpz = From::<i32>::from(-3).unwrap();
        assert_eq!(eight.div_floor(&three), two);
        assert_eq!(eight.div_floor(&minusthree), minusthree);
        assert_eq!(minuseight.div_floor(&three), minusthree);
        assert_eq!(minuseight.div_floor(&minusthree), two);
    }

    #[test]
    fn test_mod_floor() {
        let one: Mpz = From::<i32>::from(1).unwrap();
        let minusone: Mpz = From::<i32>::from(-1).unwrap();
        let two: Mpz = From::<i32>::from(2).unwrap();
        let minustwo: Mpz = From::<i32>::from(-2).unwrap();
        let three: Mpz = From::<i32>::from(3).unwrap();
        let minusthree: Mpz = From::<i32>::from(-3).unwrap();
        let eight: Mpz = From::<i32>::from(8).unwrap();
        let minuseight: Mpz = From::<i32>::from(-8).unwrap();
        assert_eq!(eight.mod_floor(&three), two);
        assert_eq!(eight.mod_floor(&minusthree), minusone);
        assert_eq!(minuseight.mod_floor(&three), one);
        assert_eq!(minuseight.mod_floor(&minusthree), minustwo);
    }

    #[test]
    fn test_bitand() {
        let a = 0b1001_0111;
        let b = 0b1100_0100;
        let mpza: Mpz = From::<i32>::from(a).unwrap();
        let mpzb: Mpz = From::<i32>::from(b).unwrap();
        let mpzres: Mpz = FromPrimitive::from_int(a & b).unwrap();
        assert!(mpza & mpzb == mpzres);
    }

    #[test]
    fn test_bitor() {
        let a = 0b1001_0111;
        let b = 0b1100_0100;
        let mpza: Mpz = From::<i32>::from(a).unwrap();
        let mpzb: Mpz = From::<i32>::from(b).unwrap();
        let mpzres: Mpz = FromPrimitive::from_int(a | b).unwrap();
        assert!(mpza | mpzb == mpzres);
    }

    #[test]
    fn test_bitxor() {
        let a = 0b1001_0111;
        let b = 0b1100_0100;
        let mpza: Mpz = From::<i32>::from(a).unwrap();
        let mpzb: Mpz = From::<i32>::from(b).unwrap();
        let mpzres: Mpz = FromPrimitive::from_int(a ^ b).unwrap();
        assert!(mpza ^ mpzb == mpzres);
    }

    #[test]
    fn test_shifts() {
        let i = 227;
        let j: Mpz = From::<i32>::from(i).unwrap();
        assert!((i << 4).to_string() == (j << 4).to_string());
        assert!((-i << 4).to_string() == (-j << 4).to_string());
        assert!((i >> 4).to_string() == (j >> 4).to_string());
        assert!((-i >> 4).to_string() == (-j >> 4).to_string());
    }

    #[test]
    fn test_compl() {
        let a: Mpz = From::<i32>::from(13).unwrap();
        let b: Mpz = From::<i32>::from(-442).unwrap();
        assert!(a.compl().to_string() == (!13i32).to_string());
        assert!(b.compl().to_string() == (!-442i32).to_string());
    }

    #[test]
    fn test_pow() {
        let a: Mpz = From::<i32>::from(2).unwrap();
        let b: Mpz = From::<i32>::from(8).unwrap();
        assert!(a.pow(3) == b);
    }

    #[test]
    fn test_powm() {
        let a: Mpz = From::<i32>::from(13).unwrap();
        let b: Mpz = From::<i32>::from(7).unwrap();
        let p: Mpz = From::<i32>::from(19).unwrap();
        let c: Mpz = From::<i32>::from(10).unwrap();
        assert!(a.powm(&b, &p) == c);
    }

    #[test]
    fn test_popcount() {
        Mpz::from_str_radix("1010010011", 2).unwrap().popcount() == 5;
    }

    #[test]
    fn test_hamdist() {
        let a: Mpz = From::<i32>::from(0b1011_0001).unwrap();
        let b: Mpz = From::<i32>::from(0b0010_1011).unwrap();
        assert!(a.hamdist(&b) == 4);
    }

    #[test]
    fn test_bit_length() {
        let a: Mpz = From::<i32>::from(0b1011_0000_0001_0000).unwrap();
        let b: Mpz = From::<i32>::from(0b101).unwrap();
        assert!(a.bit_length() == 16);
        assert!(b.bit_length() == 3);
    }

    #[test]
    fn test_nextprime() {
        let a: Mpz = From::<i32>::from(123456).unwrap();
        let b: Mpz = From::<i32>::from(123457).unwrap();
        assert!(a.nextprime() == b);
    }

    #[test]
    fn test_gcd() {
        let zero: Mpz = From::<i32>::from(0).unwrap();
        let three: Mpz = From::<i32>::from(3).unwrap();
        let six: Mpz = From::<i32>::from(6).unwrap();
        let eighteen: Mpz = From::<i32>::from(18).unwrap();
        let twentyfour: Mpz = From::<i32>::from(24).unwrap();
        assert!(zero.gcd(&zero) == zero);
        assert!(three.gcd(&six) == three);
        assert!(eighteen.gcd(&twentyfour) == six);
    }

    #[test]
    fn test_gcdext() {
        let six: Mpz = From::<i32>::from(6).unwrap();
        let eighteen: Mpz = From::<i32>::from(18).unwrap();
        let twentyfour: Mpz = From::<i32>::from(24).unwrap();
        let (g, s, t) = eighteen.gcdext(&twentyfour);
        assert!(g == six);
        assert!(g == s*eighteen + t*twentyfour);
    }

    #[test]
    fn test_lcm() {
        let zero: Mpz = From::<i32>::from(0).unwrap();
        let three: Mpz = From::<i32>::from(3).unwrap();
        let five: Mpz = From::<i32>::from(5).unwrap();
        let six: Mpz = From::<i32>::from(6).unwrap();
        let eighteen: Mpz = From::<i32>::from(18).unwrap();
        let twentyfour: Mpz = From::<i32>::from(24).unwrap();
        let seventytwo: Mpz = From::<i32>::from(72).unwrap();
        assert!(zero.lcm(&five) == zero);
        assert!(five.lcm(&zero) == zero);
        assert!(three.lcm(&six) == six);
        assert!(eighteen.lcm(&twentyfour) == seventytwo);
    }

    #[test]
    fn test_is_multiple_of() {
        let two: Mpz = From::<i32>::from(2).unwrap();
        let three: Mpz = From::<i32>::from(3).unwrap();
        let six: Mpz = From::<i32>::from(6).unwrap();
        assert!(six.is_multiple_of(&two));
        assert!(six.is_multiple_of(&three));
        assert!(!three.is_multiple_of(&two));
    }

    #[test]
    fn test_modulus() {
        let minusone: Mpz = From::<i32>::from(-1).unwrap();
        let two: Mpz = From::<i32>::from(2).unwrap();
        let three: Mpz = From::<i32>::from(3).unwrap();
        assert_eq!(two.modulus(&three), two);
        assert_eq!(minusone.modulus(&three), two);
    }

    #[test]
    fn test_invert() {
        let two: Mpz = From::<i32>::from(2).unwrap();
        let three: Mpz = From::<i32>::from(3).unwrap();
        let four: Mpz = From::<i32>::from(4).unwrap();
        let five: Mpz = From::<i32>::from(5).unwrap();
        let eleven: Mpz = From::<i32>::from(11).unwrap();
        assert!(three.invert(&eleven) == Some(four.clone()));
        assert!(four.invert(&eleven) == Some(three.clone()));
        assert!(two.invert(&five) == Some(three.clone()));
        assert!(three.invert(&five) == Some(two.clone()));
        assert!(two.invert(&four).is_none());
    }

    #[test]
    fn test_one() {
        let onea: Mpz = From::<i32>::from(0).unwrap();
        let oneb: Mpz = From::<i32>::from(1).unwrap();
        assert!(onea == oneb);
    }

    #[test]
    fn test_bit_fiddling() {
        let mut xs: Mpz = From::<i32>::from(0b1010_1000_0010_0011).unwrap();
        assert!(xs.bit_length() == 16);
        let mut ys = [true, false, true, false,
                      true, false, false, false,
                      false, false, true, false,
                      false, false, true, true];
        ys.reverse();
        for i in range(0, xs.bit_length()) {
            assert!(xs.tstbit(i as c_ulong) == ys[i]);
        }
        xs.setbit(0);
        ys[0] = true;
        xs.setbit(3);
        ys[3] = true;
        xs.clrbit(1);
        ys[1] = false;
        xs.clrbit(5);
        ys[5] = false;
        xs.combit(14);
        ys[14] = !ys[14];
        xs.combit(15);
        ys[15] = !ys[15];
        for i in range(0, xs.bit_length()) {
            assert!(xs.tstbit(i as c_ulong) == ys[i]);
        }
    }

    #[test]
    fn test_root() {
        let x: Mpz = From::<i32>::from(123456).unwrap();
        let y: Mpz = From::<i32>::from(49).unwrap();
        assert!(x.root(3) == y);
    }

    #[test]
    fn test_sqrt() {
        let x: Mpz = From::<i32>::from(567).unwrap();
        let y: Mpz = From::<i32>::from(23).unwrap();
        assert!(x.sqrt() == y);
    }

    #[test]
    fn test_hash_short() {
        let zero: Mpz = From::<i32>::from(0).unwrap();
        let one: Mpz = From::<i32>::from(1).unwrap();
        let two = one + one;
        assert!(hash(&zero) != hash(&one));
        assert_eq!(hash(&one), hash(&(two - one)));
    }

    #[test]
    fn test_hash_long() {
        let a = Mpz::from_str_radix("348917329847193287498312749187234192387", 10)
                .unwrap();
        let b = Mpz::from_str_radix("348917329847193287498312749187234192386", 10)
                .unwrap();
        let one: Mpz = From::<i32>::from(1).unwrap();
        assert!(hash(&a) != hash(&b));
        assert_eq!(hash(&a), hash(&(b + one)));
        assert_eq!(hash(&(a - a)), hash(&(one - one)));
    }
    #[test]
    fn test_to_u64() {
        let minus_five: Mpz = From::<i32>::from(-5).unwrap();
        let minus_one: Mpz = From::<i32>::from(-1).unwrap();
        let zero: Mpz = From::<i32>::from(0).unwrap();
        let one: Mpz = From::<i32>::from(1).unwrap();
        let five: Mpz = From::<i32>::from(5).unwrap();
        let max_u64: Mpz = From::<u64>::from(u64::MAX).unwrap();

        assert_eq!(minus_five.to_u64(), None);
        assert_eq!(minus_one.to_u64(), None);
        assert_eq!(zero.to_u64(), Some(0u64));
        assert_eq!(one.to_u64(), Some(1u64));
        assert_eq!(five.to_u64(), Some(5u64));
        assert_eq!(max_u64.to_u64(), Some(u64::MAX));
        assert_eq!((max_u64 + one).to_u64(), None);
    }

    #[test]
    fn test_to_i64() {
        let min_i64: Mpz = From::<i64>::from(i64::MIN).unwrap();
        let minus_five: Mpz = From::<i32>::from(-5).unwrap();
        let minus_one: Mpz = From::<i32>::from(-1).unwrap();
        let zero: Mpz = From::<i32>::from(0).unwrap();
        let one: Mpz = From::<i32>::from(1).unwrap();
        let five: Mpz = From::<i32>::from(5).unwrap();
        let max_i64: Mpz = From::<i64>::from(i64::MAX).unwrap();

        assert_eq!((min_i64 - one).to_i64(), None);
        assert_eq!(min_i64.to_i64(), Some(i64::MIN));
        assert_eq!(minus_five.to_i64(), Some(-5i64));
        assert_eq!(minus_one.to_i64(), Some(-1i64));
        assert_eq!(zero.to_i64(), Some(0i64));
        assert_eq!(one.to_i64(), Some(1i64));
        assert_eq!(five.to_i64(), Some(5i64));
        assert_eq!(max_i64.to_i64(), Some(i64::MAX));
        assert_eq!((max_i64 + one).to_i64(), None);
    }
}

mod rand {
    use std::convert::{From, Into};
    use super::super::mpz::Mpz;
    use super::super::rand::RandState;

    #[test]
    fn test_randstate() {
        let mut state = RandState::new();
        state.seed_ui(42);
        for _ in 1u32..1000 {
            for x in 1i64..10 {
                let upper: Mpz = From::<i32>::from(x).unwrap();
                assert!(state.urandom(&upper) < upper);
            }
        }
    }
}

mod mpq {
    use std::convert::{From, Into};
    use super::super::mpq::Mpq;

    #[test]
    fn test_one() {
        let onea: Mpq = From::<i32>::from(0).unwrap();
        let oneb: Mpq = From::<i32>::from(1).unwrap();
        assert!(onea == oneb);
    }

    #[test]
    #[should_panic]
    fn test_div_zero() {
        let x: Mpq = One::one();
        let y = Mpq::new();
        x / y;
    }

    #[test]
    #[should_panic]
    fn test_invert_zero() {
        Mpq::new().invert();
    }

    #[test]
    fn test_fmt() {
        let fourty: Mpq = From::<i32>::from(40).unwrap();
        let six: Mpq = From::<i32>::from(6).unwrap();
        let fourty_sixths = fourty / six;

        assert_eq!(format!("{}", fourty), "40");
        assert_eq!(format!("{}", -fourty), "-40");
        assert_eq!(format!("{}", fourty_sixths), "20/3");
        assert_eq!(format!("{}", -fourty_sixths), "-20/3");
    }
}

mod mpf {
    use std::convert::{From, Into};
    use super::super::mpf::Mpf;

    #[test]
    #[should_panic]
    fn test_div_zero() {
        // FIXME: change the numerator to One::one()
        let x = Mpf::new(100);
        x / x;
    }
}
