-- | ignore

import "dist"
import "cbrng"

-- Straight port from `diku-dk/cpprandom/random_tests.fut`.
module mktest (dist: cbrng_distribution) = {
  module engine = dist.engine
  module num = dist.num

  def test (n: i64) (cfg: dist.configuration) =
    let xs = map (dist.rand cfg) (iota n)
    let mean = num.(reduce (+) (i32 0) xs / i64 n)
    in mean
}

-- Straight port from `diku-dk/cpprandom/random_tests.fut`.
module mktest_f (dist: cbrng_distribution) (R: real with t = dist.num.t) = {
  module engine = dist.engine
  module num = dist.num

  def test (n: i64) (cfg: dist.configuration) =
    let xs = map (dist.rand cfg) (iota n)
    let mean = num.(reduce (+) (i32 0) xs / i64 n)
    let stddev =
      R.(xs
         |> map (\x -> (x - mean))
         |> map num.((** i32 2))
         |> sum
         |> (/ i64 n)
         |> sqrt)
    in (R.round mean, R.round stddev)
}

-- ==
-- entry: test_rademacher
-- compiled random input { 100i64 i64 }   output { 0.0_f32 1.0_f32 }
-- compiled random input { 1_000i64 i64 } output { 0.0_f32 1.0_f32 }

module test_rademacher_f =
  mktest_f (rademacher_distribution f32 u32 i64 squares32) f32

entry test_rademacher n k =
  test_rademacher_f.test n k

-- ==
-- entry: test_normal
-- compiled random input { 100i64 i64 }   output { 0.0_f32 1.0_f32 }
-- compiled random input { 1_000i64 i64 } output { 0.0_f32 1.0_f32 }

module test_normal_f =
  mktest_f (normal_distribution f32 u32 i64 squares32) f32

entry test_normal n k =
  test_normal_f.test n (k, {mean = 0_f32, stddev = 1_f32})

-- ==
-- entry: test_uniform_real
-- compiled random input { 10_00i64 i64 }   output { 50.0_f32 2.0_f32 }
-- compiled random input { 10_000i64 i64 } output { 50.0_f32 2.0_f32 }

module test_uniform_real_f =
  mktest_f (uniform_real_distribution f32 u32 i64 squares32) f32

entry test_uniform_real n k =
  -- TODO: Why is this (1, 100) and not (0, 100); taken from the diku-dk code?
  test_uniform_real_f.test n (k, (1.0_f32, 100.0_f32))

-- ==
-- entry: test_uniform_int
-- compiled random input { 100i64   i64 } output { 50i64 }
-- compiled random input { 1_000i64 i64 } output { 50i64 }

module test_uniform_int_i =
  mktest (uniform_int_distribution i64 u32 i64 squares32)

entry test_uniform_int n k =
  -- TODO: Why is this (1, 100) and not (0, 100); taken from the diku-dk code?
  test_uniform_int_i.test n (k, (1i64, 100i64))
