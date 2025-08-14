import "../cbrng-fut/cbrng"
import "../cbrng-fut/distribution"

module type sketch = {
	module dist: cbrng_distribution
	val sketch : dist.distribution -> (n: i64) -> (b: i64) -> [n][b]dist.num.t
}

local module mk_sketch (D: cbrng_distribution): sketch = {
	module dist = D
	def sketch param n b =
		tabulate (n * b) (D.rand param) |> unflatten
}

module sketch_rademacher (R: real) = mk_sketch (rademacher_distribution R u32 i64 squares32)
module sketch_gaussian (R: real) = mk_sketch (gaussian_distribution R u32 i64 squares32)
