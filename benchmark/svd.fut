-- | ignore
import "../lib/github.com/jonesz/rand-linalg/svd"

module SVD = mk_one_sided_jacobi_slow f32

-- ==
-- entry: bench_svd
-- compiled random input { [10][10]f32 }
-- compiled random input { [100][100]f32 }
-- compiled random input { [250][250]f32 }
entry bench_svd A = SVD.svd A
