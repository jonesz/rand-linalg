import "../lib/github.com/jonesz/rand-linalg/tsqr"
import "../lib/github.com/diku-dk/linalg/qr"

module GS = mk_gram_schmidt f32

-- ==
-- entry: bench_tsqr
-- compiled random input { 4i64 [400][10]f32  }
-- compiled random input { 4i64 [4000][10]f32 }
-- compiled random input { 8i64 [400][10]f32  }
-- compiled random input { 8i64 [4000][10]f32 }
-- compiled random input { 16i64 [400][10]f32  }
-- compiled random input { 16i64 [4000][10]f32 }
entry bench_tsqr k A = tsqr k A

-- ==
-- entry: bench_gs
-- compiled random input { [400][10]f32  }
-- compiled random input { [4000][10]f32 }
entry bench_gs A = GS.qr A
