-- | ignore
import "../lib/github.com/jonesz/rand-linalg/qr/tsqr"
import "../lib/github.com/diku-dk/linalg/qr"

module QR_ECON = qr_econ f32
module GS = mk_gram_schmidt f32
module HR_BLOCKED = mk_block_householder f32

-- ==
-- entry: bench_qr_econ
-- compiled random input { [1000][20]f32  }
-- compiled random input { [10000][20]f32 }
entry bench_qr_econ A = QR_ECON.qr A

-- ==
-- entry: bench_qr_econ_batch
-- compiled random input { [100][1000][20]f32  }
-- compiled random input { [100][10000][20]f32 }
entry bench_qr_econ_batch As = map (QR_ECON.qr) As

-- ==
-- entry: bench_gs
-- compiled random input { [1000][20]f32  }
-- compiled random input { [10000][20]f32 }
entry bench_gs A = GS.qr A

-- ==
-- entry: bench_gs_batch
-- compiled random input { [100][1000][20]f32  }
-- compiled random input { [100][10000][20]f32 }
entry bench_gs_batch As = map (GS.qr) As

-- ==
-- entry: bench_block_house
-- compiled random input { [1000][20]f32  }
-- compiled random input { [10000][20]f32 }
entry bench_block_house A = (HR_BLOCKED.qr 10i64) A

-- ==
-- entry: bench_block_house_batch
-- compiled random input { [100][1000][20]f32  }
-- compiled random input { [100][10000][20]f32 }
entry bench_block_house_blocked As = map (HR_BLOCKED.qr 10i64) As
