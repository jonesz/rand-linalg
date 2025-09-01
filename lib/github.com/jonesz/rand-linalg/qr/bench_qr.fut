-- | ignore
import "qr"

module QR = mk_householder_thin_qr f32

-- ==
-- entry: bench_householder_thin_qr_20
-- compiled random input { [100][20]f32  }
-- compiled random input { [1000][20]f32 }
-- compiled random input { [10000][20]f32 }
entry bench_householder_thin_qr_20 A = QR.qr () A

-- ==
-- entry: bench_householder_thin_qr_50
-- compiled random input { [100][50]f32  }
-- compiled random input { [1000][50]f32 }
-- compiled random input { [10000][50]f32 }
entry bench_householder_thin_qr_50 A = QR.qr () A
