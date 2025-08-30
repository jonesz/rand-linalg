-- | ignore
import "tsqr"

module QR = qr_econ f32

-- ==
-- entry: bench_qr_econ_20
-- compiled random input { [100][20]f32  }
-- compiled random input { [1000][20]f32 }
-- compiled random input { [10000][20]f32 }
entry bench_qr_econ_20 A = QR.qr A

-- ==
-- entry: bench_qr_econ_50
-- compiled random input { [100][50]f32  }
-- compiled random input { [1000][50]f32 }
-- compiled random input { [10000][50]f32 }
entry bench_qr_econ_50 A = QR.qr A
