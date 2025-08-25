all: bench_qr

.PHONY: bench_qr clean

bench_qr:
	futhark bench --backend=multicore qr.fut

clean:
	$(RM) *.c *.expected *.actual
	$(RM) qr
