.PHONY: test clean 

test: test_sketch_equivalence test_ste_eps test_ste_chebyshev test_distribution

test_sketch_equivalence: sketch_equivalence.fut
	futhark test $^

test_ste_eps: ste_eps.fut
	futhark test $^

test_ste_chebyshev: ste_chebyshev.fut
	futhark test $^

test_distribution: test_distribution.fut
	futhark test $^

test_svd: test_svd.fut
	futhark test $^

clean: rm_sketch_equivalence rm_ste_eps rm_ste_chebyshev rm_test_distribution rm_test_svd

rm_sketch_equivalence:
	$(RM) sketch_equivalence sketch_equivalence.c sketch_equivalence.fut.*actual sketch_equivalence.fut.*expected

rm_ste_eps:
	$(RM) ste_eps ste_eps.c ste_eps.fut.test_*

rm_ste_chebyshev:
	$(RM) ste_chebyshev ste_chebyshev.c ste_chebyshev.fut.test_*

rm_test_distribution:
	$(RM) test_distribution test_distribution.c test_distribution.fut.test_*

rm_test_svd:
	$(RM) test_svd test_svd.c test_svd.fut.test_*
