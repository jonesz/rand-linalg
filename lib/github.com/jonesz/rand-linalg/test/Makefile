.PHONY: test clean 

test: test_sketch_equivalence test_ste_eps test_ste_chebyshev

test_sketch_equivalence: sketch_equivalence.fut
	futhark test $^

test_ste_eps: ste_eps.fut
	futhark test $^

test_ste_chebyshev: ste_chebyshev.fut
	futhark test $^

clean: rm_sketch_equivalence rm_ste_eps rm_ste_chebyshev

rm_sketch_equivalence:
	$(RM) sketch_equivalence sketch_equivalence.c sketch_equivalence.fut.test_*

rm_ste_eps:
	$(RM) ste_eps ste_eps.c ste_eps.fut.test_*

rm_ste_chebyshev:
	$(RM) ste_chebyshev ste_chebyshev.c she_chebyshev.fut.test_*
