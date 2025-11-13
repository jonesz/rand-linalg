.PHONY: test
test:
	$(MAKE) -C lib/github.com/jonesz/rand-linalg/test

.PHONY: doc
doc:
	futhark doc -o doc/ lib/github.com/jonesz/rand-linalg

.PHONY: clean
clean:
	$(MAKE) -C lib/github.com/jonesz/rand-linalg/test clean
	$(RM) -rf doc/
