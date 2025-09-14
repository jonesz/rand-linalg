.PHONY: doc clean test

all: test

test:
	$(MAKE) -C lib/github.com/jonesz/rand-linalg/qr
	$(MAKE) -C lib/github.com/jonesz/rand-linalg/test

doc:
	futhark doc -o doc/ lib/github.com/jonesz/rand-linalg

clean:
	$(MAKE) -C lib/github.com/jonesz/rand-linalg/qr clean
	$(MAKE) -C lib/github.com/jonesz/rand-linalg/test clean
	$(RM) -rf doc/
