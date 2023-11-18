all: configure cope symlinks

configure:
	-mkdir -p ${PREFIX} 2>/dev/null
	cd ./src/cope/ && ./configure --prefix=${PREFIX}

install: symlinks
	-mkdir -p ${PREFIX}
	find . -type l|while read a;do cp --copy-contents -LR  "$$a" ${PREFIX}/;done

cope:
	cd ./src/cope/ && ${MAKE}
  

symlinks:
	-@ln -s ./src/cope/src/cope 2> /dev/null
	-@ln -s ./src/kmerfreq/kmerfreq 2> /dev/null

test:
	cd ./TEST && sh test.sh

clean:
	cd ./src/cope && ${MAKE} clean
	cd ./src/cope/src && ${MAKE} clean
	-rm cope kmerfreq

allclean: clean
	cd ./TEST && rm -r -v kmerfreq mode0 mode1 mode2
