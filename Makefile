.PHONY: all install uninstall clean

all: install

install: clean
	pip install .

uninstall: clean
	pip uninstall e2epy

clean:
	$(RM) -rf build e2epy.egg-info
