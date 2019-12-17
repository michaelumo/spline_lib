all: install
	cd sample;make

clean:
	$(RM) *~ *.o *.csv
	cd sample; make clean

install:
	cp spline.h $(HOME)/include/

uninstall:
	rm $(HOME)/include/spline.h
