all: install
	cd sample;make

clean:
	$(RM) *~ *.o;
	cd sample; make clean

install:
	cp spline.h $(HOME)/include/spline.h

uninstall:
	rm $(HOME)/include/spline.h
