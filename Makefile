reedsol.o: reedsol.c Makefile
	gcc -O -c -o $@ $< -DLIB -I/usr/local/include -L/usr/local/lib
