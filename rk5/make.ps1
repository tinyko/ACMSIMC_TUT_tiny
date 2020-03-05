gcc -c -fPIC satlut.c -o satlut.o
gcc -c -fPIC rk5.c -o rk5.o
gcc -static -shared rk5.o -o librk5.so -lsatlut