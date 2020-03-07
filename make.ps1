cd rk4
gcc -c -fPIC satlut.c -o satlut.o
gcc -c -fPIC rk4.c -o rk4.o
gcc -shared rk4.o -o librk4.so -lsatlut
cd ..