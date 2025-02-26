### 2-3-4 Way Problem

#### Using an approach to Aggreage rows into (mostly) unique aggregates, from which the problem can be solved using more conventional 2-3-4 Sum


## Setup:
### This program includes 2 C files, one to generate test files ("A.txt",..."C.txt") with at least on valid solution, as well as a program to identify these rows (both take arguments of K, L (lambda), N)

```
gcc -O3 -march=native -funroll-loops -ffast-math -flto generate.c -o generate
./generate 500, 10000, 10000

gcc -O3 -march=native -funroll-loops -ffast-math -flto three_sum.c -o solve
./solve 500, 10000, 10000

```
