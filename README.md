### 2-3-4 Way Problem

#### Uses an approach to aggregate rows into (mostly, idk not mathematically proved) unique aggregates, from which the problem can be solved using a conventional 3-4 Sum approach. Aggregates are the dot product of row and a coefficient set, currently: increasing prime set.


## Setup:
#### generate.c: Builds 3 files each with (N) rows with at least one set of rows matching specifications (K,L)

#### three_sum.c: Solves problem of size (N,K,L)

```
gcc -O3 -march=native -funroll-loops -ffast-math -flto generate.c -o generate
./generate 500, 10000, 10000

gcc -O3 -march=native -funroll-loops -ffast-math -flto three_sum.c -o solve
./solve 500, 10000, 10000

```
