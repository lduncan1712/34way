### 2-3-4 Way Problem

#### Using an approach to Aggreage rows into (mostly) unique aggregates, from which the problem can be solved using more conventional 2-3-4 Sum


## Setup:
#### Generate.c: Builds 3 files each with (N) rows with at least one set of rows matching specifications (K,L)

### three_sum.c: Solves problem of size (N,K,L)

```
gcc -O3 -march=native -funroll-loops -ffast-math -flto generate.c -o generate
./generate 500, 10000, 10000

gcc -O3 -march=native -funroll-loops -ffast-math -flto three_sum.c -o solve
./solve 500, 10000, 10000

```
