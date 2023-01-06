# Chapter2: Case study

## The shortcut problem

- v0: a baseline impletation

- v0-omp: add omp to v0

- v1: define `t` as the transpose of `d` (I think it should be done by math lib? )

- v1-omp: add omp to v1

- v2: Instruction-level parallelism

    for `min()` [split in 4 in `v2.cc`]

    ```c++
    v0=std::min(v0,z0);
    v1=std::min(v1,z1);
    ```

- v3: vectorization

- v4: data reuse