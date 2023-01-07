# Multithreading with OpenMP

## Introduction

### parallel regions

aim: create multiple threads

![parallel regions](https://ppc.cs.aalto.fi/ch3/parallel_region.png)

```c++
a();
#pragma omp parallel
{
    c(1);
    c(2);
}
z();
```
result: all the threads do the same work

### critical sections

aim: at most one thread is executeing code that is inside the critical section

![critical sections](https://ppc.cs.aalto.fi/ch3/critical_section.png)

```c++
a();
#pragma omp parallel
{
    c(1);
    #pragma omp critical
    {
        c(2);
    }
    c(3);
    c(4);
}
z();
```

### shared / private data

- shared: variable that is declared outside a `parallel` region, one copy only

- private: variable that is declared inside a `parallel` region, each thread has its own copy

- if any thread ever writes to a shared variable, make sure it is safe

```c++
static void critical_example(int v) {

    // Shared variables
    int a = 0;
    int b = v;

    #pragma omp parallel
    {
        // Private variable - one for each thread
        int c;

        // Reading from "b" is safe: it is read-only
        // Writing to "c" is safe: it is private
        c = b;

        // Reading from "c" is safe: it is private
        // Writing to "c" is safe: it is private
        c = c * 10;

        #pragma omp critical
        {
            // Any modifications of "a" are safe:
            // we are inside a critical section
            a = a + c;
        }
    }

    // Reading from "a" is safe:
    // we are outside parallel region
    std::cout << a << std::endl;
}
```

result: `critical_example(1)` -> return 40

## for loops

![for loop](https://ppc.cs.aalto.fi/ch3/parallel_for_1.png)

```c++
a();
#pragma omp parallel for
for (int i = 0; i < 10; ++i) {
    c(i);
}
z();
```

equals to :

```c++
a();
#pragma omp parallel
{
    #pragma omp for
    for (int i = 0; i < 10; ++i) {
        c(i);
    }
}
z();
```

## nowait

aim: ignore synchronization

![nowait](https://ppc.cs.aalto.fi/ch3/parallel_for_nowait.png)

```c++
a();
#pragma omp parallel
{
    b();
    #pragma omp for nowait
    for (int i = 0; i < 10; ++i) {
        c(i);
    }
    d();
}
z();
```

## shedult

## nested