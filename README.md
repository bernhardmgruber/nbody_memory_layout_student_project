# n-body memory layout exploration

This code demonstrates an n-body simulation using 3 different memory layouts to store particles.

A particle has the following floating-point properties:

 - position x
 - position y
 - position z
 - velocity x
 - velocity y
 - velocity z
 - mass

Each step of the simulation contains two subroutine, update and move.
Update iterates twice over all particles and computes their pairwise interactions, updating a particle's velocity based on the influence of other particles.
Move iterates once over all particles and updates a particle's position based on its velocity.

The same simulation is implemented using 3 different memory layouts:

 - AoS (Array of Structs)
 - SoA (Struct of Arrays)
 - AoSoA (Array of Structs of Arrays) with a compile time size for the inner array (`LANES` constant)


## Building and running

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
./nbody
```

Make sure you use g++ 9.3.0 or newer or clang++ 10 or newer.
The CMakeLists.txt asks the compiler to emit AVX512 instructions.
If that does not work on your machine, replace `-mavx512f` by `-march=native`.

## Tasks

 1. Compile and run the given sample code. Compare the runtime performance of the 3 different memory layouts.
 2. Have a quick look at the dissassembly: https://godbolt.org/z/1jqo3z.
    Search for `vsqrtss` and `vsqrtps`, which are a scalar and a vectorized square root instructions generated from the call to `sqrt(...)` in the program.
    Notice the suffixes `ss` and `ps` on the various instructions, where `ss` indicates scalar instructions and `ps` indicates vector instructions.
    Vector instructions run on multiple floats at the same time.
    They are crucial for performant CPU code and require a certain memory layout.
 3. Write CUDA versions of the given nbody simulations with different memory layouts inside your buffers.
    One CUDA thread would compute 1 particle, but you are free to try differently ;)
    Measure their runtime.
    Is there a similar runtime difference between the 3 memory layouts on the CPU version?
    What is the optimal value for `LANES` for the AoSoA?
 4. Use a profiler like `ncu` to profile your CUDA versions.
    Find explanations for the runtime differences of the memory layouts.
 5. Write a CUDA variant of the n-body that uses shared memory for caching.
    This version should be faster than the previous ones. Why?
    What is the optimal number of particles to store in shared memory?
 6. Does it matter if we even use different memory layouts (AoS, SoA or AoSoA) for the shared memory as well?
    Which combination of global and shared memory layout runs the fastest?
 7. Bonus: how can you make sure your simulation is still correct? Mistakes happen ;)
