## 2024-05-23 - Julia Small Vector Optimization
**Learning:** Broadcasting (`@.`) in Julia has measurable overhead for very small vectors (e.g., length 6) compared to explicit loops.
**Action:** For tight inner loops with small, fixed-size arrays, prefer explicit loops (or StaticArrays.jl if dependencies allow) over broadcasting.
