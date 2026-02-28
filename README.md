# DISS-Implementation
Implementation of the Dual inverse scaling and squaring algorithm for the matrix logarithm, as introduced by M. Fasi and B. Iannazzo.

Code for the seminar report for the TU Berlin Seminar Numerical Linear Algebra in the Winter semester 25/26.

# Usage
- The logarithm is evaluated by calling `diss.m` on a square matrix.
- Measurements evaluting accuracy and time, and comparing it to matlabs `logm`, can be generated with `generate_measurements.m`.
- Reference solutions used in that evaluation are generated with `create_reference_solutions.m`, this however requires that the Advanpix multiprecision computation toolbox is installed. This step can be skipped, as the pre-generated reference solutions are included in the `data/` folder.
- The plots used in the report can be generated with `generate_plots.m`
