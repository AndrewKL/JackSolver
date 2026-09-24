# JackSolver

JackSolver (`JackRHFSolver`) is a small **Restricted Hartree–Fock (RHF)** quantum chemistry solver written in Java. Given the positions and atomic numbers of a set of atoms and a total electron count, it runs a self-consistent field (SCF) calculation and reports the electronic, nuclear repulsion, and total energies of the molecule.

It is a Java port and extension of the Fortran IV program in Appendix B of *Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory* by Attila Szabo and Neil S. Ostlund. The numbered `step #N` comments in `JackRHFSolver.HFSolver` follow the SCF procedure described on page 146 of that book. If you are learning this material, the book is well worth reading alongside the code.

Written by Andrew Long, advised by Jason A. C. Clyburne, with additional help from Cory Pye.

## What it does

1. Builds a contracted Gaussian (STO-nG, default STO-3G) 1s orbital centered on each atom.
2. Computes the one-electron integrals: the overlap matrix **S**, kinetic energy **T**, and nuclear attraction **V<sub>nuc</sub>**, and forms **H<sub>core</sub> = T + V<sub>nuc</sub>**.
3. Diagonalizes **S** and builds the transformation matrix **X** by canonical orthogonalization (**X = U s<sup>-1/2</sup>**).
4. Precomputes every two-electron (four-center) integral `(μν|λσ)`, running in parallel on a thread pool sized to the number of CPU cores.
5. Starts from a zero density matrix **P** (the core-Hamiltonian guess) and iterates the SCF loop:
   - builds **G** from **P** and the two-electron integrals, then the Fock matrix **F = H<sub>core</sub> + G**
   - transforms **F' = X<sup>†</sup> F X** and diagonalizes it to get **C'** and the orbital energies ε
   - back-transforms **C = X C'** and forms a new density matrix, mixing 90% new with 10% old to damp oscillation
   - stops when the RMS change in **P** falls below `1e-4`, or after 500 iterations
6. Prints the electronic energy, the nuclear repulsion energy, the total energy, and (at verbose level 2) the Mulliken population matrix **PS**.

All quantities are in **atomic units**: coordinates in bohr, energies in hartree.

## Project layout

```
src/main/java/jackSolver/
  JackRHFSolver.java     SCF driver, main(), and matrix helpers
  Integrals.java         Overlap, kinetic, nuclear-attraction and 4-center integrals, F0 Boys function
  GaussianOrbital.java   Contracted Gaussian orbitals, STO-nG data, geometry helpers
  Atom.java              Atom position + atomic number (Z)
  Vec3D.java             Simple 3D vector
  ErrorFunction.java     erf() approximation (used by the F0 Boys function)
  BasisSet.java, Basis.java, STO3G.java
                         Early, unused scaffolding for more general basis sets
src/test/java/jackSolver/
  JackRHFSolverTest.java, GaussianOrbitalTest.java   JUnit 4 tests
build.gradle.kts         Gradle build (Java 21 toolchain)
gradlew, gradle/         Gradle wrapper
```

The solver uses the [CERN Colt](https://dst.lbl.gov/ACSSoftware/colt/) library (`DenseDoubleMatrix2D`, `Algebra`, `EigenvalueDecomposition`) for its linear algebra. Gradle downloads Colt from Maven Central (`colt:colt:1.2.0`).

## Requirements

- **JDK 21.** If your default `java` is a different version, point `JAVA_HOME` at a JDK 21 install. For example, with Homebrew on macOS: `export JAVA_HOME=/opt/homebrew/opt/openjdk@21`.
- You don't need to install Gradle. The included `./gradlew` wrapper downloads the right version.

## Building and running

```sh
./gradlew build      # compile and run the tests
./gradlew run        # run the built-in example (HeH+ at R = 1.4632 bohr)
./gradlew test       # run only the tests
```

On Windows, use `gradlew.bat`.

The example in `main` reproduces the Szabo & Ostlund HeH⁺ calculation. It should print a total energy of about **−2.8607 hartree**.

To get a standalone distribution with launch scripts, run `./gradlew installDist` and then `build/install/JackSolver/bin/JackSolver`.

## Using the solver from code

```java
import jackSolver.Atom;
import jackSolver.JackRHFSolver;

// Atom(x, y, z, atomicNumber) with coordinates in bohr
Atom[] h2 = {
    new Atom(0, 0, 0,      1),
    new Atom(0, 0, 1.4632, 1),
};

// HFSolver(atoms, numberOfElectrons, verbose)
//   verbose 0 = quiet, 1 = progress info, 2 = also print every matrix
double totalEnergy = JackRHFSolver.HFSolver(h2, 2, 1);   // ≈ -1.114 hartree
```

`HFSolver` returns the total energy (electronic plus nuclear repulsion). It returns `0` immediately if the electron count is odd, because restricted HF needs every electron paired.

## Tests

`JackRHFSolverTest` and `GaussianOrbitalTest` are JUnit 4 tests. They cover the individual pieces (integrals, orthogonalization, eigenvalue sorting, density matrix, convergence, nuclear repulsion) and complete runs checked against known energies:

| Test | System | Expected total energy (hartree) |
|------|--------|----------------------------------|
| `testHFSolver` | HeH⁺, R = 1.4632 bohr | −2.86066 |
| `test2HFSolver` | H₂, R = 1.4632 bohr | −1.11401 |
| `testH2GasTestHFSolver` | 8 widely spaced H₂ molecules + 1 He (stress test) | printed only |

```sh
./gradlew test
```

Some tests only print output and don't assert anything, such as `test3HFSolver`, `testH2GasTestHFSolver` and `testMemory`.

## Limitations

This is an educational and research prototype, not a general-purpose quantum chemistry package:

- **Only H (Z = 1) and He (Z = 2) are supported.** The STO-nG orbital exponents are scaled by Slater ζ values that exist only for these two elements (ζ = 1.24 for H and 2.0925 for He, as in Szabo & Ostlund). Any other atomic number will throw an `ArrayIndexOutOfBoundsException`.
- **One s-type (1s) basis function per atom.** There are no p/d functions, so only the F<sub>0</sub> Boys function is implemented. `STO3G.java` is an unfinished step toward carbon and larger basis sets.
- **Closed-shell only** (restricted HF, even electron count).
- The four-center integrals are stored as a full K⁴ array with no symmetry reduction, so memory and time grow quickly with the number of atoms. The parallel integral step gives up waiting after 5 minutes.
- `erf` uses a Numerical Recipes Chebyshev approximation, which has a fractional error of about 1.2 × 10⁻⁷.
- `Integrals` and `JackRHFSolver` print progress output even when `verbose = 0`.

## Credits and licenses

- The algorithm and reference data come from Szabo & Ostlund, *Modern Quantum Chemistry* (Dover).
- **Colt** is © CERN and is distributed under its own license.
- `ErrorFunction.java` is adapted from Sedgewick & Wayne's *Introduction to Programming in Java* sample code.
