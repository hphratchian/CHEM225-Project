# CHEM225-Project

### Updated 5/8/19 



#### Phase 1 | Write a program that prints the alpha and beta string determinants for singles substitutions.

##### Additional Items to Code

- [x] Add code to read nAlpha, nBeta, and nBasis directly from command line. (Done 4/13/19)
- [x] Format the output for alpha and beta strings (Done 4/13/19)
- [x] Fix the conditional statement that requires nBasis to be greater than either nAlpha or nBeta. (Done 4/12/19)

#### Phase 2 | Modify the program so that it uses the array containing nDet # of singles substitutions as the reference for the doubles substitutions. 

##### Additional Items to Code

- [x] Create another subroutine for doubles. (Done 5/7/19)

Compile using: pgfortran -i8 -o filename.exe filename.f90

-i8 Makes default integer and logical variables 8 bytes long (same as the -integer_size  64  option). The default is -integer_size 32.

###### Problems: 
-Segmentation Fault : occurs when nBasis is significantly larger than either nAlpha or nBeta, works up to /.exe 65 65 66

#### Phase 3 | Form Integrals Based on SC Rules

- [ ] Pass IStrings : Find a way to form all combinations of possible matrix elements using the reference, singles, and doubles determinants. Figure out how to pass the Istrings array from singles and doubles subroutines into allocatable arrays defined in the program. Loop over all iDet. (Not Done/In Progress)
- [ ] XOR : Check for differences in orbital occupations. Not done bitwise, using integer strings. (Not Done/In Progress)
- [ ] Spin Symmetry : Spin blocks and spin flips
#### Phase 4 | Solve S^2 Matrix elements Using MO overlap Matrices

-[ ] Generate MO coefficients and get overlaps. Compute using Gaussian--H2 minimal basis @ 2 Å (Not Done)
-[ ] Get nAlpha, nBeta, nBasis, and form overlaps to solve $S^2$ matrix elements
-[ ] Diagonalize and get $\langle S^2 \rangle$ values.

#### Phase 5 | Create Spin Projector
  .

  .

  .

  .

  .
