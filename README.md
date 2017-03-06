# uncmnd
UNCMND - Multi-dimensional Optimization Routine.  Program for minimizing a function of several variables (with no constraints.)

This program is a translation of the FORTRAN routine UNCMND written by Stephen Nash, George Mason University.
From the book "Numerical Methods and Software"
D. Kahaner, C. Moler, and S. Nash
Prentice Hall, 1988

To distinguish this program from others, an '_ak1' suffix has been appended to its name.

The underlying multi-dimensional optimization routine, UNCMND, is incorporated into the posted program to solve the General Rosenbrock
Function for N = 10. It has been tested with several other cases, and small code blocks used in these test cases have been left in the
program, but commented out. For example, for an Objective Function of the form F = A xy + B x^2 y + C x y^2. In such a case,
the coefficients A, B, and C would be saved in the coefficient matrix, which would then be passed into the Objective Function as
required. For a live, working, example of such a function, readers are invited to see a JavaScript version of this program that solves for four dimensions at the following URL:

http://www.akiti.ca/MultiOptRosen1.html

Several other changes have been made to the original version of UNCMND:

1) the options to supply analytic functions for the derivative array and Hessian matrix have been removed; instead, the derivative
array is computed by finite difference and the Hessian matrix is computed by secant updates. The user has no choice.

2) the error code info = 0 has been removed. I saw no way for the program to ever return with info having a value of 0.
I am certain this is a typo, so thought it best to not even include in the program the option info = 0.

3) A few boolean variables have been introduced to keep track of whether or not it is the first time through the program. If so,
some calculations can be avoided.
