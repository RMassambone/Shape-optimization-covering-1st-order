# Shape-optimization-covering-1st-order
Source codes related to the paper:  E. G. Birgin, A. Laurain, R. Massambone, and A. G. Santana, A shape  optimization approach to the problem of covering a two-dimensional  region with minimum-radius identical balls, SIAM Journal on Scientific  Computing 43, pp. A2047-A2078, 2021

Algencan 3.1.1 is required to run the code. Algencan 3.1.1 can be
downloaded in: https://www.ime.usp.br/~egbirgin/tango/. The
distribution includes instructions to generate the library
libalgencan.a.

Compilation instructions for the code in this folder can be found in
'runall' and 'runamerica' files. Both files require an environment
variable ALGENCAN pointing to the folder were Algencan was installed,
like, for example,

export ALGENCAN=$HOME/algencan-3.1.1

Files covering-a-single-ball.f90 and covering-two-tangent-balls.f90
correspond to two simple degenerate problems.

Files covering-star-testing-h.f90 and covering-twosq-testing-h.f90
correspond to the numerical experiments the evaluate the influence of
the discretization step h in the obtained cover and the numerical
performance of the method.

Files covering-america.f90, modamerica.f90, and geometry.f90
correspond the the problem of covering a sketch of the map of
America. Fire modamerica.f90 was borrowed from the Algencan 3.1.1
distribution; while file geometry.f90 correspond to the package
GEOMETRY and it was downloaded from:
https://people.sc.fsu.edu/~jburkardt/f_src/geometry/geometry.html.
File runamerica shows how to compile and run these code.

Files covering-disconnected.f90, covering-heart.f90,
covering-rings.f90, covering-soap.f90, covering-star.f90, and
covering-twosq.f90 correspond, respectively, to the covering of the
regions ``disconnected'', ``heart'', the 3 different ``ring regions'',
``soap'', ``peaked star'', and ``two squares''. Instructions to
compile and run them can be find in file runall.
