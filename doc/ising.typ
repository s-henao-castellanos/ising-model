#import "@preview/revtyp:0.15.0": revtable, revtyp

#show: revtyp.with(
  journal: "PRX",

  // Paper title
  title: [
    Simulation of the isotropic 2D Ising model
  ],

  // Author list
  authors: (
    (name: "Santiago Henao Castellanos", at: "uniandes", email: "s.henao10@uniandes.edu.co"),
  ),
  affiliations: (
    uniandes: [#link("https://uniandes.edu.co/")[Universidad de los Andes], Bogotá, Colombia],
  ),


  // Paper abstract
  //abstract: include "abstract.typ"
  abstract: [
    In this article we explain the constructions and results of a simulation of the 2D isotropic Ising system using the Metropolis Hastings Algorithm, implemented in the Python programming language. We use an heuristic function to determine sequence convergence, and we compare our results to the theoretical solution for energy, magnetization and heat capacity, showing coincidence within error intervals.
    The code to reproduce this article can be found in #link("https://github.com/s-henao-castellanos/ising-model")[this public Github repository].
  ],

  // Header, footer, etc.
  header: (
    //title: [PHYSICAL REVIEW ACCELERATORS AND BEAMS *00*, 000000 (0000)],
    left: (even: none, odd: none),
    right: (even: none, odd: none),
    //rule: false,
  ),
  footer: (
    title-left: none,
    title-right: none,
    center: [
      page #context counter(page).display()
      of #context counter(page).final().last()
    ],
  ),
  footnote-text: [
    Licensed under the terms of the
    #link(
      "https://creativecommons.org/licenses/by/4.0/",
    )[Creative Commons Attribution 4.0 International]
    license.
  ],
  wide-footnotes: true,

  // Writing utilities
  //show-line-numbers: true,
)


= Introduction

The Ising model @Ising1925 consists of stationary spins in a network interacting with their neighbours magnetically, with the Hamiltonian:
$
H = - sum_(chevron.l i,j chevron.r) J_(i,j) s_i s_j 
$
The spins $s_i$ can only point up or down ($s_i = plus.minus 1$), and the sum is taken over each spin $i$ and all of its first order neighbours. The interaction strength $J_(i,j)$ can be used to model anisotropy ---for instance, in a rectangular non-square lattice, but it is generally taken to be constant. We will set $J=1$. The intuition behind this system is that spins want to be parallel to each other, thus minimizing their energy, but a thermal bath will randomly flip them, causing disorder.


Historically, the model was given by Wilhelm Lenz to his student, Ernst Ising, in 1920, and he solved the one-dimensional case, where no phase transition occurs. A theoretical solution of the two-dimensional case was due to Onsager @Onsager1944, which proved that there is a phase transition when there are infinite spins, a result better understood thanks to the later work of Kramers, Wannier @Kramers1941 @Kramers1941b and Yang @Yang1952. A general solution with an external field is yet to be found in the 2D case, and the 3D case is been deemed untractable.

The simulation of this system is done via the Metropolis Hastings algorithm @Metropolis1953, using the canonical ensemble for the statistical description. We follow mostly the description of this algorithm presented in @NewmanBook.

Even with a finite number of spins, a considerable growth on the heat capacity of the lattice can be seen, at $
beta_C = 1/2 "arcsinh"(1) approx 0.4406867
$
We will explore the behaviour of the system in a range of temperatures that includes this critical point.

= Implementation details

== Mathematical assumptions

We assume $J=1$ and use $beta = 1/(k_B T)$ as our temperature parameter. We denote by $N$ the number of spins in each direction, so the total number of spins will be $N^2$, and use periodic boundary conditions: on the file $j$ $s_(N,j)$ is identified to $s_(0,j)$, and the same in the columns. 

We always start each spin as a random number $-1$ or $+1$ each with probability $1/2$. The default PRNG of Numpy is used.

As a detail for comparison, the theoretical behaviour of the system according to Onsager @Onsager1944 is the following
$
  m = M/N = (1-sinh(2 beta)^(-4))^(1/8) 
$
$
  epsilon = -1/tanh(2 beta) (1+(2/pi) (2 tanh(beta)^2-1) K(xi)
$
where $m$ is the mean magnetization, $epsilon$ the mean energy, $K(xi)$ the first complete elliptic integral and $xi=4( sinh(2beta))^2/(cosh(2beta)^4)$ its parameter.



== The code


As Python is a rather slow programming language, but it is exceedingly common in the physical sciences, we use the Numba package to speed up the calculations. This package compiles the numerical-heavy code into LLVM instructions, which gives considerable boosts of speed in comparison con python or even Numpy code. 

We calculate only the local energy and magnetization changes from interation to iteration in order to avoid unnecesary looping of the spins array.

For a given temperature $\beta$, the function `simulate_ising(grid,β,steps,thin)` simulates the model given an initial grid (random or not) for `steps` number of iterations, and thin the result by only saving the energy and magnetization for one of each `thin` steps. A taste of this procedure can be seen on @fig:iterations, where several simulations from random initial states are shown as semi-transparent lines, to highlight the stochasticity of the phenomenon. 

In principle, we should use an MCMC convergence test to detect dynamically whether a simulation has converged (i.e. reached thermal equilibrium) for the given temperature, but instead we opted for the simple approach of fixing the number of iterations. We can, in principle, spend the number of iterations in each temperature, but the systems takes more iterations to settle the clores we are to the critical temperature. Thus, we designed an empiric iteration function that can be seen in @fig:iteration-plan.

Additionaly, we used a $5sigma$ clipping procedure to calculate the statistics (mean and variance) in order to smoothen the stochastic nature of the simulations.


#figure(
  placement: top,
  image("img/iterations.pdf"),
  caption: [Convergence of mean energy and magnetization in the simulation. The temperature used for this visualization was $beta=0.8$.],
) <fig:iterations>



#figure(
  placement: top,
  image("img/iteration_plan.pdf"),
  caption: [Number of interations forced for each tempearture. ],
) <fig:iteration-plan>



= Results


We simulated 900 spins ($N=30$) in a set of 100 temperatures $beta$ ranging from $0.1$ to $0.9$, starting from a random disposition, and each temperature continuing from the last temperature final iteration in order to minimize the computation time. 

The results of the energy can be seen in @fig:E, where an excellent agreement can be seen. As the number of spins increases, there should be a discontinuity on the blue line, at the critical temperature.
The mean magnetization on @fig:m shows a more pronouced behaviour at the critical point, and also the comparison with the simulation points is more bleak. This is most likely caused by the behaviour seen in @fig:iterations, where the magnetization is notoriously more stochastic in ntature than the energy.

The main result is, however, the heat cpacity, shown in @fig:Cv, calculated as the variance of the energy, showing the characteristic lambda-like shape in the critical point, and showing remarkable agreement with the theory.



#figure(
  placement: top,
  image("img/E.pdf"),
  caption: [Energy of the simulation as a function of temperature. ],
) <fig:E>


#figure(
  placement: bottom,
  image("img/m.pdf"),
  caption: [Magnetization of the simulation as a function of temperature. ],
) <fig:m>



#figure(
  placement: top,
  image("img/Cv.pdf"),
  caption: [Heat cpacity of the simulation  as a function of temerature.],
) <fig:Cv>



/*















#figure(
  placement: auto,
  image("img/m.pdf"),
  caption: [
  ],
) <fig:m>




#figure(
  placement: auto,
  image("img/Cv.pdf"),
  caption: [
  ],
) <fig:Cv>

*/
= Conclusion

We sucessfully simulated the Ising model in a 2D isotropic lattice, and showed an approximate qualitative behaviour consistent with a phase transition at the critical point. Or course, a true phase transition is only present as $N arrow infinity$, but a computer simulation can only hope to capture an approximation of that. 

As intersting future work, the inclusion of non-square rectangular lattices, or perhaps even hexagonal models, can be an intersting way to explore this computational landscape. Also, the coincidence with the model can be further studied with the possible inclusion of second or third order neighbours.


#bibliography("bib.bib")

/*
// Workaround until balanced columns are available
// See https://github.com/typst/typst/issues/466
#place(
  bottom,
  scope: "parent",
  float: true,
  clearance: 0pt, // TODO: increase clearance for manual column balancing
  [],
)
*/