# PLUMED Masterclass 22.11: Variationally Enhanced Sampling

## Author
Omar Valsson  
Department of Chemistry, University of North Texas   
Webpage: www.valsson.info  
July 4, 2022

## Aims

In this Masterclass, we will discuss how to run variationally enhanced sampling simulations with PLUMED. We will also understand how to analyze the results.

## Objectives

Once you have completed this Masterclass you will be able to:

- Run variationally enhanced sampling simulations biasing one and two CVs
- Assess the convergence of the simulation
- Perform reweighting from variationally enhanced sampling simulations

## Setting up PLUMED

For this masterclass you will need versions of PLUMED (with the VES module enabled) and GROMACS that are compiled using the MPI library. All the exercises were tested with PLUMED version 2.8.0 and GROMACS 2020.6. For example, to obtain these versions, you can follow the instructions at [this link](https://github.com/plumed/masterclass-2022).

The data needed to execute the exercises of this Masterclass can be found on [GitHub](https://github.com/valsson-group/masterclass-22-11).
You can clone this repository locally on your machine using the following command:

```
git clone https://github.com/valsson-group/masterclass-22-11.git
```

These files include the GROMACS and PLUMED files needed to run this tutorial, along with Jupyter notebook `Analysis.ipynb` that you can use to plot and analyze the results of the simulations. 

## Summary of Theory

Here, we will briefly summarize the theory behind variationally enhanced sampling (VES). For an full overview of VES, you should read the [original paper](https://doi.org/10.1103/PhysRevLett.113.090601), a  review in a [book chapter on VES](https://doi.org/10.1007/978-3-319-44677-6_50), or a recent [enhanced sampling review](https://doi.org/10.33011/livecoms.4.1.1583) that includes discussion about VES.

VES is based on the the following functional of the bias potential:

$$\Omega [V]  =
\frac{1}{\beta} \log
\frac
{\int d\mathbf{s} \ e^{-\beta \left[ F(\mathbf{s}) + V(\mathbf{s})\right]}}
{\int d\mathbf{s} \ e^{-\beta F(\mathbf{s})}}
+
\int d\mathbf{s} \ p_{\mathrm{tg}}(\mathbf{s}) \ V(\mathbf{s}),$$

where $\mathbf{s}$ are the CVs that we are biasing,
$p_{\mathrm{tg}}(\mathbf{s})$ is a predefined probability distribution that we will refer
to as the target distribution, and $F(\mathbf{s})$ is the free energy
surface. This functional can be shown to be convex and to have a minimum at:

$$V(\mathbf{s}) = -F(\mathbf{s})-{\frac {1}{\beta}} \log {p_{\mathrm{tg}}(\mathbf{s})}.$$

The last equation states that when we minimize the functional $\Omega [V]$,
we can obtain the free energy surface from the bias potential (and the target distribution).
We can choose the target distribution $p_{\mathrm{tg}}(\mathbf{s})$ at will and it is
the CV distribution that we obtain when minimizing $\Omega [V]$.

We put the variational principle to practice by expanding $V(\mathbf{s})$
in some basis set:

$$V(\mathbf{s}) = \sum\limits_{i} \alpha_i \ f_i(\mathbf{s}),$$

where $f_i(\mathbf{s})$ are the basis functions and the $\boldsymbol\alpha$ are the coefficients in the expansion.

We then need to find the coefficients $\boldsymbol\alpha$ that minimize $\Omega
[V]$. In principle one could use any optimization algorithm. In practice
the algorithm that has become the default choice for VES is the so-called
[averaged stochastic gradient descent algorithm](http://papers.nips.cc/paper/4900-non-strongly-convex-smooth-stochastic-approximation-with-convergence-rate-o1n.pdf).
In this algorithm the $\boldsymbol\alpha$ are evolved iteratively
according to:

$$\boldsymbol\alpha^{(n+1)} = \boldsymbol\alpha^{(n)}-\mu
 \left[
\nabla\Omega(\bar{\boldsymbol\alpha}^{(n)})+
H(\bar{\boldsymbol\alpha}^{(n)})[\boldsymbol\alpha^{(n)}-\bar{\boldsymbol\alpha}^{(n)}]
\right]$$

where $\mu$ is the step size,
$\bar{\boldsymbol\alpha}^{(n)}$ is the running average of $\boldsymbol\alpha^{(n)}$ at iteration $n$, and
$\nabla\Omega(\bar{\boldsymbol\alpha}^{(n)})$ and
$H(\bar{\boldsymbol\alpha}^{(n)})$ are the gradient and Hessian of $\Omega[V]$ evaluated at the running
average at iteration $n$, respectively.
The behavior of the coefficients will become apparent in the examples below.

As said above, we can choose the target distribution $p_{\mathrm{tg}}(\mathbf{s})$ at will.
The most simple choice would be a uniform target distribution. However, it has been found more optimal to
employ the so-called well-tempered distribution  (see [here](http://doi.org/10.1021/acs.jctc.5b00076)):

$$p_{\mathrm{tg}}(\mathbf{s}) =
\frac{[ P(\mathbf{s}) ]^{1/\gamma}}
{\int d\mathbf{s}\ [ P(\mathbf{s}) ]^{1/\gamma}}$$

where $\gamma$ is the so-called bias
factor and $P(\mathbf{s})$ is the unbiased CV distribution. Therefore the
well-tempered distribution is the unbiased distribution with enhanced fluctuations
and lowered barriers. This is the same distribution as sampled in well-tempered
metadynamics. The advantages of this distribution are that the features of the
FES (metastable states) are preserved and that the system is not forced to sample regions of high
free energy (that can represent un-physical configurations)
as it would if we had chosen the uniform target distribution.
There is a caveat though, the well-tempered $p_{\mathrm{tg}}(\mathbf{s})$ depends on 
$F(\mathbf{s})$ that is the function that we are trying to calculate.
One way to approach this problem is to calculate $p_{\mathrm{tg}}(\mathbf{s})$
self-consistently (see [here](http://doi.org/10.1021/acs.jctc.5b00076)), for instance at iteration $k$:

$$p^{(k+1)}(\mathbf{s})=\frac{e^{-(\beta/\gamma) F^{(k+1)}(\mathbf{s})}}{\int d\mathbf{s} \ e^{-(\beta/\gamma) F^{(k+1)}(\mathbf{s})}}$$

where:

$$F^{(k+1)}(\mathbf{s})=-V^{(k)}(\mathbf{s}) - \frac{1}{\beta} \log p^{(k)}(\mathbf{s})$$

Normally $p^{(0)}(\mathbf{s})$ is taken to be uniform.
Therefore the target distribution evolves in time until it becomes stationary
when the simulation has converged. It has been shown that in some cases the
convergence is faster using the well-tempered target distribution than using
the uniform $p(\mathbf{s})$  (see [here](http://doi.org/10.1021/acs.jctc.5b00076)).

## The System: Association/Dissociation of NaCl in Aqueous Solution

In this tutorial, we will consider the association/dissociation of NaCl in aqueous solution. The system consists of 1 Na atom, 1 Cl atom, and 107 water molecules for a total of 323 atoms. In an effort to speed up the simulations, we employ a rather small water box, and thus need to employ smaller cutoffs than usually used. Therefore, this simulation setup should not be used in production runs. Typically, the run should take around 15-20 minutes to run on a laptop using two MPI processes. By running the simulations on a cluster, you reduce the simulation time.

The most relevant CV for this system is the distance between the Na and Cl atoms
that is defined in PLUMED as
```plumed
# Distance between Na and Cl atoms
dist:  DISTANCE ATOMS=322,323
```

Furthermore, the NaCl association/dissociation is coupled to the collective motion
of the solvent. To measure that, we will use a CV that measures the solvation of the
Na atom. For this, we employ the coordination number of the Na atom with respect to
the oxygens of the water molecules that we define in PLUMED as
```plumed
# Solvation of Na atom
COORDINATION ...
 GROUPA=322
 GROUPB=1-321:3
 SWITCH={RATIONAL R_0=0.315 D_MAX=0.5 NN=12 MM=24}
 NLIST
 NL_CUTOFF=0.55
 NL_STRIDE=10
 LABEL=coord
... COORDINATION
```

We will also limit CV space exploration by employing an upper wall on the distance between
Na and Cl atoms that is defined in PLUMED as
```plumed
#HIDDEN
# Distance between Na and Cl atoms
dist: DISTANCE ATOMS=322,323
#ENDHIDDEN

UPPER_WALLS ...
   ARG=dist
   AT=0.6
   KAPPA=4000.0
   LABEL=uwall
... UPPER_WALLS
```

## Exercise 1: Biasing with One Collective Variable

We will start by performing a simulation where we bias the Na-Cl distance.

Every VES simulation has three key ingredients:

- Basis set
- Target distribution
- Optimization algorithm

For the basis set, we will use the recently introduced [wavelet-based basis set](https://doi.org/10.1021/acs.jctc.2c00197) that are localized basis functions that have been shown to perform better than the previously used Chebyshev or Legendre polynomials. We will employ the least asymmetric variant of these wavelets or so-called symlets (as indicated by the `TYPE=SYMLETS` keyword). We will use an order 10 of the symlets or Sym10 (as indicated by the `ORDER=10` keyword). Furthermore information about the wavelets can be found in the reference above.

We need to select the range on which the bias potential is expanded. Here we will use the range from 0.2 nm to 0.7 nm (as indicated by the `MINIMUM=0.2` and `MAXIMUM=0.7` keywords). We also need to select the number of basis functions, and here we will use 26 basis functions (as indicated by the `NUM_BF=26` keyword). The PLUMED action corresponding to this setup is given by
```plumed
#HIDDEN
# Distance between Na and Cl atoms
dist: DISTANCE ATOMS=322,323

# Solvation of Na atom
COORDINATION ...
  LABEL=coord
  GROUPA=322
  GROUPB=1-321:3
  SWITCH={RATIONAL R_0=0.315 D_MAX=0.5 NN=12 MM=24}
  NLIST
  NL_CUTOFF=0.55
  NL_STRIDE=10
... COORDINATION
#ENDHIDDEN

# Basisset for Na-Cl distance
BF_WAVELETS ...
  LABEL=bf1
  TYPE=SYMLETS
  ORDER=10
  NUM_BF=26
  MINIMUM=0.2
  MAXIMUM=0.7
  TAILS_THRESHOLD=0.01
... BF_WAVELETS
```

For the target distribution, we employ a well-tempered distribution with a bias factor
of 10.
```plumed
#HIDDEN
# Distance between Na and Cl atoms
dist: DISTANCE ATOMS=322,323

# Solvation of Na atom
COORDINATION ...
  LABEL=coord
  GROUPA=322
  GROUPB=1-321:3
  SWITCH={RATIONAL R_0=0.315 D_MAX=0.5 NN=12 MM=24}
  NLIST
  NL_CUTOFF=0.55
  NL_STRIDE=10
... COORDINATION

# Basisset for Na-Cl distance
BF_WAVELETS ...
  LABEL=bf1
  TYPE=SYMLETS
  ORDER=10
  NUM_BF=26
  MINIMUM=0.2
  MAXIMUM=0.7
  TAILS_THRESHOLD=0.01
... BF_WAVELETS
#ENDHIDDEN

# Target distribution
td: TD_WELLTEMPERED BIASFACTOR=10
```

Then we define the VES bias potential using the `VES_LINEAR_EXPANSION` action
```plumed
#HIDDEN
# Distance between Na and Cl atoms
dist: DISTANCE ATOMS=322,323

# Solvation of Na atom
COORDINATION ...
  LABEL=coord
  GROUPA=322
  GROUPB=1-321:3
  SWITCH={RATIONAL R_0=0.315 D_MAX=0.5 NN=12 MM=24}
  NLIST
  NL_CUTOFF=0.55
  NL_STRIDE=10
... COORDINATION

# Basisset for Na-Cl distance
BF_WAVELETS ...
  LABEL=bf1
  TYPE=SYMLETS
  ORDER=10
  NUM_BF=26
  MINIMUM=0.2
  MAXIMUM=0.7
  TAILS_THRESHOLD=0.01
... BF_WAVELETS

# Target distribution
td: TD_WELLTEMPERED BIASFACTOR=10
#ENDHIDDEN

VES_LINEAR_EXPANSION ...
  LABEL=ves
  ARG=dist
  BASIS_FUNCTIONS=bf1
  TEMP=300.0
  GRID_BINS=300
  OPTIMIZATION_THRESHOLD=0.000001
  TARGET_DISTRIBUTION=td
... VES_LINEAR_EXPANSION
```

Finally, we need to define the optimization algorithm. The standard is the averaged stochastic gradient descent (`OPT_AVERAGED_SGD`).
We need to define two parameters: the stride and the step size.
The stride (given by the keyword `STRIDE`) is the number of MD steps in which samples
are collected to calculate the gradient and hessian of $\Omega [V]$. The
step size (given by the keyword `STEPSIZE`) is the step by which the coefficients
are evolved at every optimization steps, given by $\mu$ in the equation above.
It has become
traditional to choose a stride of around 500-2000 MD steps. It must be noted that we
are not looking for an accurate estimation of the gradient, since for this we
would need to sample all the CV space. The step size in the
optimization has a strong connection with the height of typical barriers in
the system. The larger the barriers, the larger the step size needed such that
the bias can grow fast enough to overcome them. For this example we have
chosen a stride of 500 steps (i.e., 1 ps) and a step size of 5.0 kJ/mol.
We also need to choose how often we update the target distribution (given by the
keyword `TARGETDIST_STRIDE`) and we do this every 100 bias potential updates
(i.e., every 100 ps in the current case).

```plumed
#HIDDEN
# Distance between Na and Cl atoms
dist: DISTANCE ATOMS=322,323

# Solvation of Na atom
COORDINATION ...
  LABEL=coord
  GROUPA=322
  GROUPB=1-321:3
  SWITCH={RATIONAL R_0=0.315 D_MAX=0.5 NN=12 MM=24}
  NLIST
  NL_CUTOFF=0.55
  NL_STRIDE=10
... COORDINATION

# Basisset for Na-Cl distance
BF_WAVELETS ...
  LABEL=bf1
  TYPE=SYMLETS
  ORDER=10
  NUM_BF=26
  MINIMUM=0.2
  MAXIMUM=0.7
  TAILS_THRESHOLD=0.01
... BF_WAVELETS

# Target distribution
td: TD_WELLTEMPERED BIASFACTOR=10

VES_LINEAR_EXPANSION ...
  LABEL=ves
  ARG=dist
  BASIS_FUNCTIONS=bf1
  TEMP=300.0
  GRID_BINS=300
  OPTIMIZATION_THRESHOLD=0.000001
  TARGET_DISTRIBUTION=td
... VES_LINEAR_EXPANSION
#ENDHIDDEN

OPT_AVERAGED_SGD ...
  LABEL=opt
  BIAS=ves
  STRIDE=500
  STEPSIZE=5.0
  FES_OUTPUT=100
  BIAS_OUTPUT=500
  COEFFS_OUTPUT=10
  TARGETDIST_STRIDE=100
  TARGETDIST_OUTPUT=500
... OPT_AVERAGED_SGD
```
The other parameters are related to the outputting frequency of various output files.

Finally, we need to define the `PRINT` action to output all the variables
```plumed
#HIDDEN
# Distance between Na and Cl atoms
dist: DISTANCE ATOMS=322,323

# Solvation of Na atom
COORDINATION ...
  LABEL=coord
  GROUPA=322
  GROUPB=1-321:3
  SWITCH={RATIONAL R_0=0.315 D_MAX=0.5 NN=12 MM=24}
  NLIST
  NL_CUTOFF=0.55
  NL_STRIDE=10
... COORDINATION

# Basisset for Na-Cl distance
BF_WAVELETS ...
  LABEL=bf1
  TYPE=SYMLETS
  ORDER=10
  NUM_BF=26
  MINIMUM=0.2
  MAXIMUM=0.7
  TAILS_THRESHOLD=0.01
... BF_WAVELETS

# Target distribution
td: TD_WELLTEMPERED BIASFACTOR=10

VES_LINEAR_EXPANSION ...
  LABEL=ves
  ARG=dist
  BASIS_FUNCTIONS=bf1
  TEMP=300.0
  GRID_BINS=300
  OPTIMIZATION_THRESHOLD=0.000001
  TARGET_DISTRIBUTION=td
... VES_LINEAR_EXPANSION

OPT_AVERAGED_SGD ...
  LABEL=opt
  BIAS=ves
  STRIDE=500
  STEPSIZE=5.0
  FES_OUTPUT=100
  BIAS_OUTPUT=1000
  COEFFS_OUTPUT=10
  TARGETDIST_STRIDE=100
  TARGETDIST_OUTPUT=500
... OPT_AVERAGED_SGD

UPPER_WALLS ...
   ARG=dist
   AT=0.6
   KAPPA=4000.0
   LABEL=uwall
...
#ENDHIDDEN

PRINT ARG=dist,coord,ves.*,uwall.* FILE=colvar.data STRIDE=125
```

The full PLUMED input file can be found in the Exercise-1 folder. There you also find all
the relevant GROMACS input file. First, you need to run the
`generate-tpr-file.sh` script that generates the GROMACS TPR file using the parameters defined in `MD-NPT.mdp` mdp file 
and the initial geometry defined using the `StartingGeometry` variable. You can then run the simulation using
the `run.sh` script
```
./generate-tpr-file.sh
./run.sh
```
The run might take around 15-20 minutes to run using two MPI processes. You can adjust the number of MPI processes used for the simulation using the `NumProcs` variable in the `run.sh` script.

At the end of simulation, you will get several files:
- `colvar.data`: Colvar file
- `coeffs.data` : Values of the coefficients $\boldsymbol\alpha$ and $\bar{\boldsymbol\alpha}$ at different iterations.
- `bias.<bias-label>.iter-<iteration-number>.data` : Bias potential at iteration <iteration-number>.
- `fes.<bias-name>.iter-<iteration-number>.data` : FES at iteration <iteration-number>.
- `targetdistribution.<bias-name>.iter-<iteration-number>.data` : Target distribution at iteration <iteration-number>.

To assess the simulation and its convergence, you should first look at the time evolution of the
biased CV and check that it is diffusive in CV space. Second, you should look at how the free energy surfaces behave as
a function of time by looking at the `fes.<bias-name>.iter-<iteration-number>.data` files at different number of iterations (the minimum of the FES is always aligned to zero
to facilitate comparison). To do this, you need to use your favorite way to plot datafiles (e.g., Matplotlib or Gnuplot). For instance, you can use the `Analysis.ipynb` Jupyter notebook provided. 

You can also visualize the trajectory by opening it with VMD by using the command
```
vmd NaCl_StartingStructure-1.gro NaCl_VES_NPT-300K.pbc-whole.xtc
```
The `NaCl_VES_NPT-300K.pbc-whole.xtc` is the trajectory file with the periodic boundary conditions made whole. This is done within the `run.sh` script.

The `coeffs.data` file includes the values of coefficients $\boldsymbol\alpha$ and $\bar{\boldsymbol\alpha}$ at different iterations. To extract the time evolution of
a given coefficient, you can use the `ExtractCoeff.sh` script, for example to extract
coefficient 3:
```
./ExtractCoeff.sh 3 > coeff.3.data
```
This will create a file with the first column the iteration number, the second
column the averaged coefficient $\bar{\boldsymbol\alpha}$, and the third column
the instantaneous coefficients $\boldsymbol\alpha$. You should create files for
different coefficient and visualize both the second and third columns to understand
how the coefficients converge.

## Exercise 2: Reweighting a VES Simulation
Apart from estimating the FES directly from the bias potential, you can also estimate
the FES through reweighting by histogramming where each configuration is weighted by the
bias acting on it, $e^{\beta V(\mathbf{s})}$. The VES bias acting at each time step
is given by the `ves.bias` variable in the `colvar.dat` file.

When doing performing reweighting, it is better to ignore the initial part of
the simulation where the bias potential is changing more rapidly. You can use the
`trim-colvar-file.py` python script in the Exercise-2 folder to do this
```
./trim-colvar-file.py --colvar-file ../Exercise-1/colvar.data --output-file colvar_reweight.data --time-min 400
```
where here we ignore the first 400 ps of the `colvar.data` file from the Exercise 1 and create
a new file called `colvar_reweight.data`.

We can then perform the reweighting for the distance using the following PLUMED input
(`plumed_reweight.dat` in the `Exercise-2` folder)
```plumed
dist:   READ FILE=colvar_reweight.data IGNORE_TIME VALUES=dist
coord:  READ FILE=colvar_reweight.data IGNORE_TIME VALUES=coord
ves:    READ FILE=colvar_reweight.data IGNORE_TIME VALUES=ves.bias

weights: REWEIGHT_BIAS TEMP=300 ARG=ves.bias

HISTOGRAM ...
  ARG=dist
  GRID_MIN=0.2
  GRID_MAX=0.7
  GRID_BIN=60
  KERNEL=DISCRETE
  LOGWEIGHTS=weights
  LABEL=hg_dist
... HISTOGRAM

fes_dist: CONVERT_TO_FES GRID=hg_dist TEMP=300 MINTOZERO
DUMPGRID GRID=fes_dist FILE=fes-reweight.dist.data FMT=%24.16e
```

You can run this input by using the PLUMED driver
```
plumed driver --plumed plumed_reweight.dat --noatoms
```

You should compare the resulting FES (the `fes-reweight.dist.data` file)
to the results obtained directly from the bias potential in Exercise 1.

We can also obtain the FES for CVs that are not biased in the VES simulation.
For example, we can obtain the two-dimensional FES for the distance and the
solvation CV given by the coordination number CV. For this, you will need to
add the following to the plumed_reweight.dat file and repeat the PLUMED driver run

```plumed
#HIDDEN
dist:   READ FILE=colvar_reweight.data IGNORE_TIME VALUES=dist
coord:  READ FILE=colvar_reweight.data IGNORE_TIME VALUES=coord
ves:    READ FILE=colvar_reweight.data IGNORE_TIME VALUES=ves.bias

weights: REWEIGHT_BIAS TEMP=300 ARG=ves.bias
#ENDHIDDEN

HISTOGRAM ...
  ARG=dist,coord
  GRID_MIN=0.2,2.5
  GRID_MAX=0.7,7.5
  GRID_BIN=200,200
  BANDWIDTH=0.004,0.04
  LOGWEIGHTS=weights
  LABEL=hg_dist_coord
... HISTOGRAM
fes_dist_coord: CONVERT_TO_FES GRID=hg_dist_coord TEMP=300 MINTOZERO
DUMPGRID GRID=fes_dist_coord FILE=fes-reweight.dist-coord.data FMT=%24.16e
```
Note that here we use kernel density estimation with Gaussian kernels to
obtain smoother results.

You can also try to obtain the one-dimensional FES for the solvation CV by adjusting
the input above.

This will generate a two-dimensional FES that you can visualize.

## Exercise 3: Running Another Independent Simulation

To check the results, it is a good practice to run another independent simulation
using different initial conditions. You can achieve this here by changing the initial
geometry in the `generate-tpr-file.sh` script
```
StartingGeometry=NaCl_StartingStructure-2.gro
```
and regenerating the GROMACS tpr file. Do this and rerun the simulation,
check the convergence, and perform reweighting
as before. Make sure that you do this in a new clean folder that is separate from the run
in Exercise 1.

## Exercise 4: Biasing with Two Collective Variables

We will now bias also the solvation CV. To achieve this, we need first to setup
a separate basis set for the solvation CV, where again we use the symlets but
with a different range from 2.0 to 8.0
```plumed
#HIDDEN
# Distance between Na and Cl atoms
dist: DISTANCE ATOMS=322,323

# Solvation of Na atom
COORDINATION ...
  LABEL=coord
  GROUPA=322
  GROUPB=1-321:3
  SWITCH={RATIONAL R_0=0.315 D_MAX=0.5 NN=12 MM=24}
  NLIST
  NL_CUTOFF=0.55
  NL_STRIDE=10
... COORDINATION

# Basisset for Na-Cl distance
BF_WAVELETS ...
  LABEL=bf1
  TYPE=SYMLETS
  ORDER=10
  NUM_BF=26
  MINIMUM=0.2
  MAXIMUM=0.7
  TAILS_THRESHOLD=0.01
... BF_WAVELETS
#ENDHIDDEN

# Basisset for coordination number
BF_WAVELETS ...
  LABEL=bf2
  TYPE=SYMLETS
  ORDER=10
  NUM_BF=22
  MINIMUM=2.5
  MAXIMUM=7.5
  TAILS_THRESHOLD=0.01
... BF_WAVELETS
```

We also need to change the relevant keywords in the `VES_LINEAR_EXPANSION` action,
namely the `ARG`, `BASIS_FUNCTIONS`, and `GRID_BINS` keywords
```plumed
#HIDDEN
# Distance between Na and Cl atoms
dist: DISTANCE ATOMS=322,323

# Solvation of Na atom
COORDINATION ...
  LABEL=coord
  GROUPA=322
  GROUPB=1-321:3
  SWITCH={RATIONAL R_0=0.315 D_MAX=0.5 NN=12 MM=24}
  NLIST
  NL_CUTOFF=0.55
  NL_STRIDE=10
... COORDINATION

# Basisset for Na-Cl distance
BF_WAVELETS ...
  LABEL=bf1
  TYPE=SYMLETS
  ORDER=10
  NUM_BF=26
  MINIMUM=0.2
  MAXIMUM=0.7
  TAILS_THRESHOLD=0.01
... BF_WAVELETS

# Basisset for coordination number
BF_WAVELETS ...
  LABEL=bf2
  TYPE=SYMLETS
  ORDER=10
  NUM_BF=26
  MINIMUM=2.5
  MAXIMUM=7.5
  TAILS_THRESHOLD=0.01
... BF_WAVELETS


# Target distribution
td: TD_WELLTEMPERED BIASFACTOR=10
#ENDHIDDEN


VES_LINEAR_EXPANSION ...
  LABEL=ves
  ARG=dist,coord
  BASIS_FUNCTIONS=bf1,bf2
  TEMP=300.0
  GRID_BINS=300,300
  OPTIMIZATION_THRESHOLD=0.000001
  TARGET_DISTRIBUTION=td
  PROJ_ARG1=dist
  PROJ_ARG2=coord
... VES_LINEAR_EXPANSION
```
Additionally, we have turned on the calculation of the one-dimensional projection
of the FES on the two CVs (we also need to set the keyword
`FES_PROJ_OUTPUT=100` in `OPT_AVERAGED_SGD`). This is useful to assess the
convergence as this is
easier in one-dimension. We can also compare the projection to the results from
Exercise 1. These changes should be sufficient to do the simulations using two CVs.
You can see full input file in the `Exercise-4` folder.

Once you have performed this simulation, you should also try to reweight from
this simulations. Furthermore, if you have time, you should also try to perform
another independent simulation.

## Optional Exercises

The following three exercises are optional, but they will show you how different
parameters effect the results. You should base these exercises on the files from
the Exercise-1 folder and make the necessary changes. Make sure that you run these
simulations in separate folders and start from clean files from the Exercise-1 folder.

### Optional Exercise 5: Play with the Optimization Parameters

The main parameter in the optimization algorithm is the step size and
it is not always easy to choose this parameter. Luckily, the algorithm
is quite robust and will work for different step sizes.

Run different simulations using step sizes 0.5 and 50.0
and try to rationalize the behavior. Normally, when the step size is too large, the
system gets stuck in CV space and coefficients oscillate wildly. When the step size is
too small, the algorithm runs out of "steam" too fast and the simulation converges slowly.
These two extreme cases should be avoided.

### Optional Exercise 6: Uniform Target Distribution
Perform a simulation using an uniform target distribution and see how this changes
the results.

In this case, you need to change the target distribution to
```plumed
#HIDDEN
# Distance between Na and Cl atoms
dist: DISTANCE ATOMS=322,323

# Solvation of Na atom
COORDINATION ...
  LABEL=coord
  GROUPA=322
  GROUPB=1-321:3
  SWITCH={RATIONAL R_0=0.315 D_MAX=0.5 NN=12 MM=24}
  NLIST
  NL_CUTOFF=0.55
  NL_STRIDE=10
... COORDINATION

# Basisset for Na-Cl distance
BF_WAVELETS ...
  LABEL=bf1
  TYPE=SYMLETS
  ORDER=10
  NUM_BF=26
  MINIMUM=0.2
  MAXIMUM=0.7
  TAILS_THRESHOLD=0.01
... BF_WAVELETS
#ENDHIDDEN

# Target distribution
td: TD_UNIFORM
```
and remove the `TARGETDIST_STRIDE` and `TARGETDIST_OUTPUT` keywords from the \ref OPT_AVERAGED_SGD
action.

### Legendre polynomials basis function 
Perform a simulation using Legendre polynomials (`BF_LEGENDRE`) basis functions instead of the
wavelets and see how this will affect the result. As the Legendre polynomials are delocalized
basis functions, this should lead to more fluctuations in the bias potential as has been observed
in the [paper introducing the wavelets](https://doi.org/10.1021/acs.jctc.2c00197).

In this case, you need to change the basis set action in the PLUMED input to
```plumed
#HIDDEN
# Distance between Na and Cl atoms
dist: DISTANCE ATOMS=322,323

# Solvation of Na atom
COORDINATION ...
  LABEL=coord
  GROUPA=322
  GROUPB=1-321:3
  SWITCH={RATIONAL R_0=0.315 D_MAX=0.5 NN=12 MM=24}
  NLIST
  NL_CUTOFF=0.55
  NL_STRIDE=10
... COORDINATION
#ENDHIDDEN 
# Basisset for Na-Cl distance
BF_LEGENDRE ...
  LABEL=bf1
  ORDER=20
  MINIMUM=0.2
  MAXIMUM=0.7
... BF_LEGENDRE
```

## Additional comments

This should cover the basics of running VES simulations in PLUMED, but the following
comments might be of interest to some.

### Restarting
VES simulations can be easily restarted. The code will automatically output
all the file needed at the end of the simulation. To restart, you just need to
reset the simulation with your MD code in the traditional way and add a `RESET`
keyword to the top of the PLUMED input (for some codes like GROMACS, a restart is
automatically detected by PLUMED and thus this keyword is not needed).

### Multiple Walkers

VES simulations supports the usage of multiple walkers where different copies of the system share the same bias potential (i.e. coefficients) and cooperatively sample the averages needed for the gradient and Hessian. This can significantly help with convergence in difficult cases. It is of course best to start the different copies from different positions in CV space. To activate this option you just need to add the `MULTIPLE_WALKERS` keyword to the `OPT_AVERAGED_SGD` action. Note that this is only supported if the MD code support running multiple replicas connected via MPI (e.g., GROMACS or LAMMPS).


