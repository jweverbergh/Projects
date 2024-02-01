# Natural convection and mixing

This project has been done together with [Christophe Heneffe](https://github.com/ChristopheHeneffe).

The goal is to analyze the convection of a fluid in a 2D box heated from below, in these two cases : 
   - natural convection;
   - convection with a mixer.

 <br>  

The equations to solve are the following :

$$\nabla \cdot \textbf{v} = 0$$

$$\frac{D\textbf{v}}{Dt} = -\nabla P + \nu \nabla^2 \textbf{v} - \beta (T-T_{\infty}) \textbf{g} - \chi \frac{(\textbf{v} - \textbf{v}_{\textbf{s}})}{\Delta \tau}$$

$$\frac{DT}{Dt} = \alpha \nabla^2T - \chi \frac{(T-T_s)}{\Delta \tau}$$

with

- $\textbf{v}$ : velocity vector of the fluid;
- $\textbf{v}_{\textbf{s}}$ : velocity vector of the mixer;
- $T$ = temperature of the fluid;
- $T_s$ = temperature of the mixer;
- $T_{\infty}$ = temperature of the atmosphere outside the box;
- $\nu = \frac{\mu}{\rho}$ : kinematic viscosity;
- $\textbf{g}$ : gravitational accelation;
- $P$ : kinematic pressure;
- $\beta$ : fluid thermal expansion coefficient;
- $\alpha$ : thermal diffusivity;
- $\Delta \tau$ : parameter to choose, to take the mixer into account (dimension of time);
- $\chi$ : mask function to add the mixer ($\chi=1$ on the mixer, 0 otherwise).


## Run

Before compiling, make sure that you installed PETSC library. Then, modify the two first lines in the Makefile to set the path towards the PETSC library.

To compile the code, go into the Code directory, and execute the command: \
```make```

To run it, launch the executable with the parameters explained later, for example: \
```./project 0 0 0.005 10 0 0.5```

## Parameters

Parameters when we run ```./project arg1 arg2 arg3 arg4 arg5 arg6``` :
   - arg1 [int] : 1 to save results in .txt files (to plot them later), 0 otherwise
   - arg2 [int] : 1 to print parameters and some intermediate results, 0 otherwise
   - arg3 [double] : spatial step h
   - arg4 [double] : adimensional time of the end of the simulation
   - arg5 [int] : 1 to add the mixer (case 2), 0 othewise (case 1)
   - arg6 [double] : CFL number


## Simulations of the report

To obtain the simulations of the report, choose the following parameters: \
    - case 1 (no mixer) : ```./project 0 0 0.005 10000 0 0.5``` \
    - case 2 (mixer) : ```./project 0 0 0.005 10000 1 0.5``` (to be more precise, also change 0.06 to 0.08 in the definition of dt at line 456 of project.c)


For details, check out the [report](report.pdf).


## Results visualization

Remark : the simulations may take a few seconds to load on the webpage.

<p align="center">
    <img src="Images/anim_nomixer.gif" />
    <br>  
    Figure 1 : Simulation without mixer <br> Temperature | Vorticity | Norm of the velocity 
</p>


<p align="center">
    <img src="Images/anim_mixer.gif" />
    <br>  
    Figure 2 : Simulation with mixer <br> Temperature | Vorticity | Norm of the velocity
</p>






