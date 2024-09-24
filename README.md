# PD-Diffusion-Sonification-Project
This is an audio-visual representation of the Poisson Dirichlet Diffusion Process made in collaboration with Dr. Noah Forman from McMaster University. I'd also like to credit his other collaborators Rey Chou, Alex Forney, Chengning Li, and Gerandy Brito for prior work on the Julia scripts. The media is [linked here](https://www.youtube.com/shorts/s_6QK3bPJns). I also give credit to Matt Russo, whose [sonification tutorial](https://github.com/SYSTEMSounds/sonification-tutorials) I've used as a template for Sonification_Script.py.

## How it works
The Poisson Dirichlet Diffusion Process is a Continuous Markov Process that can be used to model situations such as population genetics. In fact, one can imagine the circles in the video as packs of asexually-reproducing bison roving about a vast plain. Each packs' or circles' movement is modelled by a random walk. The growth of the circles can described as a pack frowing larger until dying off at a random time (modelled by the Galton-Watson Branching Process). Moreover, each pack will have a member leaving at random, and starting a new one that behaves the same as it's "mother" pack. 

Note that in this project, our focus was not to model the behavior of bison herds; the above was just a useful analogy. Rather, we simulated random numbers with the Poisson Dirichlet Diffusion Process to create an animation and a musical representation of the process. Below is a brief summary of how the codebase works.

The Julia scripts: 
* ScaffSpind.jl 
  * Generates circle size, 2-dimensional coordinates of the circles, "lifespan", and "parents" of each circle 
* Skewer_NF 
  * Assigns a discrete time to a group of circle sizes, 2-dimensional coordinates, and parent status. This results in time frames that inform us where the circle is and whether if it is alive
* Create_CSV.jl 
  * Organizes the time frames into CSV files

R Script:
* SST_Animate.R 
  * Uses ggplot2 to create individual frames from each CSV file to be uploaded to a GIF-making site. 

Python Script: 
* Sonification_Script.py 
  * Creates a midi file. For each CSV file, it maps circle size to volume and one of the two coordinates to pitch.

## How to use
Here are some languages and packages that should be installed:

* For Julia:
  * Distributions
  * Random
  * DataFrames
  * CSV 
* For R: 
  * ggplot2 
  * ggforce 
* Python: 
  * numpy
  * Pandas
  * MIDI-UTIL
  * audiolazy 

The workflow will be available soon (at the moment, we are working on a Bash script to streamline this process):
<!---**Part 1**
1. In the Julia REPL, import the required packages with `using <Package_Name>`. If package isn't there, use `using Pkg; Pkg.add("Package Name")`
2. Then import ScaffSpind.jl using `include("path\to\ScaffSpind.jl")`.
3. Use `SST; Clr`. The former stores a) 
-->


