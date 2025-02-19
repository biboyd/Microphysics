Takes in two of pressure,temperature,density then solves for the third for different amount of the A=23 Urca pair.

Currently takes in a C12, O16 and a sum of Na23 and Ne23 mass fractions.

Output: CSV file with following info

    * `rho`:  density in g/cm^3
    * `T`:  temperature in K
    * `p`: pressure in cgs units
    * `C12`: mass frac of carbon 12 
    * `O16`: mass frac of oxygen 16 
    * `Ne23`: mass frac of neon 23 
    * `Na23`: mass frac of sodium 23
    * `Ye`: Electron fraction 
    * `S`: entropy  in cgs
    * `e`: Internal energy in cgs
    * `eta`: the degeneracy parameter. (electron chemical potential - rest mass) / (k_b T)   where k_b is the boltzmann constant
