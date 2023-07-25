# ASSE23
## Aerospace Structures Course Project
- A stiffened panel is to be analysed regarding its structural strength and stability by means of an existing FE model. The panel is a skin segment of an aircraft wing between two spars and two ribs, which is stiffened by T-shaped stringers, as illustrated in Figure 1.
- All detailed dimensions, properties and boundary conditions can be obtained from within the FE model.
- Three representative load cases are considered. The internal section loads applied on the model boundaries are modelled as a combination of a longitudinal, transverse and shear load fraction, as seen in Figure 2. These fractions are scaled in each load case.
## How to use the scripts? 
In this repo, there exist two python files, namely Task1 and Task2. Each script is used for differnet task although they're nearly the same. 
### The process is as follows
1. Hypermesh Model
2. Getting info about 2D Stress Element in mid plane for each load case in:
   1. XX Stresses in X direction
   2. YY Stresses in Y direction
   3. XY Shear Stresses
   4. Von-Misses Stresses
3. All data are in .txt format and they're then put in their respective directory
4. The output will be a text file called LC{x}_Final, where x is the load case you're working on
