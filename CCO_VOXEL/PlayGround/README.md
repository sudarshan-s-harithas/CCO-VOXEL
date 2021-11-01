## Playground

This folder consists of a simple python implementation of CEM logic for a single Query point and a single point obstacle. Execute the *Test_logic.py* program the obstacle is placed at the location (1,1) on a 2D plane, as expected with every CEM iteration the cost reduces and the  distribution over constraint violation tends towards the dirac function. 

Furthermore the *mmd_cost.py* provides for various functions to evaluate the MMD cost at a single point, a trajectory or a batch of trajectories, additionally to improve the execution speed the reduced set method has been used.   