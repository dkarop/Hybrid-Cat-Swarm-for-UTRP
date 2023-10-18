# Hybrid-Cat-Swarm-for-UTRP
A matlab implemenation of Hybrid Cat Swarm Optimization Algorithm ment to solve UTRP using Mandl's Data

UTRP is a type of optimization problem that aims to find the most efficient routes of a fleet of vehicles (such as buses or trains) for their movement in the urban fabric, 
while taking into account factors such as traffic congestion, passenger demand, vehicle availability as well as the cost of moving them. 
This is an NP-hard (Nondeterministic Polynomial) problem, meaning that it is computationally expensive to find an exact solution for large path networks. 
This is also the reason why meta-heuristic algorithms are used to find approximate solutions. In this paper, a hybrid cat swarm optimization algorithm is developed to solve the urban transport routing problem.
It is worth emphasizing that Cat Swarm Optimization is not a traditional optimization algorithm, nor the most efficient one, but it has been applied to routing problems with good results and is relatively easy to implement. 
The results are compared against the experimental data of a Swiss urban bus network, which is the most widespread reference data set in the literature. 
Comparing the results with corresponding CSO-based publications shows that the performance of the proposed algorithm is superior to existing CSO techniques. 
In particular, the proposed one makes use of CDC /SRD, in a static and dynamic way, local optima avoidance techniques are implemented using reversing of paths and node addition, 
while an attempt has been made to hybridize using a Hill Climbing algorithm. 
