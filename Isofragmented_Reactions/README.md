The file runMethod.py contains the code needed to run the full model on the users device. 

The MILP solver can be run using the NEOS solver which will require the user to input their email on line 8 in the file. Alternatively the user can use their own solver such as Gurobi or CPLEX. In order to do this the user will need to uncomment the appropriate part of the findRxns function located around line 130. 
