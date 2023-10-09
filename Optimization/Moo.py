'''
MIT License

Copyright (c) 2023 Sophia Apostolidou

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

# PARALLELIZATION USING THREADPOOLEXCECUTOR

from TEST_R_ALLLAYERS import *
from TEST_L_ALLLAYERS import *
from TEST_A_ALLLAYERS import f3

from pymoo.optimize import minimize
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.algorithms.soo.nonconvex.ga import GA
from pymoo.core.sampling import Sampling
from pymoo.core.callback import Callback
from pymoo.core.problem import Problem
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from pymoo.visualization.scatter import Scatter

import time
from datetime import timedelta

import sys

# import necessary QGIS modules
from qgis.core import *

from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor
from pymoo.core.evaluator import Evaluator
from pymoo.operators.crossover.sbx import SimulatedBinaryCrossover as sbx
from pymoo.operators.mutation.pm import PolynomialMutation as pm

from scipy.optimize import linear_sum_assignment


# Supply path to where is qgis installed
qgs = QgsApplication([], False)
qgs.initQgis()

# Set the project path
path = 'D:\My Documents\Egna_projekt\Textsättning\Testdata\London\London20211006_Clip.qgz'
project = QgsProject.instance()
project.read(path)


# What to print after each generation
class MyCallback(Callback):
    def __init__(self, problem, group_name, t_layer, g_layers, o_layers, b_layers):
        super().__init__()
        self.problem = problem
        self.group_name = group_name
        self.t_layer = t_layer
        self.g_layers = g_layers
        self.o_layers = o_layers
        self.b_layers = b_layers
        self.start_time = time.time()

    def notify(self, algorithm):
        elapsed_time = (time.time() - self.start_time) / 60
        print(f'Elapsed time after generation {algorithm.n_gen}: {elapsed_time:.2f} minutes')
        
        # The current best solution can be accessed with algorithm.pop.get("X")[0]
        x = algorithm.opt.get("X")[0]

        x_coords_new = x[:len(self.problem.x_coords)]
        y_coords_new = x[len(self.problem.x_coords):]
        
        # Call the f1, f2t and f2i functions with the current coordinates
        f1_value = f1(self.group_name, x_coords_new, y_coords_new)
        f2t_value = f2t(x_coords_new, y_coords_new, self.group_name, self.t_layer)
        f3_value = f3(x_coords_new, y_coords_new, self.g_layers, self.o_layers, self.b_layers)
        
        # Print the f1, f2t, and f3 values
        print("f1: ", f1_value)
        print("f2t: ", f2t_value)
        print("f3: ", f3_value)
        '''
        # Plot the positions of the features
        plt.scatter(self.problem.x_coords, self.problem.y_coords, color='blue', label='Original')
        plt.scatter(x_coords_new, y_coords_new, color='red', label='Optimized')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title('Comparison of feature positions at generation {}'.format(algorithm.n_gen))
        plt.legend()
        plt.show()
        '''
        # Reset the start time for the next generation
        self.start_time = time.time()


# Initialize the population of first gen with the exact same initial positions of the icons
# as I already have a good initial solution and I want to start the optimization from this point.

class MySampling(Sampling):
    def __init__(self, x_coords, y_coords):
        super().__init__()
        self.x_coords = x_coords
        self.y_coords = y_coords

    def _do(self, problem, n_samples, **kwargs):
        X = np.empty((n_samples, problem.n_var))
        for i in range(n_samples):
            if i == 0:
                X[i] = np.concatenate([np.array(self.x_coords), np.array(self.y_coords)])
            else:
                X[i] = np.concatenate([np.random.uniform(low=problem.xl[:len(self.x_coords)], high=problem.xu[:len(self.x_coords)]),
                                       np.random.uniform(low=problem.xl[len(self.x_coords):], high=problem.xu[len(self.x_coords):])])
        return X
    
class MyProblem(Problem):
    def __init__(self, x_coords, y_coords, group_name, t_layers, g_layers, o_layers, b_layers, 
                    f1_value_ul, f2t_value_ul, f3_value_ul):
        
        self.x_coords = x_coords
        self.y_coords = y_coords
        self.group_name = group_name
        self.t_layers = t_layers
        self.g_layers = g_layers
        self.o_layers = o_layers
        self.b_layers = b_layers
        
        # Upper limits for the functions
        self.f1_value_ul = f1_value_ul        
        self.f2t_value_ul = f2t_value_ul
        self.f3_value_ul = f3_value_ul
        
        super().__init__(n_var=2*len(x_coords),
                         n_obj=3, 
                         n_constr = 3, 
                         xl=np.concatenate([np.array(x_coords) - 40, np.array(y_coords) - 40]),
                         xu=np.concatenate([np.array(x_coords) + 40, np.array(y_coords) + 40]))
        
    def _evaluate(self, X, out, *args, **kwargs):
        # X is a 2D array, each row is a different individual
        F = np.full([len(X), 3], np.inf) # Initialize the objectives to a high value
        G = np.full([len(X), 3], np.inf) # Initialize the constraints


        for i, x in enumerate(X):
            # Split the decision variable into x and y parts
            x_coords_new = x[:len(self.x_coords)]
            y_coords_new = x[len(self.x_coords):]

            f1_val = f1(self.group_name, x_coords_new, y_coords_new)
            f2t_val = f2t(x_coords_new, y_coords_new, self.group_name, self.t_layers)
            f3_val = f3(x_coords_new, y_coords_new, self.g_layers, self.o_layers, self.b_layers)

            # Call the f1, f2t and f2i functions with the new coordinates
            F[i, 0] = f1_val
            F[i, 1] = f2t_val
            F[i, 2] = f3_val
            
            
            # Evaluate the constraints for functions
            # The algorithm will try to minimize the functions' values, 
            # but if it can't, it will ensure that the values don't exceed the initial values by more than 0.05.
            
            G[i, 0] = f1_val - self.f1_value_ul 
            G[i, 1] = f2t_val - self.f2t_value_ul    
            G[i, 2] = f3_val - self.f3_value_ul  
            
        out["F"] = F
        out["G"] = G

# Variables initialization
group_name = '60center'

t_layers = [project.mapLayersByName('clip_LMF_Landmark_Building_T')[0], project.mapLayersByName('clip_L_ LMF_Road_Names_T')[0]]
g_layers = getLayer(group_name)
o_layers = getLayer('Icons')
b_layers = getLayer('Background areas')[4:-1]

coords = getpositions(getLayer(group_name))
x_coords = [coord[0] for coord in coords]
y_coords = [coord[1] for coord in coords]


# Allow some error around the initial values
error_margin = 0.05

# Assuming f1_value_initial is the initial value of f1
f1_value_initial = f1(group_name, x_coords, y_coords)
f1_value_ul = f1_value_initial + error_margin

# For other functions
f2t_value_initial = f2t(x_coords, y_coords, group_name, t_layers)
f2t_value_ul = f2t_value_initial + error_margin


f3_value_initial = f3(x_coords, y_coords, g_layers, o_layers, b_layers)
f3_value_ul = f3_value_initial + error_margin


problem = MyProblem(x_coords, y_coords, group_name, t_layers, g_layers, o_layers, b_layers, 
                    f1_value_ul, f2t_value_ul, f3_value_ul)

# Initialize the sampling method
sampling = MySampling(x_coords, y_coords)

Evaluator().executor = ThreadPoolExecutor(max_workers=4)

# Initialize the NSGA2 algorithm
algorithm = NSGA2(pop_size=30, sampling=sampling, eliminate_duplicates=True) #crossover=sbx(prob=0.9, eta=15), mutation=pm(eta=20, prob=1.0)

#callback = MyCallback()
callback=MyCallback(problem, group_name, t_layers, g_layers, o_layers, b_layers)

# Change the behavior to print the entire array
np.set_printoptions(threshold=np.inf)

# Record the start time
start_time = time.time()


# Specify the path where the output file will be stored
#file_path = r"C:\Users\sophi\Desktop\SOPHIA\Master\4th semester\Opt Results\AllResults.txt"

# Create (or overwrite) the FinalResults.txt file at the specified path and set it as the default output.
#original_stdout = sys.stdout  # Save the original standard output
#sys.stdout = open(file_path, 'w')


# Run the optimization
res = minimize(problem,
               algorithm,
               ("n_gen", 61),
               callback=callback,
               verbose=True,
               seed = 1) 



# Calculate the elapsed time
elapsed_time = time.time() - start_time
f_time = str(timedelta(seconds=elapsed_time)) 


#sys.stdout.close()  # Close the file
#sys.stdout = original_stdout  # Restore the standard output to the terminal

#print(f'All results saved to {file_path}')

# Print the elapsed time
print('Elapsed time: ', f_time)

# Print the resulting decision variables and objective value
print("Best solution found for functions: \nf1 = %s\nf2t = %s\nf3 = %s" % (res.F[0][0], res.F[0][1], res.F[0][2]))

# Print all Pareto optimal solutions
print('Optimized coords:',res.X, '\n')

# Save the coordinates of the best solution in a txt file
np.savetxt(r'C:\Users\sophi\Desktop\SOPHIA\Master\4th semester\Opt Results\Opt Coords.txt', res.X, fmt='%.8f')
print('Optimized coords saved', '\n')

# Save the Pareto front data in a txt file
np.savetxt(r'C:\Users\sophi\Desktop\SOPHIA\Master\4th semester\Opt Results\ParetoFront.txt', res.F, fmt='%.8f')
print('Pareto front data saved')


#___________________________________________________CHECK BOUNDS_________________________________________________________

# Check if the new coordinates are inside the bounds
all_within_bounds = True


for solution in res.X:
    for i, value in enumerate(solution):
        if not (problem.xl[i] <= value <= problem.xu[i]):
            all_within_bounds = False
            print(f"Solution {solution} is out of bounds for variable {i}. Value: {value}, Bounds: ({problem.xl[i]}, {problem.xu[i]})")
            break
    if not all_within_bounds:
        break

if all_within_bounds:
    print("All solutions are within the bounds!")
else:
    print("Some solutions are out of bounds!")





#___________________________________________________PARETO FRONT_________________________________________________________

#3D
# get the Pareto front (assuming res is your Result object)
pareto_F = res.F

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')  # Create a 3D subplotf

# Scatter the data points
ax.scatter(pareto_F[:, 0], pareto_F[:, 1], pareto_F[:, 2], color="blue", s=30)

ax.set_xlabel('f1: Disturbance')
ax.set_ylabel('f2: Legibility (text labels)')
ax.set_zlabel('f3: Association')

# Set the x, y, and z limits
ax.set_xlim([0, 0.5])
ax.set_ylim([0, 0.5])
ax.set_zlim([0, 0.5])

plt.show()

#2D
# Plot 3 different Pareto front charts

# Plot for f1-f2t
plt.figure()
plt.scatter(pareto_F[:, 0], pareto_F[:, 1], color="blue", s=30)
plt.xlabel('f1: Disturbance')
plt.ylabel('f2: Legibility (text labels)')
plt.xlim([0, 0.5])
plt.ylim([0, 0.5])
plt.title('Pareto Front: f1 vs f2t')
plt.show()

# Plot for f1-f3
plt.figure()
plt.scatter(pareto_F[:, 0], pareto_F[:, 2], color="blue", s=30)
plt.xlabel('f1: Disturbance')
plt.ylabel('f3: Association')
plt.xlim([0, 0.5])
plt.ylim([0, 0.5])
plt.title('Pareto Front: f1 vs f3')
plt.show()

# Plot for f2t-f3
plt.figure()
plt.scatter(pareto_F[:, 1], pareto_F[:, 2], color="blue", s=30)
plt.xlabel('f2: Legibility (text labels)')
plt.ylabel('f3: Association')
plt.xlim([0, 0.5])
plt.ylim([0, 0.5])
plt.title('Pareto Front: f2t vs f3')
plt.show()
