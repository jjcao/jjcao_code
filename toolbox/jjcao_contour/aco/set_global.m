% Define global variables
global num_iterations;
global num_ants;
global initial_pheromone;
global pheromone_influence;
global pheromone_persistency;
global pheromone_deposit;
global minimum_pheromone;
global descriptor_influence;
global distance_sim_gaussian;
global descriptor_sim_gaussian;

global order_preserving;
global contours;
global open_contour;
global to_one;

global feature_points;

global verbose;
global show_progress;
global show_results;

% Set default values
num_iterations = 1000;
num_ants = 1;
initial_pheromone = 1.0;
pheromone_influence = 0.3;
pheromone_persistency = 0.95;
pheromone_deposit = 0.01;
minimum_pheromone = 0.1;
descriptor_influence = 0.3;
distance_sim_gaussian = 0.1;
descriptor_sim_gaussian = 0.1;

order_preserving = 1;
contours = 1;
open_contour = 0;
to_one = 1;

feature_points = [];

verbose = 1;
show_progress = 0;
show_results = 1;
