import sys, os
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from pyga_package.pyga import PyGA


# ------------------------------------------------
#           This script demonstrates
#       how to use PyGA to find an optimal
#       solution for the knapsack problem
# ------------------------------------------------


# ------------------------------------------------
# In this example the single-point crossover method is used.
# This method selects a random crossover point and swaps the genes
# of the parents to create two offspring.
# 
# The tournament selection method is used to select the parents for crossover.
# This method selects a subset of the population (tournament) and chooses the best
# organism from that subset. The tournament size is set to 5 in this example.
# The recombination process is repeated until the given selection size is reached,
# which is set to 5 in this example.
# 
# A small mutation rate of 0.01 (1%) is used to introduce genetic diversity.
# ------------------------------------------------



# Define the parameters for the genetic algorithm
max_weight = 35
population_size = 500
generations = 1000
mutation_probability = 0.01
recombination_method = PyGA.RECOMBINATION_SINGLE_CROSSOVER
selection_method = PyGA.SELECTION_TOURNAMENT
tournament_size = 5
selection_size = 5
dynamic_mutation = True

# Initialze the array to track the fitness score
fitness_scores = []


# Define the items that can be placed in the knapsack
items = {
    "book": {"weight": 1, "value": 600},
    "binoculars": {"weight": 2, "value": 500},
    "map": {"weight": 2, "value": 400},
    "compass": {"weight": 1, "value": 250},
    "water": {"weight": 3, "value": 350},
    "sandwich": {"weight": 1, "value": 160},
    "glucose": {"weight": 1, "value": 60},
    "tin": {"weight": 2, "value": 45},
    "banana": {"weight": 3, "value": 60},
    "watch": {"weight": 1, "value": 10},
    "umbrella": {"weight": 1, "value": 40},
    "sunglasses": {"weight": 1, "value": 20},
    "notebook": {"weight": 1, "value": 30},
    "pencil": {"weight": 1, "value": 10},
    "gold_bar": {"weight": 8, "value": 1000},
}


# Define the fitness function
def fitness_function(properties: dict)-> int:
    """
    This method calculates the fitness of an organism based on the total value of the items in the knapsack.
    If the total weight of the items exceeds the maximum weight of the knapsack, the fitness is set to 1.

    ### Parameters
    - properties (`dict`): The properties of the organism.

    ### Returns
    - `int`: The fitness value.
    """

    total_weight = 0
    total_value = 0
    global max_weight

    # Calculate the total weight and value of the selected items
    for item, gene_value in properties.items():
        gene_value = round(gene_value)
        if gene_value > 0:
            total_weight += items[item]["weight"] * gene_value
            total_value += items[item]["value"] * gene_value

    # If the total weight exceeds the maximum weight, set the fitness to 1
    if total_weight > max_weight:
        return 1
    return total_value # The fitness is the total value of the items



# Initialize the genetic algorithm
ga = PyGA(fitness_function) # Pass the fitness function
ga.recombination_method(recombination_method)
ga.selection_method(selection_method)
ga.selection_size = selection_size
ga.tournament_size = tournament_size
ga.mutation_probability = mutation_probability
ga.dynamic_mutation = dynamic_mutation


# Add the items as genes
for item in items:

    # The items are represented as single bits
    # within the chromosome of an organism.
    # The value of the bit indicates whether the item is selected or not.
    # 1 -> Selected, 0 -> Not selected

    item_name = item
    item_weight = items[item]["weight"]
    ga.add_gene(
        property = item_name, # The name of the item (for easy identification later)
        bottom_value = 0, # The minimum value that the gene can take
        top_value = 5, # The maximum value that the gene can take (max. amount of items)
        bitsize = 3, # The number of bits needed to represent the gene..
        # (2^3 = 8, which is more than the maximum amount of items)
    )


# Generate the genesis population
ga.generate_population(population_size)
for generation in range(generations):

    # Evaluate the fitness of the population
    ga.evaluate_population()

    if generation % 100 == 0:
        champion_organism = ga.champion_organism()
        print(f"Generation {generation}")
        print(f"\tHighest Fitness: {champion_organism.fitness}")
        print()

    # Generate the next generation
    ga.generate_next_population()

# Print the information of the best organism
champion_organism = ga.champion_organism()
champion_organism_description = ga.describe_organism(champion_organism, True)
print("Champion Organism:")
print(champion_organism_description)

# Retrieve the list with highest fitness scores
fitness_scores = ga.highest_fitnesses

# Retrieve the indices of where the dynamic mutations took place
mutation_indices = ga.dynamic_mutation_indices
mutation_fitnesses = np.array(fitness_scores)[mutation_indices]

# Show the fitness scores
plt.plot(fitness_scores)
plt.scatter(mutation_indices, mutation_fitnesses, color='red')
plt.xlabel("Generation")
plt.ylabel("Fitness Score")
plt.grid()
plt.show()