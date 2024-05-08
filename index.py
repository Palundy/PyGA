from GeneticAlgorithm import GeneticAlgorithm
from Organism import Organism
from Selection import Selection

import matplotlib.pyplot as plt
    


# Knapsack problem items
# 
items = {
    "bowling_ball": {"weight": 10, "value": 100},
    "guitar": {"weight": 1, "value": 150},
    "laptop": {"weight": 3, "value": 200},
    "tv": {"weight": 15, "value": 300},
    "painting": {"weight": 2, "value": 70},
    "books": {"weight": 1, "value": 30},
    "plants": {"weight": 1, "value": 20},
    "candles": {"weight": 1, "value": 10000},
    "pillows": {"weight": 1, "value": 15},
    "blankets": {"weight": 1, "value": 20},
    "headphones": {"weight": 1, "value": 50},
    "camera": {"weight": 2, "value": 100},
    # Add some radical items
    "diamonds": {"weight": 1, "value": 1000},
    "gold": {"weight": 1, "value": 500},
    "silver": {"weight": 1, "value": 200},
    "platinum": {"weight": 1, "value": 300},
    "copper": {"weight": 1, "value": 100},
}


def fitness_function(properties: dict):
    """
    This function calculates the fitness of an organism based on the items in the knapsack.

    ### Parameters
    - properties (`dict`): The properties of the organism.

    ### Returns
    - float: The fitness of the organism.
    """

    # Initialize the total weight and value of the knapsack
    total_weight = 0
    total_value = 0

    # Calculate the total weight and value of the knapsack
    for item, gene in properties.items():
        if gene == 1:
            total_weight += items[item]["weight"]
            total_value += items[item]["value"]

    # Penalize the fitness if the total weight exceeds the maximum weight
    max_weight = 20
    if total_weight > max_weight:
        total_value -= (total_weight - max_weight) * 10
    return total_value



# Initialize the genetic algorithm
GA = GeneticAlgorithm(fitness_function)

# Configure the chromosome of the organisms
for item in items:
    # Represent all the items in the chromosome
    GA.add_gene(item, 0, 2, 1)

# Configure the genetic algorithm
GA.recombination_method(GeneticAlgorithm.RECOMBINATION_MULTI_POINT_CROSSOVER)
GA.multi_point_num_points = 4
GA.selection_method(GeneticAlgorithm.SELECTION_TOURNAMENT)
GA.generate_population(100)
GA.selection_size = 5
GA.tournament_size = 3

mean_fitness = []
best_fitness = []
champion_organism_fitness = 0
champion_organism = None
for i in range(100):

    # Evaluate the fitness of the population
    GA.evaluate_population()

    # Print the best organism
    best_organism = GA.best_organism()

    if champion_organism_fitness == 0 or champion_organism_fitness < best_organism.fitness:
        champion_organism_fitness = best_organism.fitness
        champion_organism = best_organism
        
    # Append the values to the list
    mean_fitness.append(GA.mean_fitness())
    best_fitness.append(best_organism.fitness)

    # Generate the next population
    # this method recombines and mutates the population
    GA.generate_next_population()

# Show the stats of the best organism
print(f"Best Organism: {GA.describe_organism(champion_organism)}")

plt.plot(range(len(mean_fitness)), mean_fitness, label="Mean Fitness")
plt.plot(range(len(best_fitness)), best_fitness, label="Best Fitness")
plt.xlabel("Generation")
plt.ylabel("Fitness")
plt.legend()
plt.show()