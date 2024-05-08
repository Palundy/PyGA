from Organism import Organism

import random


class Recombination:

    # --------------------------------------------------
    # Methods
    # --------------------------------------------------
    def recombination_single_crossover(parent1: Organism, parent2: Organism)-> list:
        """
        This method performs a single-point crossover recombination between two parent organisms.
        In single-point crossover, a random crossover point is selected, and the genes of the parents
        are exchanged to create two new offspring.

        ### Parameters
        - parent1 (`Organism`): The first parent organism.
        - parent2 (`Organism`): The second parent organism.

        ### Returns
        - `list`: A list of two offspring organisms.
        """

        # Select a random crossover point
        crossover_point = random.randint(1, parent1.chromosome_length - 1)

        # Perform single-point crossover
        offspring1 = Organism(None, parent1.chromosome[:crossover_point] + parent2.chromosome[crossover_point:])
        offspring2 = Organism(None, parent2.chromosome[:crossover_point] + parent1.chromosome[crossover_point:])

        return [offspring1, offspring2]
    
    
    def recombination_two_point_crossover(parent1: Organism, parent2: Organism)-> list:
        """
        This method performs a two-point crossover recombination between two parent organisms.
        In two-point crossover, two random crossover points are selected, and the genes of the parents
        are exchanged to create two new offspring.

        ### Parameters
        - parent1 (`Organism`): The first parent organism.
        - parent2 (`Organism`): The second parent organism.

        ### Returns
        - `list`: A list of two offspring organisms.
        """

        # Select two random crossover points
        crossover_points = sorted(random.sample(range(1, parent1.chromosome_length - 1), 2))

        # Perform two-point crossover
        offspring1 = Organism(None, parent1.chromosome[:crossover_points[0]] + parent2.chromosome[crossover_points[0]:crossover_points[1]] + parent1.chromosome[crossover_points[1]:])
        offspring2 = Organism(None, parent2.chromosome[:crossover_points[0]] + parent1.chromosome[crossover_points[0]:crossover_points[1]] + parent2.chromosome[crossover_points[1]:])

        return [offspring1, offspring2]
    

    def recombination_uniform_crossover(parent1: Organism, parent2: Organism, probability: float)-> list:
        """
        This method performs a uniform crossover recombination between two parent organisms.
        In uniform crossover, each gene is selected from one of the parents with a given probability.

        ### Parameters
        - parent1 (`Organism`): The first parent organism.
        - parent2 (`Organism`): The second parent organism.
        - probability (`float`): The probability of selecting a gene from the first parent.

        ### Returns
        - `list`: A list of two offspring organisms.
        """

        if probability is None:
            raise ValueError("The uniform crossover probability must be specified. Use `GeneticAlgorithm.uniform_crossover_probability = probabiltiy` to set the probability.")
        if probability < 0 or probability > 1 or probability:
            raise ValueError("The uniform crossover probability must be a value between 0 and 1. Update the value with `GeneticAlgorithm.uniform_crossover_probability = probability`.")

        # Perform uniform crossover
        offspring1 = Organism(None, [parent1.chromosome[i] if random.random() < probability else parent2.chromosome[i] for i in range(parent1.chromosome_length)])
        offspring2 = Organism(None, [parent2.chromosome[i] if random.random() < probability else parent1.chromosome[i] for i in range(parent1.chromosome_length)])

        return [offspring1, offspring2]
    

    def recombination_multi_point_crossover(parent1: Organism, parent2: Organism, num_points: int)-> list:
        """
        This method performs a multi-point crossover recombination between two parent organisms.
        In multi-point crossover, multiple random crossover points are selected, and the genes of the parents
        are exchanged to create two new offspring.

        ### Parameters
        - parent1 (`Organism`): The first parent organism.
        - parent2 (`Organism`): The second parent organism.
        - num_points (`int`): The number of crossover points.

        ### Returns
        - `list`: A list of two offspring organisms.
        """

        if num_points == 0:
            raise ValueError("The number of crossover points must be greater than 0. Update the value with `GeneticAlgorithm.multi_point_crossover_points = num_points`.")
        if num_points >= parent1.chromosome_length:
            raise ValueError("The number of crossover points must be less than the chromosome length. Update the value with `GeneticAlgorithm.multi_point_crossover_points = num_points`.")

        # Select random crossover points
        crossover_points = sorted(random.sample(range(1, parent1.chromosome_length - 1), num_points))

        # Perform multi-point crossover
        offspring1_chromosome = []
        offspring2_chromosome = []

        # Initialize the parent index
        parent_index = 0

        # Iterate over the chromosome
        for i in range(parent1.chromosome_length):
            # Check if the parent index is odd or even
            if parent_index % 2 == 0:
                offspring1_chromosome.append(parent1.chromosome[i])
                offspring2_chromosome.append(parent2.chromosome[i])
            else:
                offspring1_chromosome.append(parent2.chromosome[i])
                offspring2_chromosome.append(parent1.chromosome[i])

            # Check if the current index is a crossover point
            if i in crossover_points:
                parent_index += 1


        # Create offspring organisms
        offspring1 = Organism(None, offspring1_chromosome)
        offspring2 = Organism(None, offspring2_chromosome)

        return [offspring1, offspring2]
    # +++++++++++++++++++++++++++++++++++++++++++++++++++





    # --------------------------------------------------
    # Private Methods
    # --------------------------------------------------
    # +++++++++++++++++++++++++++++++++++++++++++++++++++