from .Organism import Organism

import random
import numpy as np


class Recombination:

    # --------------------------------------------------
    # Methods
    # --------------------------------------------------
    def recombination_single_crossover(parent1: Organism, parent2: Organism, mutation_probability: float)-> list:
        """
        This method performs a single-point crossover recombination between two parent organisms.
        In single-point crossover, a random crossover point is selected, and the genes of the parents
        are exchanged to create two new offspring.

        ### Parameters
        - parent1 (`Organism`): The first parent organism.
        - parent2 (`Organism`): The second parent organism.
        - mutation_probability (`float`): The probability of mutating the offspring.

        ### Returns
        - `list`: A list of two offspring organisms.
        """

        # Select a random crossover point
        crossover_point = random.randint(1, parent1.chromosome_length - 1)

        # Perform single-point crossover
        offspring1 = Organism(None, parent1.chromosome[:crossover_point] + parent2.chromosome[crossover_point:])
        offspring2 = Organism(None, parent2.chromosome[:crossover_point] + parent1.chromosome[crossover_point:])
        
        return Recombination.__mutate_offspring([offspring1, offspring2], mutation_probability)
    
    
    def recombination_two_point_crossover(parent1: Organism, parent2: Organism, mutation_probability: float)-> list:
        """
        This method performs a two-point crossover recombination between two parent organisms.
        In two-point crossover, two random crossover points are selected, and the genes of the parents
        are exchanged to create two new offspring.

        ### Parameters
        - parent1 (`Organism`): The first parent organism.
        - parent2 (`Organism`): The second parent organism.
        - mutation_probability (`float`): The probability of mutating the offspring.

        ### Returns
        - `list`: A list of two offspring organisms.
        """

        # Select two random crossover points
        crossover_points = sorted(random.sample(range(1, parent1.chromosome_length - 1), 2))

        # Perform two-point crossover
        offspring1 = Organism(None, parent1.chromosome[:crossover_points[0]] + parent2.chromosome[crossover_points[0]:crossover_points[1]] + parent1.chromosome[crossover_points[1]:])
        offspring2 = Organism(None, parent2.chromosome[:crossover_points[0]] + parent1.chromosome[crossover_points[0]:crossover_points[1]] + parent2.chromosome[crossover_points[1]:])

        return Recombination.__mutate_offspring([offspring1, offspring2], mutation_probability)
    

    def recombination_uniform_crossover(parent1: Organism, parent2: Organism, probability: float, mutation_probability: float)-> list:
        """
        This method performs a uniform crossover recombination between two parent organisms.
        In uniform crossover, each gene is selected from one of the parents with a given probability.

        ### Parameters
        - parent1 (`Organism`): The first parent organism.
        - parent2 (`Organism`): The second parent organism.
        - probability (`float`): The probability of selecting a gene from the first parent.
        - mutation_probability (`float`): The probability of mutating the offspring.

        ### Returns
        - `list`: A list of two offspring organisms.
        """

        if probability is None:
            raise ValueError("The uniform crossover probability must be specified. Use `GeneticAlgorithm.uniform_crossover_probability = probabiltiy` to set the probability.")
        if probability < 0 or probability > 1:
            raise ValueError(f"The uniform crossover probability must be a value between 0 and 1. Thed current value is '{probability}' Update the value with `GeneticAlgorithm.uniform_crossover_probability = probability`.")

        # Perform uniform crossover
        offspring1 = Organism(None, [parent1.chromosome[i] if random.random() < probability else parent2.chromosome[i] for i in range(parent1.chromosome_length)])
        offspring2 = Organism(None, [parent2.chromosome[i] if random.random() < probability else parent1.chromosome[i] for i in range(parent1.chromosome_length)])

        return Recombination.__mutate_offspring([offspring1, offspring2], mutation_probability)
    

    def recombination_multi_point_crossover(parent1: Organism, parent2: Organism, num_points: int, mutation_probability: float)-> list:
        """
        This method performs a multi-point crossover recombination between two parent organisms.
        In multi-point crossover, multiple random crossover points are selected, and the genes of the parents
        are exchanged to create two new offspring.

        ### Parameters
        - parent1 (`Organism`): The first parent organism.
        - parent2 (`Organism`): The second parent organism.
        - num_points (`int`): The number of crossover points.
        - mutation_probability (`float`): The probability of mutating the offspring.

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

        return Recombination.__mutate_offspring([offspring1, offspring2], mutation_probability)
    

    def recombination_2d_spatial_crossover(parent1: Organism, parent2: Organism, grid_shape: tuple, subgrid_shape: tuple, mutation_probability: float, gene_lengths: list)-> list:
        """
        This method performs a 2D spatial crossover recombination between two parent organisms.
        In 2D spatial crossover, a subgrid is selected from each parent, and the genes of the parents
        are exchanged to create two new offspring.

        This method is useful for problems where the position of the genes in the chromosome is important,
        such as grid-based problems.

        ### Parameters
        - parent1 (`Organism`): The first parent organism.
        - parent2 (`Organism`): The second parent organism.
        - grid_shape (`tuple`): The shape of the grid.
        - subgrid_shape (`tuple`): The shape of the subgrid.
        - mutation_probability (`float`): The probability of mutating the offspring.
        - gene_lengths (`list`): A list consisting of the lengths of each gene in the chromosome.

        ### Returns
        - `list`: A list of two offspring organisms.
        """

        # Define the grid parameters and subgrid parameters
        n, m = grid_shape
        k, l = subgrid_shape

        # Retrieve the chromosomes of the parent organisms
        chromosome1 = parent1.chromosome
        chromosome2 = parent2.chromosome


        # Check whether the chromosome can be represented as a 2D grid
        if len(gene_lengths) != n * m:
            raise ValueError(f"The number of genes in the chromosome ({len(gene_lengths)}) does not match the grid shape ({n}x{m}). Spatial crossover requires that the number of genes in the chromosome matches the the area of the grid")


        # Initialize a 2D grid with the indices
        # of all the genes in the chromosome
        gene_grid1 = np.array([], dtype=str)
        gene_grid2 = np.array([], dtype=str)

        # Fill the gene grids
        for i in range(len(gene_lengths)):

            # Retrieve the gene length
            gene_length = gene_lengths[i]

            # Retrieve the genes (slice of the chromosome)
            gene1 = "".join([str(gene) for gene in chromosome1[:gene_length]])	
            gene2 = "".join([str(gene) for gene in chromosome2[:gene_length]])

            # Slice the genes off the chromosome
            chromosome1 = chromosome1[gene_length:]
            chromosome2 = chromosome2[gene_length:]

            # Append the genes to the gene grids
            gene_grid1 = np.append(gene_grid1, gene1)
            gene_grid2 = np.append(gene_grid2, gene2)

        # Reshape the genes grids to match the grid shape
        gene_grid1 = gene_grid1.reshape(grid_shape)
        gene_grid2 = gene_grid2.reshape(grid_shape)


        # Select a random subgrid
        x = random.randint(0, n - k + 1)
        y = random.randint(0, m - l + 1)


        # Slice the subgrids from the parent gene grids
        subgrid1 = gene_grid1[y:y + l, x:x + k].copy() # From the first parent
        subgrid2 = gene_grid2[y:y + l, x:x + k].copy() # From the second parent

        # Perform the crossover (swapping the gene subgrids)
        offspring_gene_grid1 = gene_grid1.copy()
        offspring_gene_grid2 = gene_grid2.copy()
        offspring_gene_grid1[y:y + l, x:x + k] = subgrid2
        offspring_gene_grid2[y:y + l, x:x + k] = subgrid1

        # Reflatten the genes grids 
        offspring_gene_grid1 = offspring_gene_grid1.flatten().tolist()
        offspring_gene_grid2 = offspring_gene_grid2.flatten().tolist()
        # These lists now contain the genes (as strings)
        
        #   Convert the genes back to a list of integers (the chromosome)
        offspring_chromosome1 = []
        offspring_chromosome2 = []
        for i in range(len(offspring_gene_grid1)):
            offspring_chromosome1 += list(offspring_gene_grid1[i])
            offspring_chromosome2 += list(offspring_gene_grid2[i])


        # Create offspring organisms
        offspring1 = Organism(None, offspring_chromosome1)
        offspring2 = Organism(None, offspring_chromosome2)

        return Recombination.__mutate_offspring([offspring1, offspring2], mutation_probability)
    # +++++++++++++++++++++++++++++++++++++++++++++++++++





    # --------------------------------------------------
    # Private Methods
    # --------------------------------------------------
    def __mutate_offspring(offspring: list, mutation_rate: float):
        """
        This method mutates the offspring organisms by flipping the value of a random gene.
        This is based on the mutation probability (GeneticAlgorithm._mutation_rate)

        ### Parameters
        - offspring (`list`): A list of offspring organisms.
        - mutation_rate (`float`): The probability of mutating a gene.

        ### Returns
        - `list`: A list of mutated offspring organisms.
        """

        # Iterate over the offspring organisms
        organism: Organism
        for organism in offspring:
            # Iterate over the genes in the chromosome
            for i in range(organism.chromosome_length):
                # Check if the gene should be mutated
                if random.random() < mutation_rate:
                    # Flip the value of the gene
                    organism.chromosome[i] = 1 if organism.chromosome[i] == 0 else 0
        return offspring
    # +++++++++++++++++++++++++++++++++++++++++++++++++++