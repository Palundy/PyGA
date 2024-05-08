from Organism import Organism
from Selection import Selection
from Recombination import Recombination

import random


# Ideas and Notes!
# - Add parameter to make a gene only be interpreted as a integer or float
# - Add parameter to make a gene only be interpreted as a positive or negative number


class GeneticAlgorithm:

    # --------------------------------------------------
    # Constants
    # --------------------------------------------------
    RECOMBINATION_SINGLE_CROSSOVER = 0
    RECOMBINATION_TWO_POINT_CROSSOVER = 1
    RECOMBINATION_UNIFORM_CROSSOVER = 2
    RECOMBINATION_MULTI_POINT_CROSSOVER = 3
    SELECTION_ROULETTE_WHEEL = 0
    SELECTION_RANK = 1
    SELECTION_TOURNAMENT = 2
    # +++++++++++++++++++++++++++++++++++++++++++++++++++





    # --------------------------------------------------
    # Variables
    # --------------------------------------------------
    _fitness_function = None
    _population_size = 0
    _population_initialised = False
    _population_evaluated = False
    _population = []
    _gene_count = 0
    _genes = {}
    _chromosome_length = 0
    _recombination_method = None
    _selection_method = None
    _selection_size = 0
    _generation_count = 0
    _tournament_size = 0
    _multi_point_num_points = 0
    __uniform_crossover_probability = None
    # +++++++++++++++++++++++++++++++++++++++++++++++++++





    # --------------------------------------------------
    # Methods
    # --------------------------------------------------
    def __init__(self, fitness_function: callable):
        self._fitness_function = fitness_function
        pass


    def add_gene(self, alias: str, bottom_value: int, top_value: int, bitsize: int)-> "GeneticAlgorithm":
        """
        This method adds a gene to the genetic algorithm.

        ### Parameters
        - alias (`str`): The alias of the gene.
        - bottom_value (`int`): The bottom value of the gene.
        - top_value (`int`): The top value of the gene. 
        - bitsize (`int`): The bitsize of the gene. This represents the amount of bits used to represent the gene. Given that bitsize is `n`, the gene can have `2ⁿ` possible values. Minimum value of the bitsize 1.

        ### Returns
        - `GeneticAlgorithm`: This instance of the genetic algorithm.

        ### Exceptions
        - `InvalidGeneValue`: Raised when an invalid gene value is set.
        """

        # Check whether the values are valid
        if bottom_value < 0:
            raise GeneticAlgorithm.InvalidGeneValue("The bottom_value cannot be lower than 0.")
        if top_value < 0:
            raise GeneticAlgorithm.InvalidGeneValue("The top_value cannot be lower than 0.")
        if bitsize < 1:
            raise GeneticAlgorithm.InvalidGeneValue("The bitsize cannot be lower than 1.")
        if bottom_value > top_value:
            raise GeneticAlgorithm.InvalidGeneValue("The bottom_value cannot be greater than the top_value.")
        
        # Check whether there is a gene with the same alias
        if alias in self._genes:
            raise GeneticAlgorithm.InvalidGeneValue("The gene alias already exists.")
        
        
        # The bit positions of the gene should be calculated
        # this depends on the number of genes that have been added before this one
        bottom_index = self._chromosome_length
        top_index = self._chromosome_length + bitsize


        # Add the gene to the genetic algorithm
        self._genes[alias] = (bottom_index, top_index, bottom_value, top_value, bitsize)

        # Update the chromosome length
        self._chromosome_length += bitsize

        # Return this instance of the GeneticAlgorithm
        return self
    

    def recombination_method(self, method: int)-> "GeneticAlgorithm":
        """
        This method sets the recombination method to be used
        by the genetic algorithm.

        ### Parameters
        - method (`int`): The recombination method to be used.

        ### Returns
        - `GeneticAlgorithm`: This instance of the genetic algorithm.

        ### Exceptions
        - `InvalidRecombinationMethod`: Raised when an invalid recombination method is set.
        """

        # Check whether the method is valid
        if method not in [
            GeneticAlgorithm.RECOMBINATION_SINGLE_CROSSOVER,
            GeneticAlgorithm.RECOMBINATION_TWO_POINT_CROSSOVER,
            GeneticAlgorithm.RECOMBINATION_UNIFORM_CROSSOVER,
            GeneticAlgorithm.RECOMBINATION_MULTI_POINT_CROSSOVER
        ]:
            raise GeneticAlgorithm.InvalidRecombinationMethod("Invalid recombination method.")
        
        # Set the recombination method
        self._recombination_method = method

        # Return this instance of the GeneticAlgorithm
        return self
    

    def selection_method(self, method: int)-> "GeneticAlgorithm":
        """
        This method sets the selection method to be used
        by the genetic algorithm.

        ### Parameters
        - method (`int`): The selection method to be used.

        ### Returns
        - `GeneticAlgorithm`: This instance of the genetic algorithm.

        ### Exceptions
        - `InvalidSelectionMethod`: Raised when an invalid selection method is set.
        """

        # Check whether the method is valid
        if method not in [
            GeneticAlgorithm.SELECTION_ROULETTE_WHEEL,
            GeneticAlgorithm.SELECTION_RANK,
            GeneticAlgorithm.SELECTION_TOURNAMENT
        ]:
            raise GeneticAlgorithm.InvalidSelectionMethod("Invalid selection method.")
        
        # Set the selection method
        self._selection_method = method

        # Return this instance of the GeneticAlgorithm
        return self


    def generate_population(self, population_size: int)-> "GeneticAlgorithm":
        """
        This method generates a population of organisms.

        ### Parameters
        - population_size (`int`): The size of the population.

        ### Returns
        - `GeneticAlgorithm`: This instance of the genetic algorithm.

        ### Exceptions
        - `GeneticAlgorithmException`: Raised when the population has already been initialised.
        - `InvalidGeneValue`: Raised when no genes have been added to the genetic algorithm.
        """

        if self._population_initialised:
            raise GeneticAlgorithm.GeneticAlgorithmException("The population has already been initialised.")
        if len(self._genes) == 0:
            raise GeneticAlgorithm.InvalidGeneValue("No genes have been added to the genetic algorithm. Use the `GeneticAlgorithm.add_gene()` method to add genes.")

        # Generate the population
        for i in range(population_size):
            self._population.append(Organism(self._chromosome_length))
        
        # Set the current population size
        self._population_size = population_size

        # Update the population initialised flag
        self._population_initialised = True

        # Update the generation count
        self._generation_count += 1

        # Return this instance of the GeneticAlgorithm
        return self
    

    def generate_next_population(self)-> "GeneticAlgorithm":
        """
        This method generates the next population of organisms.
        
        The selection method that has been set is used to select
        organisms from the current population. And the selection size
        that has been set is used to determine the number of organisms
        to reproduce the next population.
        
        The recombination method that has been set is used to recombine the selected
        organisms.

        This method replaces the current population with the new population.

        ### Returns
        - `GeneticAlgorithm`: This instance of the genetic algorithm.
        """

        if not self._population_initialised:
            raise GeneticAlgorithm.GeneticAlgorithmException("The population has not been initialised. Use the `GeneticAlgorithm.generate_population()` method to generate a population.")
        if not self._population_evaluated:
            raise GeneticAlgorithm.GeneticAlgorithmException("The population has not been evaluated. Use the `GeneticAlgorithm.evaluate_population()` method to evaluate the fitness of the population.")
        if len(self._genes) == 0:
            raise GeneticAlgorithm.InvalidGeneValue("No genes have been added to the genetic algorithm. Use the `GeneticAlgorithm.add_gene()` method to add genes.")
        if self._selection_method is None:
            raise GeneticAlgorithm.GeneticAlgorithmException("No selection method has been set. Use the `GeneticAlgorithm.selection_method()` method to set a selection method.")
        if self._recombination_method is None:
            raise GeneticAlgorithm.GeneticAlgorithmException("No recombination method has been set. Use the `GeneticAlgorithm.recombination_method()` method to set a recombination method.")
        
        # Select organisms to reproduce
        selected_organisms = self.__select()
        
        # Recombine the organisms
        new_population = []
        new_population_size = 0
        while new_population_size < self._population_size:
            for i in range(0, len(selected_organisms), 2):
                organism1 = selected_organisms[i]
                if i + 1 < len(selected_organisms):
                    organism2 = selected_organisms[i + 1]

                offspring = self.__recombine(organism1, organism2)
                # Check whether both offspring can be added to the new population
                if new_population_size + 2 <= self._population_size:
                    new_population.extend(offspring)
                    new_population_size += 2
                else:
                    # Add the first offspring to the new population
                    new_population.append(offspring[0])
                    new_population_size += 1
                    break
        
        # Set the new population
        self._population = new_population

        # Update the population evaluated flag
        self._population_evaluated = False

        # Update the generation count
        self._generation_count += 1

        # Return this instance of the GeneticAlgorithm
        return self
    

    def evaluate_population(self)-> "GeneticAlgorithm":
        """
        This method evaluates the fitness of the population.

        ### Returns
        - `GeneticAlgorithm`: This instance of the genetic algorithm.
        """

        if not self._population_initialised:
            raise GeneticAlgorithm.GeneticAlgorithmException("The population has not been initialised. Use the `GeneticAlgorithm.generate_population()` method to generate a population.")
        if len(self._genes) == 0:
            raise GeneticAlgorithm.InvalidGeneValue("No genes have been added to the genetic algorithm. Use the `GeneticAlgorithm.add_gene()` method to add genes.")

        # Evaluate the fitness of the population
        for organism in self._population:
            organism.fitness = self.__fitness(organism)

        # Update the population evaluated flag
        self._population_evaluated = True
        
        # Return this instance of the GeneticAlgorithm
        return self
    

    def mean_fitness(self)-> float:
        """
        This method calculates the mean fitness of the population.

        ### Returns
        - float: The mean fitness of the population.
        """

        if not self._population_evaluated:
            raise GeneticAlgorithm.GeneticAlgorithmException("The population has not been evaluated. Use the `GeneticAlgorithm.evaluate_population()` method to evaluate the fitness of the population.")

        # Calculate the mean fitness of the population
        return sum([organism.fitness for organism in self._population]) / self._population_size
    

    def best_organism(self)-> Organism:
        """
        This method returns the best organism in the population.

        ### Returns
        - `Organism`: The best organism in the population.
        """

        if not self._population_evaluated:
            raise GeneticAlgorithm.GeneticAlgorithmException("The population has not been evaluated. Use the `GeneticAlgorithm.evaluate_population()` method to evaluate the fitness of the population.")

        # Find the best organism in the population
        return max(self._population, key=lambda organism: organism.fitness)
    

    def describe_organism(self, organism: Organism)-> str:
        """
        This method describes an organism, by listing the properties
        of the organism in a readable format.

        ### Parameters
        - organism (`Organism`): The organism to describe.

        ### Returns
        - str: A string representation of the organism.
        """
        
        # Convert the genes to properties
        properties = {}

        # Retrieve the genes
        genes = organism.chromosome

        # Retrieve the gene properties
        gene_value = 0
        for gene_alias, gene in self._genes.items():
            # Unpack the gene
            bottom_index, top_index, bottom_value, top_value, bitsize = gene

            # Slice the gene from the chromosome
            gene_part = genes[bottom_index:top_index]

            # Convert the gene (list) to a binary string
            gene_part = ''.join([str(gene) for gene in gene_part])

            # Convert the binary string to an integer
            gene_value = int(gene_part, 2)

            # Calculate the value of the gene
            gene_value = bottom_value + (top_value - bottom_value) * gene_value / (2**bitsize)
        
            # Add the gene value to the properties
            properties[gene_alias] = gene_value

        # Format the properties
        return_string = "[Organism]\n"
        for key, value in properties.items():
            return_string += f"     {key}: {value}\n"
        return_string += f"     Fitness: {organism.fitness}\n"

        return return_string
    # +++++++++++++++++++++++++++++++++++++++++++++++++++





    # --------------------------------------------------
    # Private Methods
    # --------------------------------------------------
    def __select(self)-> list:
        """
        This method selects organisms from the population.
        The selection method that has been set is used to select
        organisms from the population.

        ### Returns
        - list: A list of selected organisms.
        """

        # Check whether the population size has been set
        if self._selection_size == 0:
            raise GeneticAlgorithm.GeneticAlgorithmException("The selection size has not been set. Use the `GeneticAlgorithm.selection_size()` method to set the selection size.")

        selection = []
        for i in range(self._selection_size):
            if self._selection_method == GeneticAlgorithm.SELECTION_ROULETTE_WHEEL:
                selection.append(Selection.selection_roulette_wheel(self._population))
                continue
            
            if self._selection_method == GeneticAlgorithm.SELECTION_RANK:
                selection.append(Selection.selection_rank_based(self._population))
                continue
            
            if self._selection_method == GeneticAlgorithm.SELECTION_TOURNAMENT:
                selection.append(Selection.selection_tournament(self._population, self._tournament_size))
                continue
                # Note to self, a method to set the tournament size should be added

            raise GeneticAlgorithm.GeneticAlgorithmException("Invalid selection method.")
        
        return selection
    
            
    def __recombine(self, organism1: Organism, organism2: Organism)-> list[Organism, Organism]:
        """
        This method recombines two organisms.
        The recombination method that has been set is used
        to recombine the two organisms.

        ### Parameters
        - organism1 (`Organism`): The first parent organism.
        - organism2 (`Organism`): The second parent organism.

        ### Returns
        - `list`: A list of two offspring organisms.
        """

        if self._recombination_method is None:
            raise GeneticAlgorithm.GeneticAlgorithmException("No recombination method has been set. Use the `GeneticAlgorithm.recombination_method()` method to set a recombination method.")
        if not self._population_initialised:
            raise GeneticAlgorithm.GeneticAlgorithmException("The population has not been initialised.")
        if len(self._genes) == 0:
            raise GeneticAlgorithm.InvalidGeneValue("No genes have been added to the genetic algorithm.")

        # Recombine the organisms
        if self._recombination_method == GeneticAlgorithm.RECOMBINATION_SINGLE_CROSSOVER:
            return Recombination.recombination_single_crossover(organism1, organism2)
        
        if self._recombination_method == GeneticAlgorithm.RECOMBINATION_TWO_POINT_CROSSOVER:
            return Recombination.recombination_two_point_crossover(organism1, organism2)
        
        if self._recombination_method == GeneticAlgorithm.RECOMBINATION_UNIFORM_CROSSOVER:
            return Recombination.recombination_uniform_crossover(organism1, organism2, self._uniform_crossover_probability)
        
        if self._recombination_method == GeneticAlgorithm.RECOMBINATION_MULTI_POINT_CROSSOVER:
            return Recombination.recombination_multi_point_crossover(organism1, organism2, self._multi_point_num_points)
        
        raise GeneticAlgorithm.GeneticAlgorithmException("Invalid recombination method.")
    

    def __fitness(self, organism: Organism)-> any:
        """
        This method calculates the fitness of an organism.
        The user given fitness function is used to calculate
        the fitness of the organism.

        ### Parameters
        - organism (`Organism`): The organism to calculate the fitness for.

        ### Returns
        - `any`: Anything that gets returned by the user given fitness function.
        """
        # Convert the genes to properties
        properties = {}

        # Retrieve the genes
        genes = organism.chromosome

        # Retrieve the gene properties
        gene_value = 0
        for gene_alias, gene in self._genes.items():
            # Unpack the gene
            bottom_index, top_index, bottom_value, top_value, bitsize = gene

            # Slice the gene from the chromosome
            gene_part = genes[bottom_index:top_index]

            # Convert the gene (list) to a binary string
            gene_part = ''.join([str(gene) for gene in gene_part])

            # Convert the binary string to an integer
            gene_value = int(gene_part, 2)

            # Calculate the value of the gene
            gene_value = bottom_value + (top_value - bottom_value) * gene_value / (2**bitsize)
        
            # Add the gene value to the properties
            properties[gene_alias] = gene_value

        # Calculate the fitness
        return self._fitness_function(properties)
    # +++++++++++++++++++++++++++++++++++++++++++++++++++





    # --------------------------------------------------
    # Getters and Setters
    # --------------------------------------------------
    @property
    def selection_size(self):
        """
        This property returns the size of the selection.

        ### Returns
        - int: The size of the selection.
        """
        return self._selection_size
    
    
    @selection_size.setter
    def selection_size(self, value: int):
        """
        This property sets the size of the selection.

        ### Parameters
        - value (`int`): The size of the selection.
        """
        self._selection_size = value


    @property
    def population(self):
        """
        This property returns a string representation of the population.

        ### Returns
        - str: A string representation of the population.
        """
        return_value = ""
        # Show the chromosomes of all the organisms
        # with the index of the organism and the fitness at the end
        for i, organism in enumerate(self._population):
            return_value += f"{i}:  {''.join([str(gene) for gene in organism.chromosome])}  -  F({organism.fitness})\n"
        return return_value
    
    @property
    def organisms(self):
        """
        This property returns the organisms in the population.

        ### Returns
        - list: The organisms in the population.
        """
        return self._population
    

    @property
    def tournament_size(self)-> int:
        """
        This property returns the size of the tournament for the tournament selection method.

        ### Returns
        - int: The size of the tournament.
        """
        return self._tournament_size
    

    @tournament_size.setter
    def tournament_size(self, value: int):
        """
        This property sets the size of the tournament for the tournament selection method.

        ### Parameters
        - value (`int`): The size of the tournament.
        """
        self._tournament_size = value


    @property
    def multi_point_num_points(self)-> int:
        """
        This property returns the number of crossover points for the multi-point crossover method.

        ### Returns
        - int: The number of crossover points.
        """
        return self._multi_point_num_points
    
    
    @multi_point_num_points.setter
    def multi_point_num_points(self, value: int):
        """
        This property sets the number of crossover points for the multi-point crossover method.

        ### Parameters
        - value (`int`): The number of crossover points.
        """
        self._multi_point_num_points = value


    @property
    def uniform_crossover_probability(self)-> float:
        """
        This property returns the probability of selecting a gene from the first parent in uniform crossover.

        ### Returns
        - float: The probability of selecting a gene from the first parent.
        """
        return self.__uniform_crossover_probability
    

    @uniform_crossover_probability.setter
    def uniform_crossover_probability(self, value: float):
        """
        This property sets the probability of selecting a gene from the first parent in uniform crossover.

        ### Parameters
        - value (`float`): The probability of selecting a gene from the first parent.
        """
        self.__uniform_crossover_probability = value
    # +++++++++++++++++++++++++++++++++++++++++++++++++++





    # --------------------------------------------------
    # Exceptions
    # --------------------------------------------------
    class GeneticAlgorithmException(Exception):
        """
        This exception is raised when an error occurs
        in the genetic algorithm.
        """
        pass
    class InvalidRecombinationMethod(Exception):
        """
        This exception is raised when an invalid recombination
        method is set for the genetic algorithm.
        """
        pass
    class InvalidSelectionMethod(Exception):
        """
        This exception is raised when an invalid selection
        method is set for the genetic algorithm.
        """
        pass
    class InvalidGeneValue(Exception):
        """
        This exception is raised when an invalid gene value
        is set for the genetic algorithm.
        """
        pass
    # +++++++++++++++++++++++++++++++++++++++++++++++++++

