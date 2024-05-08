import random


class Organism:
    """
    This class represents a organism used
    within the `GeneticAlgorithm` class.
    """

    # --------------------------------------------------
    # Variables
    # --------------------------------------------------
    _chromosome = []
    _fitness = 0.0
    # +++++++++++++++++++++++++++++++++++++++++++++++++++





    # --------------------------------------------------
    # Methods
    # --------------------------------------------------
    def __init__(self, chromosome_length: int, chromosome: None|list = None)-> None:
        """
        This method initializes the organism.

        ### Parameters
        - chromosome_length (`int`): The length of the chromosome.
        - chromosome (`None|list`): The chromosome of the organism. If None, a random chromosome is generated.
        """

        # Initialize the organism
        if (chromosome is not None):
            self._chromosome = chromosome
            self._chromosome_length = len(chromosome)
        else:
            self._chromosome_length = chromosome_length
            self._chromosome = self.__chromosome()
        self._fitness = 0.0
        pass


    def __str__(self)-> str:
        """
        This method returns a string representation of the organism.

        ### Returns
        - str: A string representation of the organism.
        """
        return_string = f"[Organism]\n"
        return_string += f"     Chromosome: {''.join([str(gene) for gene in self._chromosome])}\n"
        return_string += f"     Fitness: {self._fitness}\n"
        return return_string
    # +++++++++++++++++++++++++++++++++++++++++++++++++++





    # --------------------------------------------------
    # Private Methods
    # --------------------------------------------------
    def __chromosome(self)-> list:
        """
        Generates a random chromosome for the organism.

        ### Returns
        - list: A list of genes that represent the chromosome of the organism.
        """
        return [random.randint(0, 1) for _ in range(self.chromosome_length)]
    # +++++++++++++++++++++++++++++++++++++++++++++++++++
        




    # --------------------------------------------------
    # Getters and Setters
    @property
    def chromosome_length(self)-> int:
        """
        This property returns the length of the chromosome.

        ### Returns
        - int: The length of the chromosome.
        """
        return self._chromosome_length
    
    
    @property
    def chromosome(self)-> list:
        """
        This property returns the chromosome of the organism.

        ### Returns
        - list: The chromosome of the organism.
        """
        return self._chromosome
    

    @property
    def fitness(self)-> float:
        """
        This property returns the fitness of the organism.

        ### Returns
        - float: The fitness of the organism.
        """
        return self._fitness
    

    @fitness.setter
    def fitness(self, value: float)-> None:
        """
        This property sets the fitness of the organism.

        ### Parameters
        - value (`float`): The fitness of the organism.
        """
        self._fitness = value
    # +++++++++++++++++++++++++++++++++++++++++++++++++++





    # --------------------------------------------------
    # Exceptions
    # --------------------------------------------------
    # +++++++++++++++++++++++++++++++++++++++++++++++++++
