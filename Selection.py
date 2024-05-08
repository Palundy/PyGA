from Organism import Organism

import random


class Selection:

    # --------------------------------------------------
    # Methods
    # --------------------------------------------------
    def selection_roulette_wheel(population: list)-> "Organism":
        """
        This method selects an organism from the population using the roulette wheel selection method.
        Roulette wheel selection is a proportional selection method based on the concept of
        a roulette wheel, where each individual's selection probability is proportional to its fitness.

        The roulette wheel selection favors individuals with higher fitness values, but
        can lead to premature convergence if the fitness values are significantly different,
        or of there are only a few organisms with high fitness values.

        ### Parameters
        - population (`list`): The population of organisms.

        ### Returns
        - Organism: The selected organism.
        """

        # Calculate the total fitness of the population
        total_fitness = sum([organism.fitness for organism in population])

        # Assign a probability to each organism (based on their fitness)
        probabilities = [organism.fitness / total_fitness for organism in population]

        # Select individuals using a roulette wheel selection method
        return population[Selection.__roulette_wheel(probabilities)]


    def selection_rank_based(population: list)-> "Organism":
        """
        This method selects an organism from the population using the rank-based selection method.
        With rank-based selection, the organisms are ranked based on their fitness values.
        The selection probability of an organism is proportional to its rank.

        Rank-based selection ensures a more uniform selection pressure among individuals,
        which reduces premature convergence. However, it may be less sensitive to differences
        in fitness values, and requires more computational resources.

        ### Parameters
        - population (`list`): The population of organisms.

        ### Returns
        - Organism: The selected organism.
        """

        # Sort the population based on fitness
        ranked_population = sorted(population, key=lambda organism: organism.fitness, reverse=True)

        # Assign a selection probability to each orgamism based on its rank.
        # A linear rank-based selection probability is used.
        # The probability of selection is proportional to the rank of the organism.
        probabilities = [i / len(ranked_population) for i in range(1, len(ranked_population) + 1)]

        # Select individuals using a roulette wheel selection method
        return ranked_population[Selection.__roulette_wheel(probabilities)]
    

    def selection_tournament(population: list, tournament_size: int)-> "Organism":
        """
        This method selects an organism from the population using the tournament selection method.
        Tournament selection involves selecting a subset of the population (tournament) and
        choosing the best organism from that subset.

        The tournament selection method is less biased towards high-fitness organisms,
        which can help maintain genetic diversity and prevent premature convergence.

        ### Parameters
        - population (`list`): The population of organisms.
        - tournament_size (`int`): The size of the tournament subset.

        ### Returns
        - Organism: The selected organism.
        """

        # Randomly select a subset of the population for the tournament
        tournament = random.sample(population, tournament_size)

        # Select the best organism from the tournament
        return max(tournament, key=lambda organism: organism.fitness)
    # +++++++++++++++++++++++++++++++++++++++++++++++++++





    # --------------------------------------------------
    # Private Methods
    # --------------------------------------------------
    def __roulette_wheel(probabilities: list)-> int:
        """
        This method selects an index from a list of probabilities using the roulette wheel selection method.

        ### Parameters
        - probabilities (`list`): A list of probabilities.

        ### Returns
        - int: The selected index.
        """

        # Construct a roulette wheel by diving the interval [0, 1] into segments
        # proportional to the probability of selection for each organism
        # (All the probabilities should sum to 1)
        roulette_wheel = []
        lower_bound = 0.0
        for probability in probabilities:
            roulette_wheel.append((lower_bound, lower_bound + probability))
            lower_bound += probability

        # Spin the roulette wheel
        roulette_wheel_value = random.random()
        for i, (lower_bound, upper_bound) in enumerate(roulette_wheel):
            if lower_bound <= roulette_wheel_value < upper_bound:
                return i

        # If the roulette wheel value is at the end of the wheel, select the last organism
        return len(roulette_wheel) - 1
    # +++++++++++++++++++++++++++++++++++++++++++++++++++