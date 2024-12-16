#include"chromosome.h"
#include"mogap.h"

#include<unordered_set>

template <typename T>
std::vector<T> concatenate(const std::vector<T>& v1, const std::vector<T>& v2) {
    std::vector<T> result;
    result.reserve(v1.size() + v2.size()); 
    result.insert(result.end(), v1.begin(), v1.end());
    result.insert(result.end(), v2.begin(), v2.end()); 
    return result;
}

std::vector<Chromosome> init_population(int population_size, int NA, int NG) {
    std::vector<Chromosome> population;
    for (int i=0; i<population_size; i++) {
        Chromosome chromosome = init_random_chromosome(NA, NG);
        if (!is_valid_chromosome(chromosome)) {
            throw std::runtime_error("Invalid chromosome during initialization");
        }
        population.push_back(std::move(chromosome));
    }
    return population;
}

// Fitness evaluation
std::vector<double> fitness_evaluation(std::vector<Chromosome>& population, AirportGAP& airport) {
    std::vector<double> fitness_values;
    for (auto& chromosome : population) {
        if (!is_valid_chromosome(chromosome)) {
            throw std::runtime_error("Invalid chromosome during fitness evaluation");
        }
        Solution solution = chromosome_to_solution(chromosome);
        double fitness = airport.calculateFitness(solution);
        fitness_values.push_back(fitness);
    }
    return fitness_values;
}

// Selection
std::vector<Chromosome> selection(std::vector<Chromosome>& population, std::vector<double>& fitness_values, int n) {
    std::vector<Chromosome> selected_population;
    std::vector<std::pair<double, int>> fitness_indices;

    if (selected_population.size() != fitness_values.size()) {
        throw std::runtime_error("Mismatch in population and fitness values");
    }

    for (int i=0; i<fitness_values.size(); i++) {
        fitness_indices.emplace_back(fitness_values[i], i);
    }

    std::sort(fitness_indices.begin(), fitness_indices.end(),
              [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                  return a.first > b.first;
              });

    for (int i=0; i<n && i<fitness_indices.size(); i++) {
        selected_population.push_back(population[fitness_indices[i].second]);
    }

    return selected_population;
}

// Mutation
void mutation(std::vector<Chromosome>& population) {
    for (auto& chromosome : population) {
        apply_mutation_1(chromosome);
        apply_mutation_2(chromosome);
    }
}

// Crossover
std::vector<Chromosome> crossover(std::vector<Chromosome>& population, int n) {
    std::vector<Chromosome> offspring;

    for (int i=0; i<n-1; i++) {
        Chromosome& parent1 = population[i];
        Chromosome& parent2 = population[i+1];

        if (!is_valid_chromosome(parent1) || !is_valid_chromosome(parent2)) {
            throw std::runtime_error("Invalid parent chromosomes during crossover");
        }

        Chromosome child = perform_crossover(parent1, parent2);
        if (!is_valid_chromosome(child)) {
            throw std::runtime_error("Invalid child chromosome after crossover");
        }
        offspring.push_back(std::move(child));
    }

    return concatenate(population, offspring);
}

// Genetic Algorithm Main Loop
void GA(int NA, int NG, AirportGAP& airport, int mu, int lambda, int iterations) {

    // Initialize population
    std::vector<Chromosome> population = init_population(lambda, NA, NG);

    // Validate initial population
    for (auto& chromosome : population) {
        if (!is_valid_chromosome(chromosome)) {
            throw std::runtime_error("Invalid chromosome in initial population");
        }
    }

    std::cout << "Initial population is valid and size: " << population.size() << std::endl;

    // Evaluate fitness
    std::vector<double> fitness = fitness_evaluation(population, airport);
    std::cout << "Initial fitness values: " << fitness.size() << std::endl << " Initial population size: " << population.size() << std::endl;   

    for (int iter=0; iter<iterations; iter++) {

        std::vector<Chromosome> selected_population = selection(population, fitness, mu); // Selection

        population = std::move(selected_population);

        mutation(population); // Mutation
        
        population = crossover(population, lambda); // Crossover
    
        fitness = fitness_evaluation(population, airport); // Evaluate fitness

        std::cout << "Iteration " << iter + 1 << ": Population size = " << population.size() << ", Fitness values = " << fitness.size() << std::endl;

        double max_fitness = *std::max_element(fitness.begin(), fitness.end());
        std::cout << "Iteration " << iter + 1 << ": Best fitness = " << max_fitness << std::endl;
    }
}

int main(int argc, char *argv[]) { 

    int NA = 40;
    int NG = 10;

    AirportGAP airport(NA, NG);
    airport.generateAircraftSet(0, 1440);  // 24 hours
    GA(NA, NG, airport, 50, 200, 100);
    
    return 0; 
}
