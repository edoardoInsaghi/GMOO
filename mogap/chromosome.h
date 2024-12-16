#ifndef CHROMOSOME_H
#define CHROMOSOME_H

#include<iostream>
#include<vector>
#include<array>
#include<random>
#include<cmath>
#include<algorithm>
#include<set>
#include<iomanip>

constexpr int COLWIDTH = 2;

class Chromosome {
public:
    int NA;  // Number of planes
    int NG;  // Number of gates
    std::vector<int> relative_position_genes; // NA * NA matrix
    std::vector<int> gate_genes;              // NA * 1 vector
    std::set<int> used_gates;

    Chromosome(int NA, int NG, std::vector<int>& relative_position_genes, std::vector<int>& gate_genes);
    Chromosome(const Chromosome& other);
    Chromosome& operator=(const Chromosome& other);
    bool operator==(const Chromosome& other) const;

    void set_used_gates();
    inline int operator()(int i, int j) const {
        return relative_position_genes[i * NA + j];
    };
    inline int gate_of_plane(int j) const {
        return gate_genes[j];
    };
    inline void set(int i, int j){
        relative_position_genes[i * NA + j] = 1;
    };
    inline void unset(int i, int j){
        relative_position_genes[i * NA + j] = 0;
    };
};

Chromosome init_random_chromosome(int NA, int NG);

// Mutations
void apply_mutation_1(Chromosome& chromosome); // swaps two consequent planes in the same queue
void apply_mutation_2(Chromosome& chromosome); // swaps two planes from different queues

// Crossover
Chromosome perform_crossover(Chromosome& parent1, Chromosome& parent2); // uniform crossover

// Constraints
bool constraint_1(Chromosome& chromosome);
bool constraint_2(Chromosome& chromosome);
bool constraint_3(Chromosome& chromosome);
bool constraint_4(Chromosome& chromosome);
bool constraint_5(Chromosome& chromosome);
bool queue_gate_consistency(Chromosome& chromosome);
bool is_valid_chromosome(Chromosome& chromosome);

// Utilities
void show_chromosome_absolute_position(Chromosome& chromosome); 
void show_chromosome(Chromosome& chromosome);
int find_no_nonzero_row(Chromosome& chromosome, int i); // init chromosome, returns index of non 0 element in given row, default -1
bool no_nonzero_row(Chromosome& chromosome, int i); // init chromosome, returns true if row is all zeros (ignores element in the diagonal)
int select_random_non_assigned_plane(Chromosome& chromosome); // crossover, returns random index for which column is all zeros
int select_random_first_in_queue(Chromosome& chromosome); // init chromosome, returns the first plane on a random queue
bool is_column_free(Chromosome& chromosome, int j); // crossover, returns true if all elements in the column are zeros
int select_random_nonzero_gene(Chromosome& chromosome); // mutation 1, returns flattened intex of random non-zero element in the relative position representation
                                                        // to retrieve j: j = index % number of planes
                                                        // to retrieve i: i = index / chromosome.NA
std::array<int, 2> select_2_random_nonzero_genes(Chromosome& chromosome); // mutation 2, returns an array of two distinct elements, see above
int predecessor_in_queue(Chromosome& chromosome, int i); // {mutation 1, perform crossover}, returns predecessor of plane i in its own queue
int successor_in_queue(Chromosome& chromosome, int i); // {mutation 1, mutation 2, other utils}, returns predecessor of plane i in its own queue
std::vector<int> find_first_in_queue(Chromosome& chromosome); // other utils, returns first plane in each gate, -1 if gate is not used
int find_first_in_queue_g(Chromosome& chromosome, int g); // other utils, returns first plane in gate g
int find_last_in_queue(Chromosome& chromosome, int i); // crossover, returns last plane of gate i
int gate_num_first_in_queue(Chromosome& chromosome, int gate); // crossover, returns the number of gates in first position of gate
int find_nonzero_in_column(Chromosome& chromosome, int j); // crossover, returns index of non-zero element in column j

template<typename T>
T select_random_from_set(std::set<T>& set); // utils, selects random element from a std::set

template<typename T>
T select_random_from_set_and_remove(std::set<T>& set); // init chromosome, selects random element from a std::set and removes

template<typename T>
T select_random_from_set(std::set<T>& set) {
    int random_index = rand() % set.size();
    auto it = std::next(set.begin(), random_index);
    if (it == set.end()) {
        throw std::runtime_error("Random index out of bounds\n");
    }
    return *it;
}

template<typename T>
T select_random_from_set_and_remove(std::set<T>& set) {
    int random_index = rand() % set.size();
    auto it = std::next(set.begin(), random_index);
    T element = *it;
    set.erase(it);
    return element;
}

#endif // CHROMOSOME_H