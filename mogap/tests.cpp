#include "chromosome.h"

int main(int argc, char *argv[]) {

  int NA = 40;
  int NG = 10;
  bool valid = true;

  for (int i = 0; i < 100; i++) {
    Chromosome c = init_random_chromosome(NA, NG);
    for (int j = 0; j < 100; j++) {
      apply_mutation_1(c);
      apply_mutation_2(c);
      valid = valid && is_valid_chromosome(c);
      if (!valid) {
        throw std::runtime_error("Invalid chromosome after mutation. ❌");
      }
    }
  }
  std::cout << "All mutation tests passed. ✅" << std::endl;

  for (int i = 0; i < 100; i++) {
    Chromosome c1 = init_random_chromosome(NA, NG);
    Chromosome c2 = init_random_chromosome(NA, NG);
    for (int j = 0; j < 100; j++) {
      Chromosome c = perform_crossover(c1, c2);
      valid = valid && is_valid_chromosome(c);
      if (!valid) {
        throw std::runtime_error("Invalid chromosome after crossover. ❌");
      }
    }
  }
  std::cout << "All crossover tests passed. ✅" << std::endl;

    for (int i = 0; i < 1000; i++) {
    Chromosome c1 = init_random_chromosome(NA, NG);
    Chromosome c2 = init_random_chromosome(NA, NG);
    for (int j = 0; j < 100; j++) {
      apply_mutation_1(c1);
      apply_mutation_2(c1);
      apply_mutation_1(c2);
      apply_mutation_2(c2);
      valid = valid && is_valid_chromosome(c1);
      valid = valid && is_valid_chromosome(c2);
      if (!valid) {
        throw std::runtime_error("Invalid chromosome after repeated mutation. ❌");
      }
      Chromosome c1 = perform_crossover(c1, c2);
      Chromosome c2 = perform_crossover(c2, c1);
      valid = valid && is_valid_chromosome(c1);
      if (!valid) {
        throw std::runtime_error("Invalid chromosome after repeated crossover 1. ❌");
      }
      valid = valid && is_valid_chromosome(c2);
      if (!valid) {
        throw std::runtime_error("Invalid chromosome after repeated crossover 2. ❌");
      }
      if (!valid) {
        throw std::runtime_error("Invalid chromosome after repeated mutation-crossover. ❌");
      }
    }
  }
  std::cout << "All hybrid mutation-crossover tests passed. ✅" << std::endl;

  return 0;
}
