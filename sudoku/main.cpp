#include <array>
#include <random>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <unordered_set>
#include <omp.h>
#include <fstream>
#include <sstream>

class Sudoku {
public:
    static constexpr int SIZE = 81;
    static constexpr int ROW_SIZE = 9;
    static constexpr int NUM_BLOCKS = 9;
    static constexpr int BLOCK_SIZE = 9;

    std::array<int, SIZE> board{}; 
    std::array<bool, SIZE> fixed{}; 

    // Check if a number can be placed at a given position.
    bool isValid(int row, int col, int num) const {
        for (int i = 0; i < ROW_SIZE; ++i) {

            // Check row and column
            if (board[row * ROW_SIZE + i] == num || board[i * ROW_SIZE + col] == num)
                return false;

            // Check subgrid
            int subgridRow = (row / 3) * 3 + i / 3;
            int subgridCol = (col / 3) * 3 + i % 3;
            if (board[subgridRow * ROW_SIZE + subgridCol] == num)
                return false;
        }
        return true;
    }

    // Recursive solver
    bool solve(int index = 0) {
        if (index == SIZE)
            return true;

        if (board[index] != 0)
            return solve(index + 1);

        for (int num = 1; num <= 9; ++num) {
            int row = index / ROW_SIZE;
            int col = index % ROW_SIZE;
            if (isValid(row, col, num)) {
                board[index] = num;
                if (solve(index + 1)) return true;
                board[index] = 0;
            }
        }
        return false;
    }

    // Constructor for random Sudoku with n fixed values.
    Sudoku(int n) {
        if (n < 0 || n > SIZE) {
            throw std::invalid_argument("Fixed count must be between 0 and 81.");
        }

        solve();

        std::array<int, SIZE> indices;
        std::iota(indices.begin(), indices.end(), 0);
        std::shuffle(indices.begin(), indices.end(), std::mt19937{std::random_device{}()});

        for (int i=0; i<SIZE; ++i) {
            if (i < n) {
                fixed[indices[i]] = true;
            } 
            else {
                board[indices[i]] = 0;
            }
        }
    }

    // Constructor for random Sudoku with random initial values and given fixed values. Zero values are replaced with random numbers
    // No repetitions are allowed within blocks by construction
    Sudoku(const std::array<int, SIZE>& initialBoard, const std::array<bool, SIZE>& fixedValues) : fixed(fixedValues) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(1, 9);

        board.fill(0);

        for (int i=0; i<SIZE; ++i) {
            if (fixedValues[i]) {
                board[i] = initialBoard[i];
            }
        }

        for (int i=0; i<NUM_BLOCKS; ++i) {

            std::array<std::pair<int, int>, BLOCK_SIZE> blockIndices = getBlockIdx(i);
            std::unordered_set<int> usedInBlock;

            for (const auto& [row, col] : blockIndices) {
                int idx = row * ROW_SIZE + col;
                if (fixedValues[idx]) {
                    usedInBlock.insert(board[idx]);
                }
            }

            std::vector<int> remainingValues;
            for (int val = 1; val <= 9; ++val) {
                if (usedInBlock.find(val) == usedInBlock.end()) {
                    remainingValues.push_back(val);
                }
            }

            std::shuffle(remainingValues.begin(), remainingValues.end(), gen);

            int valueIndex = 0;
            for (const auto& [row, col] : blockIndices) {
                int idx = row * ROW_SIZE + col;
                if (!fixedValues[idx]) {
                    board[idx] = remainingValues[valueIndex++];
                }
            }
        }
    }

    static std::array<std::pair<int, int>, BLOCK_SIZE> getBlockIdx(int block) {
        std::array<std::pair<int, int>, BLOCK_SIZE> blockIndices;
        int row = block / 3;
        int col = block % 3;
        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                blockIndices[i * 3 + j] = std::make_pair(row * 3 + i, col * 3 + j);
            }
        }
        return blockIndices;
    }

    int get(int row, int col) const {
        return board[row * ROW_SIZE + col];
    }

    void set(int row, int col, int val) {
        board[row * ROW_SIZE + col] = val;
    }

    void print() const {
        for (int i=0; i<ROW_SIZE; ++i) {
            for (int j=0; j<ROW_SIZE; ++j) {
                int val = board[i * ROW_SIZE + j];
                std::cout << (val == 0 ? "." : std::to_string(val)) << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n\n";
    }

    void printFixed() const {
        for (int i=0; i<ROW_SIZE; ++i) {
            for (int j=0; j<ROW_SIZE; ++j) {
                int val = fixed[i * ROW_SIZE + j];
                std::cout << (val ? "1" : "0") << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n\n";
    }
};

int penalty(const Sudoku& sudoku) {
    int penalty = 0;
    for (int i=0; i<Sudoku::ROW_SIZE; ++i) {
        for (int j=0; j<Sudoku::ROW_SIZE-1; ++j) {
            for (int k = j+1; k < Sudoku::ROW_SIZE; ++k) {

                // Check row
                if (sudoku.get(i, j) == sudoku.get(i, k)) {
                    ++penalty;
                }
            }
        }
    }
    for (int i=0; i<Sudoku::ROW_SIZE-1; ++i) {
        for (int j=0; j<Sudoku::ROW_SIZE; ++j) {
            for (int k = i+1; k < Sudoku::ROW_SIZE; ++k) {

                // Check column
                if (sudoku.get(i, j) == sudoku.get(k, j)) {
                    ++penalty;
                }
            }
        }
    } 
    return penalty; 
}

double fitness(const Sudoku& sudoku) {
    return 1.0 / (penalty(sudoku) + 1);
}


void mutate(Sudoku& candidate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, Sudoku::NUM_BLOCKS-1);
    std::uniform_int_distribution<> num(0, Sudoku::BLOCK_SIZE-1);

    bool valid = false;
    while (!valid) {

        int block = dis(gen);
        std::array<std::pair<int, int>, Sudoku::BLOCK_SIZE> blockIndices = Sudoku::getBlockIdx(block);

        std::pair<int, int> cell1 = blockIndices[num(gen)];
        std::pair<int, int> cell2 = blockIndices[num(gen)];

        int i1 = cell1.first;
        int j1 = cell1.second;
        int i2 = cell2.first;
        int j2 = cell2.second;

        bool sameCell = (i1 == i2) && (j1 == j2);
        bool fixedNumbers = (candidate.fixed[i1 * candidate.ROW_SIZE + j1] || candidate.fixed[i2 * candidate.ROW_SIZE + j2]);
        if (sameCell || fixedNumbers) continue;

        int tmp = candidate.get(i1, j1);
        candidate.set(i1, j1, candidate.get(i2, j2));
        candidate.set(i2, j2, tmp);

        valid = true;
    }
}



Sudoku crossover(const Sudoku& sudoku1, 
                 const Sudoku& sudoku2) {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1, 2);

    std::array<int, Sudoku::SIZE> childBoard;
    std::array<bool, Sudoku::SIZE> childFixed = sudoku1.fixed;
    Sudoku child(childBoard, childFixed);
    
    for (int i=0; i<Sudoku::NUM_BLOCKS; ++i) {
        std::array<std::pair<int, int>, Sudoku::BLOCK_SIZE> blockIndices = Sudoku::getBlockIdx(i);

        if (dis(gen) == 1) {
            for (int j=0; j<Sudoku::BLOCK_SIZE; ++j) {
                int ii = blockIndices[j].first;
                int jj = blockIndices[j].second;
                int val = sudoku1.get(ii, jj);
                child.set(ii, jj, val);
            }
        } 
        else {
            for (int j=0; j<Sudoku::BLOCK_SIZE; ++j) {
                int ii = blockIndices[j].first;
                int jj = blockIndices[j].second;
                int val = sudoku2.get(ii, jj);
                child.set(ii, jj, val);
            }
        }
    }

    return child;
}


Sudoku kTournament(std::vector<Sudoku>& generation, 
                   std::vector<double>& fitnesses, 
                   int k) {

    int N = (int)generation.size();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> indexDist(0, N-1);

    std::vector<int> idx;
    for (int i=0; i<k; ++i) { // Sampling k individuals from the population
        idx.push_back(indexDist(gen));
    }

    std::sort(idx.begin(), idx.end(), [&](int a, int b) { return fitnesses[a] > fitnesses[b]; });

    Sudoku candidate = generation[idx[0]];
    return candidate;
}


int sudoku_GA(const Sudoku& sudoku, 
              const double childMutationRate, 
              const double parentMutationRate, 
              const int N, 
              const int maxIter, 
              const int k,
              const bool purgeOp,
              std::ofstream& outMax,
              std::ofstream& outAvg) {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> childMutationProb(0, 1);
    std::uniform_real_distribution<> parentMutationProb(0, 1);

    std::array<bool, Sudoku::SIZE> fixed = sudoku.fixed;
    std::vector<Sudoku> population;
    std::vector<double> fitnesses;

    for (int i=0; i<N; i++) { // init and eval population
        Sudoku candidate(sudoku.board, fixed);
        double fit = fitness(candidate);
        population.push_back(candidate);
        fitnesses.push_back(fit);
    }

    std::vector<Sudoku> Generation = population;
    std::vector<Sudoku> NewGeneration;
    std::vector<double> newFitnesses;

    auto bestIdx = std::max_element(fitnesses.begin(), fitnesses.end());
    double bestFitness = *bestIdx;
    Sudoku bestSolution = Generation[std::distance(fitnesses.begin(), bestIdx)];
    double lastBestFitness = 0.0;
    int stepsSinceLastImprovement = 0;

    int iter = 1;
    while (iter <= maxIter) { // && bestFitness < 1.0) {  used for computing average iterations across runs

        newFitnesses.clear();
        while (NewGeneration.size() < N) {

            Sudoku parent1 = kTournament(Generation, fitnesses, k);
            Sudoku parent2 = kTournament(Generation, fitnesses, k);

            Sudoku child = crossover(parent1, parent2);
            if (childMutationProb(gen) < childMutationRate) {
                mutate(child);
            }
            NewGeneration.push_back(child);
            newFitnesses.push_back(fitness(child));

            if (parentMutationProb(gen) < parentMutationRate) {
                mutate(parent1);
            }
            double fit1 = fitness(parent1);

            if (parentMutationProb(gen) < parentMutationRate) {
                mutate(parent2);
            }
            double fit2 = fitness(parent2);

            fit1 > fit2 ? NewGeneration.push_back(parent1) : NewGeneration.push_back(parent2);  
            newFitnesses.push_back(fit1 > fit2 ? fit1 : fit2);          
        } 

        Generation = NewGeneration;
        Generation.push_back(bestSolution);
        NewGeneration.clear();

        // find best candidate for current solution
        bestIdx = std::max_element(newFitnesses.begin(), newFitnesses.end());
        double currentBestScore = *bestIdx;
        Sudoku currentBest = Generation[std::distance(newFitnesses.begin(), bestIdx)];
        if (currentBestScore > bestFitness) {
            bestFitness = currentBestScore;
            bestSolution = currentBest;
        }
        int bestErr = penalty(bestSolution);

        // Purge operator, each individual performs crossover with the best candidate of the generation
        if (stepsSinceLastImprovement > 200 and purgeOp) {
            for (auto& candidate : Generation) {
                Sudoku child = crossover(candidate, bestSolution);
                NewGeneration.push_back(child);
            }
            Generation = NewGeneration;
            NewGeneration.clear();
        }
        // Evaluate improvement
        for (int i=0; i<N; i++) {
            fitnesses[i] = fitness(Generation[i]);
        }
        double avgFitness = std::accumulate(fitnesses.begin(), fitnesses.end(), 0.0) / fitnesses.size();
        bestFitness = *std::max_element(fitnesses.begin(), fitnesses.end());

        outMax << iter << "," << bestFitness << "," << bestErr << "\n";
        outAvg << iter << "," << avgFitness << "\n";

        if (bestFitness > lastBestFitness) {
            stepsSinceLastImprovement = 0;
        } else {
            stepsSinceLastImprovement++;
        }    
        lastBestFitness = bestFitness;
        iter++;
    }

    return iter;
}


std::string generateFilename(const std::string& prefix, int fixedNumbers, bool purgeOp, int threadID) {
    std::ostringstream oss;
    oss << prefix << "_fixed" << fixedNumbers
        << "_purge" << (purgeOp ? "true" : "false")
        << "_thread" << threadID << ".csv";
    return oss.str();
}



int main() {
    double childMutationRate = 0.05;
    double parentMutationRate = 0.5;
    int N = 700;
    int maxIter = 200;
    int k = 7;

    int numSim = 30; 
    std::vector<int> fixedValues = {10, 20, 30, 40};
    std::vector<bool> purgeOptions = {true, false}; 

    #pragma omp parallel for schedule(dynamic)
    for (int configIndex = 0; configIndex < fixedValues.size() * purgeOptions.size(); ++configIndex) {
        int threadID = omp_get_thread_num();

        int fixedIdx = configIndex / purgeOptions.size();
        int purgeIdx = configIndex % purgeOptions.size();

        int fixedNumber = fixedValues[fixedIdx];
        bool purgeOp = purgeOptions[purgeIdx];

        std::string maxFilename = generateFilename("max", fixedNumber, purgeOp, threadID);
        std::string avgFilename = generateFilename("avg", fixedNumber, purgeOp, threadID);

        std::ofstream outMax(maxFilename);
        std::ofstream outAvg(avgFilename);

        if (!outMax.is_open() || !outAvg.is_open()) {
            #pragma omp critical
            {
                std::cerr << "Failed to open output files for thread " << threadID << "\n";
            }
            continue;
        }

        outMax << "Iter,BestFit,BestErr\n";
        outAvg << "Iter,AvgFit\n";

        for (int sim = 0; sim < numSim; ++sim) {
            Sudoku sudoku(fixedNumber);
            sudoku_GA(sudoku, childMutationRate, parentMutationRate, N, maxIter, k, purgeOp, outMax, outAvg);
        }

        outMax.close();
        outAvg.close();

        #pragma omp critical
        {
            std::cout << "Completed: Fixed=" << fixedNumber
                      << ", Purge=" << purgeOp
                      << ", Thread=" << threadID << "\n";
        }
    }

    return 0;
}
 

/*
int main() {

    double childMutationRate = 0.05;
    double parentMutationRate = 0.5;
    int N = 700;
    int maxIter = 200;
    int k = 7;
    bool purgeOp = true;

    std::ofstream file("results.csv");
    file << "Fixed,WinRate,AvgIter\n";

    #pragma omp parallel for schedule(dynamic)
    for (int i=0; i<81; ++i) {

        double winRate = 0.0;
        double avgIter = 0.0;
        int totalWin = 0;

        for (int j=0; j<numSim; ++j) {
            Sudoku sudoku(i);
            int iter = sudoku_GA(sudoku, childMutationRate, parentMutationRate, N, maxIter, k, purgeOp);
            if (iter < maxIter) {
                winRate += 1.0 / numSim;
                avgIter += iter;
                totalWin++;
            }
        }

        if (totalWin > 0) {
            avgIter /= totalWin;
        } else {
            avgIter = maxIter;
        }

        std::cout << "Win rate for " << i << " fixed values: " << winRate << ", Average Iterations Required: " << avgIter << std::endl;
        #pragma omp critical
        {
        file << i << "," << winRate << "," << avgIter << "\n";
        }
    }

    file.close();
    return 0;
}
*/

