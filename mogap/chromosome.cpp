#include "chromosome.h"


///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// Chromosome /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

Chromosome::Chromosome(int NA, int NG, std::vector<int>& relative_position_genes, std::vector<int>& gate_genes)
    : NA(NA), NG(NG), relative_position_genes(relative_position_genes), gate_genes(gate_genes) {
    for (int& g : gate_genes) {
        if (g != -1) {
            used_gates.insert(g);
        }
    }
}

Chromosome::Chromosome(const Chromosome& other)
    : NA(other.NA),
      NG(other.NG),
      relative_position_genes(other.relative_position_genes),
      gate_genes(other.gate_genes),
      used_gates(other.used_gates) {}

Chromosome& Chromosome::operator=(const Chromosome& other) {
    if (this == &other) {
        return *this;
    }
    NA = other.NA;
    NG = other.NG;
    relative_position_genes = other.relative_position_genes;
    gate_genes = other.gate_genes;
    used_gates = other.used_gates;
    return *this;
}

bool Chromosome::operator==(const Chromosome& other) const {
    return NA == other.NA && NG == other.NG && relative_position_genes == other.relative_position_genes && gate_genes == other.gate_genes;
}

void Chromosome::set_used_gates() {
    used_gates.clear();
    for (int& g : gate_genes) {
        if (g != -1) {
            used_gates.insert(g);
        }
    }
}


Chromosome init_random_chromosome(int NA, int NG) {

    std::vector<int> relative_position_genes(NA*NA, 0);
    std::vector<int> gate_genes(NA, -1);
    Chromosome chromosome(NA, NG, relative_position_genes, gate_genes);

    std::set<int> plane_set;
    for (int i=0; i<NA; i++) {
        plane_set.insert(i);
    }

    std::set<int> gate_set;
    for (int i=0; i<NG; i++) {
        gate_set.insert(i);
    }

    int diagonal_places = 0;
    while (!plane_set.empty()) {
        if (diagonal_places == NG) { // One or more planes on the diagonal would break condition 4

            int first_in_queue = select_random_first_in_queue(chromosome);
            int g = chromosome.gate_genes[first_in_queue];

            while (!no_nonzero_row(chromosome, first_in_queue)) {
                first_in_queue = find_no_nonzero_row(chromosome, first_in_queue);
            }
            int j = select_random_from_set_and_remove(plane_set);
            chromosome.relative_position_genes[first_in_queue*NA + j] = 1;
            chromosome.gate_genes[j] = g;
            chromosome.used_gates.insert(g);
        }
        else {
            int random_plane = select_random_from_set_and_remove(plane_set);
            int random_gate = select_random_from_set_and_remove(gate_set);

            chromosome.relative_position_genes[random_plane*NA + random_plane] = 1;
            diagonal_places++;
            chromosome.gate_genes[random_plane] = random_gate;
            chromosome.used_gates.insert(random_gate);
        }
    }

    return chromosome;
}





///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// Mutations //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

void apply_mutation_1(Chromosome& chromosome) { // shift position of two successive aircraft in the same queue

    // (a ->) b -> c (-> d)
    // (a ->) c -> b (-> d) 

    int mutation_idx = select_random_nonzero_gene(chromosome);
    int c = mutation_idx / chromosome.NA; 
    int b = mutation_idx % chromosome.NA;

    int a = successor_in_queue(chromosome, b);
    bool has_successor = a != -1;

    int d = predecessor_in_queue(chromosome, c);
    bool has_predecessor = d != -1;

    std::cout << "Mutation 1: " << a << " " << b << " " << c << " " << d << std::endl;
    std::cout << "After mutation 1:" << a << " " << c << " " << b << " " << d << std::endl;

    bool is_first_in_queue = chromosome(c, c) == 1;

    if (has_predecessor != !is_first_in_queue) {
        throw std::runtime_error("Predecessor and first in queue are not consistent\n");
    }

    if (b==c) {
        return;
    }

    if (chromosome.gate_of_plane(b) != chromosome.gate_of_plane(c)) {
        throw std::runtime_error("Planes are not in the same queue wtf\n");      
    }

    int case_id = (is_first_in_queue << 1) | has_successor;
    std::cout << "Mutation 1, case id: " << case_id << std::endl;
    switch (case_id) {
        case 0: // !is_first_in_queue && !has_successor
            chromosome.unset(d, c);
            chromosome.set(d, b);
            chromosome.unset(c, b);
            chromosome.set(b, c);
            break;
        
        case 1: // !is_first_in_queue && has_successor
            chromosome.unset(d, c);
            chromosome.set(d, b);
            chromosome.unset(c, b);
            chromosome.set(b, c);
            chromosome.unset(b, a);
            chromosome.set(c, a);
            break;

        case 2: // is_first_in_queue && !has_successor
            chromosome.unset(c, c);
            chromosome.set(b, c);
            chromosome.unset(c, b);
            chromosome.set(b, b);
            break;

        case 3: // is_first_in_queue && has_successor
            chromosome.unset(c, c);
            chromosome.set(b, b);
            chromosome.unset(c, b);
            chromosome.set(b, c);
            chromosome.unset(b, a);
            chromosome.set(c, a);
            break;
    }
}


void apply_mutation_2(Chromosome& chromosome) {

    // (a1 ->) b1 -> c1
    // (a2 ->) b2 -> c2

    // (a1 ->) b2 -> c1
    // (a2 ->) b1 -> c2

    std::array<int, 2> mutation_indices = select_2_random_nonzero_genes(chromosome);
    int idx1 = mutation_indices[0];
    int idx2 = mutation_indices[1];

    int c1 = idx1 / chromosome.NA;
    int b1 = idx1 % chromosome.NA;
    int c2 = idx2 / chromosome.NA;
    int b2 = idx2 % chromosome.NA;

    int g1 = chromosome.gate_of_plane(c1);
    int g2 = chromosome.gate_of_plane(c2);
    std::cout << "Mutation 2, gates: " << g1 << " " << g2 << std::endl;
    if (g1 == g2) { return; }

    chromosome.gate_genes[b1] = g2;
    chromosome.gate_genes[b2] = g1;

    int a1 = successor_in_queue(chromosome, b1);
    int a2 = successor_in_queue(chromosome, b2);

    std::cout << "Mutation 2: " << a1 << " " << b1 << " " << c1 << " | " << a2 << " " << b2 << " " << c2 << std::endl;
    std::cout << "After mutation 2:" << a1 << " " << b2 << " " << c1 << " | " << a2 << " " << b1 << " " << c2 << std::endl;

    bool has_successor1 = a1 != -1;
    bool has_successor2 = a2 != -1;
    int has_successor = (has_successor1 << 1) | has_successor2;

    bool is_first_in_queue1 = b1 == c1;
    bool is_first_in_queue2 = b2 == c2;
    int is_first_in_queue = (is_first_in_queue1 << 1) | is_first_in_queue2;

    std::cout << "Mutation 2, case id: " << is_first_in_queue << std::endl;
    std::cout << "Mutation 2, has successor: " << has_successor << std::endl;

    switch (is_first_in_queue) {
        case 0: // !is_first_in_queue1 && !is_first_in_queue2
            chromosome.unset(c1, b1);
            chromosome.set(c1, b2);
            chromosome.unset(c2, b2);
            chromosome.set(c2, b1);

            switch (has_successor) {
                case 0: // !has_successor1 && !has_successor2
                    break;
                case 1: // !has_successor1 && has_successor2
                    chromosome.unset(b2, a2);
                    chromosome.set(b1, a2);
                    break;
                case 2: // has_successor1 && !has_successor2
                    chromosome.unset(b1, a1);
                    chromosome.set(b2, a1);
                    break;
                case 3: // has_successor1 && has_successor2
                    chromosome.unset(b1, a1);
                    chromosome.set(b2, a1);
                    chromosome.unset(b2, a2);
                    chromosome.set(b1, a2);
                    break;
            }
            break;
        
        case 1: // !is_first_in_queue1 && is_first_in_queue2
            chromosome.unset(c1, b1);
            chromosome.set(b1, b1);
            chromosome.unset(c2, c2);
            chromosome.set(c1, b2);

            switch (has_successor) {
                case 0: // !has_successor1 && !has_successor2
                    break;
                case 1: // !has_successor1 && has_successor2
                    chromosome.unset(c2, a2);
                    chromosome.set(b1, a2);
                    break;
                case 2: // has_successor1 && !has_successor2
                    chromosome.unset(b1, a1);
                    chromosome.set(c2, a1);
                    break;
                case 3: // has_successor1 && has_successor2
                    chromosome.unset(b1, a1);
                    chromosome.set(c2, a1);
                    chromosome.unset(c2, a2);
                    chromosome.set(b1, a2);
                    break;
            }
            break;

        case 2: // is_first_in_queue1 && !is_first_in_queue2
            chromosome.unset(c2, b2);
            chromosome.set(b2, b2);
            chromosome.unset(b1, b1);
            chromosome.set(c2, b1);

            switch (has_successor) {
                case 0: // !has_successor1 && !has_successor2
                    break;
                case 1: // !has_successor1 && has_successor2
                    chromosome.unset(b2, a2);
                    chromosome.set(c1, a2);
                    break;
                case 2: // has_successor1 && !has_successor2
                    chromosome.unset(c1, a1);
                    chromosome.set(b2, a1);
                    break;
                case 3: // has_successor1 && has_successor2
                    chromosome.unset(b2, a2);
                    chromosome.set(c1, a2);
                    chromosome.unset(c1, a1);
                    chromosome.set(b2, a1);
                    break;
            }
            break;

        case 3: // is_first_in_queue1 && is_first_in_queue2
            switch (has_successor) {
                case 0: // !has_successor1 && !has_successor2
                    break;
                case 1: // !has_successor1 && has_successor2
                    chromosome.unset(b2, a2);
                    chromosome.set(b1, a2);
                    break;
                case 2: // has_successor1 && !has_successor2
                    chromosome.unset(b1, a1);
                    chromosome.set(b2, a1);
                    break;
                case 3: // has_successor1 && has_successor2
                    chromosome.unset(b1, a1);
                    chromosome.set(b2, a1);
                    chromosome.unset(b2, a2);
                    chromosome.set(b1, a2);
                    break;
            }
            break;
    }
}





///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// Crossover //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

Chromosome perform_crossover(Chromosome& parent1, Chromosome& parent2) {

    if (parent1 == parent2) { return std::move(parent1); }

    int NA = parent1.NA;
    int NG = parent1.NG;
    std::vector<int> relative_position_genes(NA*NA, 0);
    std::vector<int> gate_genes(NA, -1);
    Chromosome child(NA, NG, relative_position_genes, gate_genes);

    for (int i=0; i<NA; i++) {
        for (int j=0; j<NA; j++) {
            child.relative_position_genes[i*NA + j] = parent1(i, j) && parent2(i, j);
        }
        child.gate_genes[i] = parent1.gate_of_plane(i) == parent2.gate_of_plane(i) ? parent1.gate_of_plane(i) : -1;
        child.used_gates.insert(gate_genes[i]);
    }

    for (int i=0; i<NA; i++) {

        int rand = std::rand() % 2;
        if (rand == 0) {
            int gate = parent1.gate_of_plane(i);
            child.gate_genes[i] = gate;
            child.used_gates.insert(gate);
            int num_first = gate_num_first_in_queue(child, gate);
            bool free_column = is_column_free(child, i);
            if (num_first > 1) {
                child.unset(i, i);
            }
            else if (num_first < 1 && free_column) {
                child.set(i, i);
            }
            else if (num_first < 1 && !free_column) {
                int eliminate = find_nonzero_in_column(child, i);
                child.unset(eliminate, i);
                child.set(i, i);
            }
    
        }
        else {
            int gate = parent2.gate_of_plane(i);
            child.gate_genes[i] = gate;
            child.used_gates.insert(gate);
            int num_first = gate_num_first_in_queue(child, gate);
            bool free_column = is_column_free(child, i);
            if (num_first > 1) {
                child.relative_position_genes[i*NA + i] = 0;
            }
            else if (num_first < 1 && free_column) {
                child.relative_position_genes[i*NA + i] = 1;
            }    
            else if (num_first < 1 && !free_column) {
                int eliminate = find_nonzero_in_column(child, i);
                child.unset(eliminate, i);
                child.set(i, i);
            }
        }
    }

    Chromosome dummy = child;

    for(int i=0; i<NA; i++) {
        dummy.relative_position_genes[i*NA + i] = -1;
        for (int j=0; j<NA; j++) {
            if (child(i, j) == 1 && i != j) {
                for (int m=0; m<NA; m++) {
                    dummy.relative_position_genes[m*NA + j] = -1;
                    dummy.relative_position_genes[i*NA + m] = -1;
                }
            }
        }
        if (child(i, i) == 1) {
            for (int j=0; j<NA; j++) {
                dummy.relative_position_genes[j*NA + i] = -1;
            }
        }
    }

    child.used_gates.erase(-1);
    bool full = constraint_1(child);
    while (!full) {

        int plane = select_random_non_assigned_plane(child);
        int g = child.gate_of_plane(plane);

        int predecessor1 = predecessor_in_queue(parent1, plane);
        int has_predecessor1 = predecessor1 != -1 && dummy(predecessor1, plane) == 0;
        
        int predecessor2 = predecessor_in_queue(parent2, plane);
        int has_predecessor2 = predecessor2 != -1 && dummy(predecessor2, plane) == 0;

        int has_predecessor = (has_predecessor1 << 1) | has_predecessor2;

        int i3;
        switch (has_predecessor) {
            case 0:  // !has_predecessor1 && !has_predecessor2
                i3 = find_last_in_queue(child, g);
                break;

            case 1: // !has_predecessor1 && has_predecessor2
                i3 = predecessor2;
                break;

            case 2: // has_predecessor1 && !has_predecessor2   
                i3 = predecessor1;
                break;

            case 3: // has_predecessor1 && has_predecessor2
                int rand = std::rand() % 2;
                i3 = rand == 0 ? predecessor1 : predecessor2;
                break;
        }

        child.set(i3, plane);
        dummy.set(i3, plane);

        for (int m=0; m<NA; m++) {
            dummy.relative_position_genes[m*NA + plane] = -1;
            dummy.relative_position_genes[i3*NA + m] = -1;
        }

        full = constraint_1(child);
    }

    return child;
}





///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Constraints /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

bool constraint_1(Chromosome& chromosome) {
    int sum = std::accumulate(chromosome.relative_position_genes.begin(), chromosome.relative_position_genes.end(), 0);
    return sum == chromosome.NA;
}

bool constraint_2(Chromosome& chromosome) {
    for (int i=0; i<chromosome.NA; i++) {
        int rowsum = 0;
        int ii = chromosome(i, i);
        for (int j=0; j<chromosome.NA; j++) {
            rowsum += chromosome(i, j);
        }
        if (ii > 0) {
            if (rowsum > 2) { return false; }
        }
        else if (ii == 0) {
            if (rowsum > 1) { return false; }
        }
    }
    return true;
}

bool constraint_3(Chromosome& chromosome) {
    for (int j=0; j<chromosome.NA; j++) {
        int colsum = 0;
        for (int i=0; i<chromosome.NA; i++) {
            colsum += chromosome.relative_position_genes[i*chromosome.NA + j];
        }
        if (colsum != 1) { return false; }
    }
    return true;
}

bool constraint_4(Chromosome& chromosome) {
    int diagsum = 0;
    for (int i=0; i<chromosome.NA; i++) {
        diagsum += chromosome(i, i);
    }
    if (diagsum < 1 || diagsum > chromosome.NG) { return false; } 
    return true;
};

bool constraint_5(Chromosome& chromosome) {
    for (std::set<int>::iterator it=chromosome.used_gates.begin(); it!=chromosome.used_gates.end(); it++) {
        int g = *it;
        int has_first = 0;
        for (int j=0; j<chromosome.NA; j++) {
            if (chromosome.gate_of_plane(j) == g) {  // If gate is used by the plane j
                has_first += chromosome(j, j);
            }
        }
        if (has_first != 1) { return false; }
    }
    return true;
}

bool queue_gate_consistency(Chromosome& chromosome) {
    for (int i=0; i<chromosome.NA; i++) {
        int g = chromosome.gate_of_plane(i);
        if (g == -1) { continue; }
        int first_in_queue = find_first_in_queue_g(chromosome, g);
        if (chromosome.gate_of_plane(first_in_queue) != g) {
            std::cout << "Gate of first plane in queue is not the same as the gate" << std::endl;
            return false;
            
        }
        if (first_in_queue == -1) {
            std::cout << "First in queue is -1 but gate is used in chromosome" << std::endl;
            return false;
        }
        while (first_in_queue != -1) {
            first_in_queue = successor_in_queue(chromosome, first_in_queue);
            if (first_in_queue!=-1 && chromosome.gate_of_plane(first_in_queue) != g) {
                std::cout << "Plane: " << first_in_queue << " found in queue: " << chromosome.gate_of_plane(first_in_queue) << " rather than: " << g << std::endl;
                return false;
            }
        }
    }
    return true;
}

bool is_valid_chromosome(Chromosome& chromosome) {
    bool c1 = constraint_1(chromosome);
    if (!c1) { std::cout << "Constraint 1 failed\n"; }
    bool c2 = constraint_2(chromosome);
    if (!c2) { std::cout << "Constraint 2 failed\n"; }
    bool c3 = constraint_3(chromosome);
    if (!c3) { std::cout << "Constraint 3 failed\n"; }
    bool c4 = constraint_4(chromosome);
    if (!c4) { std::cout << "Constraint 4 failed\n"; }
    bool c5 = constraint_5(chromosome);
    if (!c5) { std::cout << "Constraint 5 failed\n"; }
    bool c6 = queue_gate_consistency(chromosome);
    if (!c6) { std::cout << "Queue gate consistency failed\n"; }

    if (!c1 || !c2 || !c3 || !c4 || !c5 || !c6) {
        show_chromosome(chromosome);
        show_chromosome_absolute_position(chromosome);
        //throw std::runtime_error("Chromosome is invalid\n");
    }
    return c1 && c2 && c3 && c4 && c5;
}





///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// Utilities //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

void show_chromosome_absolute_position(Chromosome& chromosome) {
    std::vector<int> first_in_queue = find_first_in_queue(chromosome);
    for (int i=0; i<chromosome.NG; i++) {
        int plane = first_in_queue[i];
        std::cout << "Gate " << i << ": ";
        while (plane != -1) {
            std::cout << plane << " <-- ";
            plane = successor_in_queue(chromosome, plane);
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void show_chromosome(Chromosome& chromosome) {
    for (int i=0; i<chromosome.NA; i++) {
        for (int j=0; j<chromosome.NA; j++) {
            std::cout << std::setw(COLWIDTH) << chromosome(i, j) << " ";
        }
        std::cout << std::endl;
    }
    for (int i=0; i<chromosome.NA; i++) {
        std::cout << std::setw(COLWIDTH) << chromosome.gate_of_plane(i) << " ";
    }
    std::cout << "\n" << std::endl;
}

int find_no_nonzero_row(Chromosome& chromosome, int i) {
    for (int j=0; j<chromosome.NA; j++) {
        if (j == i) { continue; }
        if (chromosome(i, j) == 1) {
            return j;
        }
    }
    throw std::runtime_error("Row is all zeros.\n");
}

bool no_nonzero_row(Chromosome& chromosome, int i) {
    for (int j=0; j<chromosome.NA; j++) {
        if (j == i) { continue; }
        if (chromosome(i, j) != 0) {
            return false;
        }
    }
    return true;
}

int select_random_non_assigned_plane(Chromosome& chromosome) {
    std::set<int> unassigned_planes;
    for (int i=0; i<chromosome.NA; i++) {
        bool assigned = false;
        for (int j=0; j<chromosome.NA; j++) {
            if (chromosome(j, i) != 0) {
                assigned = true;
                break;
            }
        }
        if (!assigned) {
            unassigned_planes.insert(i);
        }
    }
    if (unassigned_planes.empty()) {
        throw std::runtime_error("All planes are assigned\n");
    }
    return select_random_from_set(unassigned_planes);
}

int select_random_first_in_queue(Chromosome& chromosome) {
    std::vector<int> first_in_queue;
    for (int i=0; i<chromosome.NA; i++) {
        if (chromosome(i, i) == 1) {
            first_in_queue.push_back(i);
        }
    }

    if (first_in_queue.empty()) {
        throw std::runtime_error("Diagonal of relative position matrix equals [0, 0, ..., 0]\n");
    }
    
    int random_index = rand() % first_in_queue.size();
    return first_in_queue[random_index];
}

bool is_column_free(Chromosome& chromosome, int j) {
    for (int i=0; i<chromosome.NA; i++) {
        if (chromosome(i, j) != 0) {
            return false;
        }
    }
    return true;
}

int select_random_nonzero_gene(Chromosome& chromosome) {
    std::vector<int> nonzero_genes;
    for (int i=0; i<chromosome.NA; i++) {
        for (int j=0; j<chromosome.NA; j++) {
            if (chromosome(i, j) != 0) {
                nonzero_genes.push_back(i*chromosome.NA + j);
            }
        }
    }
    int random_index = rand() % nonzero_genes.size();
    return nonzero_genes[random_index];
}

std::array<int, 2> select_2_random_nonzero_genes(Chromosome& chromosome) {
    std::vector<int> nonzero_genes;
    for (int i=0; i<chromosome.NA; i++) {
        for (int j=0; j<chromosome.NA; j++) {
            if (chromosome(i, j) != 0) {
                nonzero_genes.push_back(i*chromosome.NA + j);
            }
        }
    }
    int random_index_1 = rand() % nonzero_genes.size();
    int random_index_2 = rand() % nonzero_genes.size();
    while (random_index_2 == random_index_1) {
        random_index_2 = rand() % nonzero_genes.size();
    }
    return {nonzero_genes[random_index_1], nonzero_genes[random_index_2]};
}

int predecessor_in_queue(Chromosome& chromosome, int i) {
    for (int j=0; j<chromosome.NA; j++) {
        if (j == i) { continue; }
        if (chromosome(j, i) == 1) {
            return j;
        }
    }
    return -1;
}

int successor_in_queue(Chromosome& chromosome, int i) {
    for (int j=0; j<chromosome.NA; j++) {
        if (j == i) { continue; }
        if (chromosome(i, j) == 1) {
            return j;
        }
    }
    return -1;
}

std::vector<int> find_first_in_queue(Chromosome& chromosome) {
    std::vector<int> first_in_queue(chromosome.NG, -1);
    for (int i=0; i<chromosome.NA; i++) {
        if (chromosome(i, i) == 1) {
            first_in_queue[chromosome.gate_of_plane(i)] = i;
        }
    }
    return first_in_queue;
}

int find_first_in_queue_g(Chromosome& chromosome, int g) {
    for (int i=0; i<chromosome.NA; i++) {
        if (chromosome.gate_of_plane(i) == g && chromosome(i, i) == 1) {
            return i;
        }
    }
    return -1;
}

int find_last_in_queue(Chromosome& chromosome, int i) {
    int last = find_first_in_queue_g(chromosome, i);
    while (true) {
        int next = successor_in_queue(chromosome, last);
        if (next == -1) {
            return last;
        }
        last = next;
    }
}

int gate_num_first_in_queue(Chromosome& chromosome, int gate) {
    int count = 0;
    for (int i=0; i<chromosome.NA; i++) {
        if (chromosome.gate_of_plane(i) == gate && chromosome(i, i) == 1) {
            count++;
        }
    }
    return count;
}

int find_nonzero_in_column(Chromosome& chromosome, int j) {
    for (int i=0; i<chromosome.NA; i++) {
        if (chromosome(i, j) != 0) {
            return i;
        }
    }
    return -1;
}