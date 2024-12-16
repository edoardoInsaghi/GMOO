#include"mogap.h"

void show_solution(Solution& solution) {
    for(int i = 0; i < solution.size(); i++) {
        std::cout << "Gate " << i << ": ";
        for(int j = 0; j < solution[i].size(); j++) {
            std::cout << solution[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

Solution chromosome_to_solution(Chromosome& chromosome) {
    Solution solution(chromosome.NG + 1, std::vector<int>());
    for(int i : chromosome.used_gates) {
        int first = find_first_in_queue_g(chromosome, i);
        if (first == -1) {
            throw std::runtime_error("first of queue is -1 but gate is used in chromosome");
        }
        int gate = chromosome.gate_of_plane(first);
        if (gate != i) {
            throw std::runtime_error("gate of first plane in queue is not the same as the gate wtf");
        }
        while (first != -1) {
            solution[gate].push_back(first);
            first = successor_in_queue(chromosome, first);
        }
    }
    solution[chromosome.NG].push_back(chromosome.NA); // to the dummy gate we always assign the dummy plane
    return solution;
}

AirportGAP::AirportGAP(int num_aircraft, int num_gates) : NA(num_aircraft), NG(num_gates) {
    Mp = std::vector<std::vector<double>>(NA+1, std::vector<double>(NA+1, 0.0));
    Mpwd = std::vector<std::vector<double>>(NG+1, std::vector<double>(NG+1, 0.0));
    Mbtd = std::vector<std::vector<double>>(NG+1, std::vector<double>(NG+1, 0.0));
    
    gates.resize(NG + 1);  // +1 for dummy gate
    for(int i=0; i<=NG; i++) {
        gates[i].id = i;
    }
}

void AirportGAP::generateAircraftSet(double start_time, double end_time) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> time_dist(start_time, end_time);
    std::uniform_real_distribution<> ground_time_dist(30, 120);  // Ground time between 30-120 minutes
    
    aircraft.clear();
    for(int i=0; i<=NA; i++) {
        Aircraft ac;
        ac.id = i;
        ac.is_arrival = (i < NA/2);  // First half arrivals, second half departures
        
        if(ac.is_arrival) {
            ac.scheduled_time = time_dist(gen);
            ac.ground_time = ground_time_dist(gen);
            ac.planned_entering_time = ac.scheduled_time;  // Pi = Ai for arrivals
        } else {
            ac.ground_time = ground_time_dist(gen);
            ac.scheduled_time = time_dist(gen);
            ac.planned_entering_time = ac.scheduled_time - ac.ground_time;  // Pi = Di - Gi for departures
        }
        
        aircraft.push_back(ac);
    }
    
    // Generate random passenger transfer matrix
    std::uniform_int_distribution<> pass_dist(0, 200);
    for(int i=0; i<=NA; i++) {
        for(int j=0; j<=NA; j++) {
            if(i != j) {  // No transfers to same aircraft
                Mp[i][j] = pass_dist(gen);
            }
        }
    }
    
    // Generate distance matrices
    std::uniform_int_distribution<> dist_dist(50, 500);  // distances between 50-500 meters
    for(int i=0; i<=NG; i++) {
        for(int j=0; j<=NG; j++) {
            if(i != j) {
                Mpwd[i][j] = dist_dist(gen);
                Mbtd[i][j] = dist_dist(gen);
            }
        }
    }
    Mbtd[NG][NG] = 0;
}


int AirportGAP::calculateFitness(Solution& solution) {
    double TPWD = calculateTPWD(solution);
    double TBTD = calculateTBTD(solution);
    double TPWT = calculateTPWT(solution);
    
    return -(alpha * TPWD + beta * TBTD + (1 - alpha - beta) * phi * TPWT);
}

double AirportGAP::calculateTPWD(const Solution& solution) {
    double total = 0;
    for(int g=0; g<=NG; g++) {
        for(int j=0; j<solution[g].size(); j++) {
            for(int i=0; i<NA; i++) {
                int aircraft_id = solution[g][j];
                if(aircraft_id >= 0 && aircraft_id<NA) {
                    total += Mp[aircraft_id][i] * Mpwd[g][gates[aircraft[i].id].id];
                }
            }
        }
    }
    return total;
}

double AirportGAP::calculateTBTD(const Solution& solution) {
    double total = 0;
    for(int g=0; g<=NG; g++) {
        for(int j=0; j<solution[g].size(); j++) {
            for(int i=0; i<NA; i++) {
                int aircraft_id = solution[g][j];
                if(aircraft_id >= 0 && aircraft_id < NA) {
                    total += Mp[aircraft_id][i] * Mbtd[g][gates[aircraft[i].id].id];
                }
            }
        }
    }
    return total;
}

double AirportGAP::calculateTPWT(const Solution& solution) {
    double total = 0;
    std::vector<double> waiting_times = calculateWaitingTimes(solution);
    for(int i=0; i<=NA; i++) {
        double Wi = waiting_times[i];
        for(int j=0; j<NA; j++) {
            total += Wi * (Mp[i][j] + Mp[j][i]);
        }
    }
    return total;
}

std::vector<double> AirportGAP::calculateWaitingTimes(const Solution& solution) {
    std::vector<double> waiting_times(NA);
    for(int g=0; g<NG; g++) {
        for(int j=0; j<solution[g].size(); j++) {
            int aircraft_id = solution[g][j];
            if(aircraft_id >= 0 && aircraft_id < NA) {
                double Ei;
                if(j == 0) {
                    Ei = aircraft[aircraft_id].planned_entering_time;
                } else {
                    int prev_aircraft = solution[g][j-1];
                    Ei = std::max(aircraft[aircraft_id].planned_entering_time,
                                waiting_times[prev_aircraft] + aircraft[prev_aircraft].ground_time);
                }
                waiting_times[aircraft_id] = Ei - aircraft[aircraft_id].planned_entering_time;
            }
        }
    }
    return waiting_times;
}