#ifndef MOGAP_H
#define MOGAP_H

#include"chromosome.h"

typedef std::vector<std::vector<int>> Solution;

void show_solution(Solution& solution);
Solution chromosome_to_solution(Chromosome& chromosome);

struct Aircraft {
    int id;
    double planned_entering_time;  // Pi
    double ground_time;           // Gi
    bool is_arrival;              // true if arrival, false if departure
    double scheduled_time;        // Ai for arrivals, Di for departures
};

struct Gate {
    int id;
    std::vector<Aircraft*> queue;  // Qg
};

class AirportGAP {
public:
    int NA;  // Number of aircraft
    int NG;   // Number of gates
    std::vector<std::vector<double>> Mp;     // Passenger transfer matrix
    std::vector<std::vector<double>> Mpwd;   // Passenger walking distance matrix
    std::vector<std::vector<double>> Mbtd;   // Baggage transport distance matrix
    std::vector<Aircraft> aircraft;
    std::vector<Gate> gates;
    
    // Parameters for the objective function
    double alpha = 0.4;  // weight for TPWD
    double beta = 0.4;   // weight for TBTD
    double phi = 25.0;   // system parameter for waiting time

    AirportGAP(int num_aircraft, int num_gates);
    void generateAircraftSet(double start_time, double end_time);
    int calculateFitness(Solution& solution);

private:
    double calculateTPWD(const Solution& solution);
    double calculateTBTD(const Solution& solution);
    double calculateTPWT(const Solution& solution);
    std::vector<double> calculateWaitingTimes(const Solution& solution);
};


#endif // MOGAP_H