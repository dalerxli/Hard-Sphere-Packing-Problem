#include <cstdlib>
#include <queue>
#include <vector>
#include <ctime>
#include <iostream>
#include <fstream>
#include "sim_class.h"

using namespace std;
using namespace kspace;

int main() {
	int size = 4;
	double time_step_size = 10;

	cout << "Begin" << endl;
	sim simulation(100, time_step_size, 0.95, 1e-3,
					size, size, size, 1, 1, 1);
	cout << "Running Simulation" << endl;
	simulation.run();
	cout << "Simulation Complete" << endl;
	

	return 1;
}