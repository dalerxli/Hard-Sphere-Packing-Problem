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
	int size = 10;
	double time_step_size = 5;

	cout << "Begin" << endl;
	sim simulation(210e3, time_step_size, 0.99, 1e-6,
					size, size, size, 1, 1, 1);
	cout << "Running Simulation" << endl;
	simulation.run();
	cout << "Simulation Complete" << endl;
	

	return 1;
}