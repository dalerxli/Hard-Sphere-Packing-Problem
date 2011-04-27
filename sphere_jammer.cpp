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
	cout << "Begin" << endl;
	sim simulation(100e3, 0.95, 1e-6,
					3, 3, 3, 1, 1, 1);
	cout << "Running Simulation" << endl;
	simulation.run();
	cout << "Simulation Complete" << endl;
	

	return 1;
}