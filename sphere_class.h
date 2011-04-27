#ifndef sphere_class
#define sphere_class

#include <cstdlib>
#include <ctime>

using namespace std;

namespace kspace {
	class sphere {
	private:
		int id[3];
		double T;
		int events;
		double pos[3];
		double vel[3];
		
		int lower_clim[3];
		int upper_clim[3];

		double lower_lim[3];
		double upper_lim[3];

	public:
		sphere() {
			T = 0;
			events = 0;

			for(int a = 0; a < 3; a++) {
				pos[a] = 0;
				vel[a] = 0;

				lower_clim[a] = 0;
				upper_clim[a] = 0;

				lower_lim[a] = 0;
				upper_lim[a] = 0;
			}
		}

		//set methods
		void set_id(const int i, const int j, const int k) {
			id[0] = i;
			id[1] = j;
			id[2] = k;
		}

		void inc_events() {
			events++;
		}

		void set_time(const double t) {
			T = t;
		}

		void set_pos(const double x, const double y, const double z) {
			pos[0] = x;
			pos[1] = y;
			pos[2] = z;
		}

		void set_pos(const int index, const double p) {
			pos[index] = p;
		}

		void set_vel(const double x, const double y, const double z) {
			vel[0] = x;
			vel[1] = y;
			vel[2] = z;
		}

		void set_vel(const int index, const double v) {
			vel[index] = v;
		}

		void set_lower_clim(const int i, const int j, const int k) {
			lower_clim[0] = i;
			lower_clim[1] = j;
			lower_clim[2] = k;
		}

		void set_upper_clim(const int i, const int j, const int k) {
			upper_clim[0] = i;
			upper_clim[1] = j;
			upper_clim[2] = k;
		}

		void set_lower_lim(const double i, const double j, const double k) {
			lower_lim[0] = i;
			lower_lim[1] = j;
			lower_lim[2] = k;
		}

		void set_upper_lim(const double i, const double j, const double k) {
			upper_lim[0] = i;
			upper_lim[1] = j;
			upper_lim[2] = k;
		}

		//get methods
		int* get_id() {
			return id;
		}

		double get_time() const {
			return T;
		}

		int get_event_count() const {
			return events;
		}

		double get_pos(const int index) const {
			return pos[index];
		}
		

		double get_vel(const int index) const {
			return vel[index];
		}

		int get_lower_clim(const int index) const {
			return lower_clim[index];
		}

		int get_upper_clim(const int index) const {
			return upper_clim[index];
		}

		double get_lower_lim(const int index) const {
			return lower_lim[index];
		}

		double get_upper_lim(const int index) const {
			return upper_lim[index];
		}

		//validation
		bool within_limits(double *pos) {
			for(int a = 0; a < 3; a++) {
				if(pos[a] < lower_lim[a] || pos[a] > upper_lim[a]) {
					return false;
				}
			}
			return true;
		}

		//print
		void print() {
			cout.precision(3);
			cout << "[";
			for(int a = 0; a < 3; a++) {
				cout << id[a];
				if(a < 2) {
					cout << ", ";
				}
			}
			cout << "] || T: " << get_time() << " | EC: " << get_event_count() << endl;
			cout << "\tP: [";
			for(int a = 0; a < 3; a++) {
				cout << get_pos(a);
				if(a < 2) {
					cout << ", ";
				}
			}
			cout << "] | V: [";
			for(int a = 0; a < 3; a++) {
				cout << get_vel(a);
				if(a < 2) {
					cout << ", ";
				}
			}
			cout << "]" << endl;
			cout << endl;
		}

	};
}

#endif