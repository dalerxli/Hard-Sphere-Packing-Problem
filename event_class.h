#ifndef event_class
#define event_class

#include <vector>

using namespace std;

namespace kspace {
	class cEvent {
		double time; //the at which the event occurs
		int sphere_count; //the number of spheres involved
		int **ids; //the ids of the spheres involved (the first id is for the primary sphere)
		int event_counts[27]; //the number of events the associated sphere had had at the time of event creation
	public:
		
		cEvent() {
			time = 0; //set the time to zero
			sphere_count = 0; //set the number of spheres involved to zero

			//create a 2d array to store the ids in (max number of spheres involved in an event is 27)
			ids = new int*[27];
			for(int a = 0; a < 27; a++) {
				ids[a] = new int[3];
			}
		}

		~cEvent() {
			for(int a = 0; a < 27; a++) {
				delete [] ids[a];
			}

			delete [] ids;
		}

		//add methods
		void set_time(const double t) {
			time = t;
		}

		void add_sphere(const int i, const int j, const int k, const int event_count) {
			event_counts[sphere_count] = event_count;

			ids[sphere_count][0] = i;
			ids[sphere_count][1] = j;
			ids[sphere_count][2] = k;

			sphere_count++;
		}

		//get methods
		double get_time() const {
			return time;
		}

		int get_sphere_count() const {
			return sphere_count;
		}

		int* get_id(const int index) {
			return ids[index];
		}

		int get_event_count(const int index) const {
			return event_counts[index];
		}
	};

	class compare_events {
	public:
		bool operator()(cEvent *event1, cEvent *event2) {
			return (event1->get_time() > event2->get_time());
		}
	};
}

#endif