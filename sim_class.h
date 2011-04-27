#ifndef sim_class
#define sim_class

#include <iostream>
#include <ctime>
#include "sphere_class.h"
#include "event_class.h"

using namespace std;

namespace kspace {
	class sim {
	private:
		double current_time;
		double end_time;

		double D0; //the initial diameter of the spheres
		double gamma; //the growth rate of the radii of the spheres

		int dim[3];
		double cell_dim[3];

		int accepted_counter;
		int discarded_counter;
		int sev_counter;
		int bev_counter;

		vector<sphere> spheres;
		priority_queue<cEvent*, vector<cEvent*>, compare_events> event_q;
	public:
		sim(double end_time, double initial_diameter, double Growth_Rate, int width, int height, int depth, double cell_width, double cell_height, double cell_depth); //initialises the sim
		~sim();

		void run();
		
		sphere& get_sphere(const int i, const int j, const int k); //return the sphere specified
		
		//get methods
		double get_end_time() const;
		double growth_rate() const;
		double diameter(double t) const;

		//calculates all potential future events for the specified sphere and adds them to the queue
		void calculate_events(const int i, const int j, const int k); 
		double calculate_boundary_event(const int i, const int j, const int k);
		double calculate_sphere_event(const int i, const int j, const int k, const int l, const int m, const int n);

		//advances the simulations
		void advance(cEvent *ev);
		void calculate_boundary_velocity(const int i, const int j, const int k);
		void calculate_sphere_velocity(const int i, const int j, const int k, const int l, const int m, const int n);

		void make_concurrent(); //advance ALL spheres to the same time, recalculates events for all spheres
		
		double round(double v); //rounds the double v to an accuracy of 9 decimal places and return the rounded value

		//save method
		void save(char fname[128]);

		//pause
		void pause();
	};
}


using namespace kspace;

sim::sim(double end_time, double initial_diameter, double Growth_Rate, int width, int height, int depth, double cell_width, double cell_height, double cell_depth) {
	current_time = 0;
	this->end_time = end_time;

	accepted_counter = 0;
	discarded_counter = 0;
	sev_counter = 0;
	bev_counter = 0;

	//set the growth rate and the initial diameter
	D0 = initial_diameter;
	gamma = Growth_Rate;

	//empty the event queue if it has some events stored already, probably unecessary 
	if(event_q.size()) {
		event_q.empty();
	}

	dim[0] = width;
	dim[1] = height;
	dim[2] = depth;

	cell_dim[0] = cell_width;
	cell_dim[1] = cell_height;
	cell_dim[2] = cell_depth;

	//create the spheres and initialise them all
	spheres.resize(dim[0]*dim[1]*dim[2]);

	srand((int) time(NULL) + clock() + rand());
	for(int a = 0; a < dim[0]; a++) {
		for(int b = 0; b < dim[1]; b++) {
			for(int c = 0; c < dim[2]; c++) {
				sphere &s = get_sphere(a, b, c);
				s.set_id(a, b, c);
				s.set_pos((a + 0.5)*cell_dim[0], (b + 0.5)*cell_dim[1], (c + 0.5)*cell_dim[2]);
				s.set_vel((double) (rand()%1000) / 1000, (double) (rand()%1000) / 1000, (double) (rand()%1000) / 1000);
				s.set_lower_clim(a, b, c);
				s.set_upper_clim(a+1, b+1, c+1);
				s.set_lower_lim(a*cell_dim[0], b*cell_dim[1], c*cell_dim[2]);
				s.set_upper_lim((a+1)*cell_dim[0], (b+1)*cell_dim[1], (c+1)*cell_dim[2]);
			}
		}
	}
}

sim::~sim() {
	cEvent *ev;
	while(event_q.size()) {
		ev = event_q.top();
		event_q.pop();
		delete ev;
	}
}

void sim::run() {
	int stime = clock();

	cEvent *ev;

	do {
		//pause();
		if(event_q.size() > 0) {
			ev = event_q.top();
			event_q.pop();
			advance(ev);
			delete ev;
		} else {
			cout << "Calculating Events" << endl;
			for(int a = 0; a < dim[0]; a++) {
				for(int b = 0; b < dim[1]; b++) {
					for(int c = 0; c < dim[2]; c++) {
						calculate_events(a, b, c);
					}
				}
			}
			cout << "Calculations Finished" << endl;
		}

	} while(current_time < end_time);

	int ftime = clock();

	cout << "Run Time: " << ftime - stime << " | " << (double) (ftime - stime) / 1000 << endl;
}

sphere& sim::get_sphere(const int i, const int j, const int k) {
	int x(i), y(j), z(k);

	if(i < 0) {x = i + dim[0];}
	else if(i >= dim[0]) {x = i - dim[0];}

	if(j < 0) {y = j + dim[1];}
	else if(j >= dim[1]) {y = j - dim[1];}

	if(k < 0) {z = k + dim[2];}
	else if(k >= dim[2]) {z = k - dim[2];}

	return spheres.at(dim[0]*dim[1]*z + dim[0]*y + x);
}

double sim::get_end_time() const {
	return end_time;
}

double sim::growth_rate() const {
	return gamma;
}

double sim::diameter(double t) const {
	return (D0 + growth_rate()*t);
}

void sim::calculate_events(const int i, const int j, const int k) {
	double sev;
	bool create_new_event;
	vector<cEvent*> events;
	cEvent *ev;
	sphere &s = get_sphere(i, j, k);
	
	double bev = calculate_boundary_event(i, j, k); //calculate the time at which collision with boundary will occur

	//if this is greater than or equal to zero calculate the times at whill collisions with other spheres will occur
	if(bev >= 0) {
		//calculate all potential events 
		for(int a = s.get_lower_clim(0) - 1; a <= s.get_upper_clim(0); a++) {
			for(int b = s.get_lower_clim(1) - 1; b <= s.get_upper_clim(1); b++) {
				for(int c = s.get_lower_clim(2) - 1; c <= s.get_upper_clim(2); c++) {
					if(a != i || b != j || c != k) {
						sev = calculate_sphere_event(i, j, k, a, b ,c); //calculate the time at which a collision with sphere[a, b, c] will occur

						if(sev >= 0 && sev <= bev) {
							//if the time is greater than or equal to zero and less than or equal to the time at which sphere[i, j, k] 
							//collides with its cell boundary then check to see if it occurs at the same time as another sphere-sphere collision
							create_new_event = true;
							for(unsigned int q = 0; q < events.size(); q++) {
								if(round(events.at(q)->get_time() - sev) == 0) {
									//add sphere[a, b, c] to the event
									events.at(q)->add_sphere(a, b, c, get_sphere(a, b, c).get_event_count());
									create_new_event = false;
									break;
								}
							}

							if(create_new_event) {
								//if there are no other events which occur at the same time then create a new one and add it both to the queue
								//and the vector 'events'
								ev = new cEvent();
								ev->set_time(sev);
								ev->add_sphere(i, j, k, s.get_event_count());
								ev->add_sphere(a, b, c, get_sphere(a, b, c).get_event_count());
								events.push_back(ev);
								event_q.push(ev);
								ev = NULL;
							}
						}
					}
				}
			}
		}

		//add the boundary collision to the queue
		ev = new cEvent();
		ev->set_time(bev);
		ev->add_sphere(i, j, k, s.get_event_count());
		event_q.push(ev);
	}
}

double sim::calculate_boundary_event(const int i, const int j, const int k) {
	sphere &s = get_sphere(i, j, k);
	double dN1, dN2, t;
	double min_t(-1);

	for(int a = 0; a < 3; a++) {
		t = -1;
		dN1 = round(s.get_upper_lim(a) - s.get_pos(a));
		dN2 = round(s.get_lower_lim(a) - s.get_pos(a));

		if(s.get_vel(a) > 0 &&  dN1>= 0) {
			//if the velocity of the sphere in this axis is poisitve and the position of the sphere
			//is less than the upper cell boundary position in this axis then calculate the time until the collision
			//with this boundary
			t = dN1 / s.get_vel(a);
		} else if(s.get_vel(a) < 0 && dN2 <= 0) {
			//if the velocity of the sphere in this axis is negative and the position of the sphere
			//is greater than the lower cell boundary position in this axis then calculate the time until the collision
			//with this boundary
			t = dN2 / s.get_vel(a);
		}

		if(t >= 0 && min_t < 0) {
			min_t = t; //if this is the first time calculated set the it as the minimum time
		} else if(t >= 0 && min_t > 0) {
			if(min_t > t) {
				min_t = t; //if the time calcualted it is greater than zero and the minimum time is greater than zero the check to see if the new time is less than the current minimum time
			}
		}
	}

	return (min_t + s.get_time());
}

double sim::calculate_sphere_event(const int i, const int j, const int k, const int l, const int m, const int n) {
	sphere &s1 = get_sphere(i, j, k);
	sphere &s2 = get_sphere(l, m, n);

	double max_t = max<double>(s1.get_time(), s2.get_time()); //find the max time
	double dt1 = max_t - s1.get_time();
	double dt2 = max_t - s2.get_time();

	//advance both spheres to the same max time
	double advPos1[3], advPos2[3];

	advPos1[0] = s1.get_pos(0) + s1.get_vel(0)*dt1;
	advPos2[0] = s2.get_pos(0) + s2.get_vel(0)*dt2;

	advPos1[1] = s1.get_pos(1) + s1.get_vel(1)*dt1;
	advPos2[1] = s2.get_pos(1) + s2.get_vel(1)*dt2;

	advPos1[2] = s1.get_pos(2) + s1.get_vel(2)*dt1;
	advPos2[2] = s2.get_pos(2) + s2.get_vel(2)*dt2;

	//translate the advance position of sphere 2 according to its relative position defined by the id [l,m,n]
	if(l < 0) {advPos2[0] = advPos2[0] - dim[0]*cell_dim[0];} 
	else if (l >= dim[0]) {advPos2[0] = advPos2[0] + dim[0]*cell_dim[0];}

	if(m < 0) {advPos2[1] = advPos2[1] - dim[1]*cell_dim[1];} 
	else if (m >= dim[1]) {advPos2[1] = advPos2[1] + dim[1]*cell_dim[1];}

	if(n < 0) {advPos2[2] = advPos2[2] - dim[2]*cell_dim[2];} 
	else if (n >= dim[2]) {advPos2[2] = advPos2[2] + dim[2]*cell_dim[2];}

	//calculate the vector difference between the advanced positions and the velocities of the spheres
	double p[3], v[3], np[3], P(0), V;
	for(int a = 0; a < 3; a++) {
		p[a] = advPos2[a] - advPos1[a];
		v[a] = s2.get_vel(a) - s1.get_vel(a);
		P += p[a]*p[a];
	}

	double N = sqrt(P);

	//normalise p
	np[0] = p[0] / N;
	np[1] = p[1] / N;
	np[2] = p[2] / N;

	//calculate the magnitude velocity along np
	V = v[0]*np[0] + v[1]*np[1] + v[2]*np[2];

	double D_max_t = diameter(max_t); //the diameter of the spheres at t(maximum)
	double dN = round(sqrt(P) - D_max_t);

	
	if(dN < 0) {
		cout << "INTERSECTION/TOUCH: " << dN << " | V: " << V << endl;
		//s1.print();
		//s2.print();
		//pause();
	}
	

	//only calculate the time of the collision if the spheres are not intersecting and moving towards each other
	//or if they are intersecting and moving away from each other. Otherwise return -1 (indicating a no event
	if((dN >= 0 && V < growth_rate()) || (dN < 0 && V > 0)) {
	
		//calculate the time till collision, dt = t(collision) - t(maximum), by solving the quadratic equation
		double A, B, C, discriminant, E, root1(-1), root2(-1);

		A = v[0]*v[0] + v[1]*v[1] + v[2]*v[2] - growth_rate()*growth_rate(); //coefficient of the t^2 term
		B = p[0]*v[0] + p[1]*v[1] + p[2]*v[2] -  D_max_t*growth_rate(); //coefficient of the t^2 term
		C = P - D_max_t*D_max_t; //coefficient of the t^0 term;
		discriminant = B*B - A*C;

		//calculate the roots of the quadratic equation
		if(A != 0) {
			if(discriminant >= 0) {
				E = sqrt(discriminant);
				root1 = (E - B) / A;
				root2 = -(B + E) / A;
			}
		} else if(B != 0) {
			root1 = -C/B;
		}

		//return the minimum non negative root [+ t(maximum)], if none exist then return -1
		if(root1 >= 0 && root2 >= 0) {
			if(root1 <= root2) {
				return (root1 + max_t);
			} else {
				return (root2 + max_t);
			}
		} else if(root1 >= 0) {
			return (root1 + max_t);
		} else if(root2 >= 0) {
			return (root2 + max_t);
		} else {
			return -1;
		}
	} else {
		return -1;
	}
}

//advances the simulation
void sim::advance(cEvent *ev) {
	int *id1, *id2;
	bool valid_event = true;
	int sphere_count = ev->get_sphere_count();
	double dt;

	current_time = ev->get_time(); //advance the current time to the event time

	for(int a = 0; a < sphere_count; a++) {
		id1 = ev->get_id(a);
		sphere &s = get_sphere(id1[0], id1[1], id1[2]);

		if(valid_event && (s.get_event_count() - ev->get_event_count(a) != 0)) {
			valid_event = false;
		}
		
		//advance the sphere to the current time
		dt = current_time - s.get_time();
		s.set_pos(s.get_pos(0) + s.get_vel(0)*dt, 
				  s.get_pos(1) + s.get_vel(1)*dt,
				  s.get_pos(2) + s.get_vel(2)*dt);
		s.set_time(current_time);
	}

	if(valid_event) {
		if(accepted_counter%5000 == 0) {
			cout << "T: " << current_time << " | QS: " << event_q.size() << " | #E: " << accepted_counter << " | #D: " << discarded_counter << " | #SEV: " << sev_counter << " | #BEV: " << bev_counter << endl;
			make_concurrent();
			save("save.txt");
		}

		accepted_counter++;
		id1 = ev->get_id(0);

		if(sphere_count == 1) {
			calculate_boundary_velocity(id1[0], id1[1], id1[2]); //calculate the new velocity of the sphere
			bev_counter++;
		} else if(sphere_count > 1) {
			for(int a = 1; a < sphere_count; a++) {
				id2 = ev->get_id(a);
				calculate_sphere_velocity(id1[0], id1[1], id1[2], id2[0], id2[1], id2[2]); //calculate the new velocity of bothe spheres
			}
			get_sphere(id1[0], id1[1], id1[2]).inc_events();
			sev_counter++;
		}

		id1 = ev->get_id(0);
		calculate_events(id1[0], id1[1], id1[2]);

		for(int a = 1; a < sphere_count; a++) {
			id1 = ev->get_id(a); //get the id of the sphere relative to primary sphere
			id2 = get_sphere(id1[0], id1[1], id1[2]).get_id(); //get the actual id
			calculate_events(id2[0], id2[1], id2[2]); //calculate the new events for the sphere
		}
	} else {
		discarded_counter++;
	}
}

void sim::calculate_boundary_velocity(const int i, const int j, const int k) {
	sphere &s = get_sphere(i, j, k);

	//if the sphere is at (or outside) any of the cell boundaries and the velocity in that axis is away from the cell then reverse
	//the direction of its velocity in that axis
	for(int a = 0; a < 3; a++) {
		if(round(s.get_pos(a) - s.get_lower_lim(a)) <= 0) {
			if(s.get_vel(a) < 0) {
				s.set_vel(a, -1*s.get_vel(a));
			}
		} else if (round(s.get_upper_lim(a) - s.get_pos(a)) <= 0) {
			if(s.get_vel(a) > 0) {
				s.set_vel(a, -1*s.get_vel(a));
			}
		}
	}

	s.inc_events(); //increment the sphere's number of events
}

void sim::calculate_sphere_velocity(const int i, const int j, const int k, const int l, const int m, const int n) {
	sphere &s1 = get_sphere(i, j, k);
	sphere &s2 = get_sphere(l, m, n);

	double p[3], P(0), np[3], p2[3], v[3], V;

	//translate the advance position of sphere 2 according to its relative position defined by the id [l,m,n]
	if(l < 0) {p2[0] = s2.get_pos(0) - dim[0]*cell_dim[0];} 
	else if (l >= dim[0]) {p2[0] = s2.get_pos(0) + dim[0]*cell_dim[0];}
	else {p2[0] = s2.get_pos(0);}

	if(m < 0) {p2[1] = s2.get_pos(1) - dim[1]*cell_dim[1];} 
	else if (m >= dim[1]) {p2[1] = s2.get_pos(1) + dim[1]*cell_dim[1];}
	else {p2[1] = s2.get_pos(1);}

	if(n < 0) {p2[2] = s2.get_pos(2) - dim[2]*cell_dim[2];} 
	else if (n >= dim[2]) {p2[2] = s2.get_pos(2) + dim[2]*cell_dim[2];}
	else {p2[2] = s2.get_pos(2);}


	//calculate the vector difference between the advanced positions and the velocities of the spheres
	for(int a = 0; a < 3; a++) {
		p[a] = p2[a] - s1.get_pos(a);
		v[a] = s2.get_vel(a) - s1.get_vel(a);
		P += p[a]*p[a];
	}

	double N = sqrt(P);

	//normalise p
	np[0] = p[0] / N;
	np[1] = p[1] / N;
	np[2] = p[2] / N;

	//calculate the magnitude velocity along np
	V = v[0]*np[0] + v[1]*np[1] + v[2]*np[2];

	if(V < 0) {
		//if the spheres are moving towards each other swap their velocities along np
		for(int a = 0; a < 3; a++) {
			s1.set_vel(a, s1.get_vel(a) + V*np[a]);
			s2.set_vel(a, s2.get_vel(a) - V*np[a]);
		}
	} else if( V >= 0 && V < growth_rate()) {
		//if the spheres are moving apart but slower than the growth rate (or not at all) increase their velocities
		//in the np direction such that the difference in their velocities in this direction is equal to the growth rate
		//i.e. effectively they are push each other away
		for(int a = 0; a < 3; a++) {
			s1.set_vel(a, s1.get_vel(a) - growth_rate()*np[a]/2);
			s2.set_vel(a, s1.get_vel(a) + growth_rate()*np[a]/2);
		}
	}
	//if the they are moving away from each other at a speed faster than the growth rate along np then do nothing to the velocities

	s2.inc_events();
}

void sim::make_concurrent() {
	double dt;

	for(int a = 0; a < dim[0]; a++) {
		for(int b = 0; b < dim[1]; b++) {
			for(int c = 0; c < dim[2]; c++) {
				sphere &s = get_sphere(a, b, c);
				//advance the sphere to the current time
				dt = current_time - s.get_time();
				s.set_pos(s.get_pos(0) + s.get_vel(0)*dt, 
						  s.get_pos(1) + s.get_vel(1)*dt,
						  s.get_pos(2) + s.get_vel(2)*dt);
				s.set_time(current_time);
			}
		}
	}
}

//rounds 9 decimal places
double sim::round(double v) {
	double order = pow(10., 9);
	double a = v*order;
	double b = floor(a);
	double c = a-b;
	double d;

	if(c >= 0.5) {
		d = b + 1;
	} else {
		d = b;
	}

	return d / order;
}

void sim::save(char fname[128]) {
	ofstream file(fname);
	int id[3];
	if(file.is_open()) {
		cout << "Saving Data" << endl;
		for(int a = 0; a < dim[0]; a++) {
			id[0] = a;
			for(int b = 0; b < dim[1]; b++) {
				id[1] = b;
				for(int c = 0; c < dim[2]; c++) {
					id[2] = c;
					sphere &s = get_sphere(a, b, c);
					//write the id
					//cout << "Writing ID: ";
					//print_array(id, true);
					file << a << " " << b << " " << c << " ";
					//write the time and number of events
					file << s.get_time() << " " << s.get_event_count() << " ";
					//write the position
					file << s.get_pos(0) << " " << s.get_pos(1) << " " << s.get_pos(2) << " ";
					//write the velocity
					file << s.get_vel(0) << " " << s.get_vel(1) << " " << s.get_vel(2) << " ";

					if(a != dim[0] - 1 || b != dim[1] - 1 || c != dim[2] - 1) {
						file << "\n";
					}
				}
			}
		}

		file.close();
		cout << "Save Complete" << endl;
	} else {
		cout << "File Not Open" << endl;
	}
}

void sim::pause() {
	char kin;
	cin >> kin;
}

#endif