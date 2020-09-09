#include "configuration.h"
#include "input.h"
#include "particle.h"
#include <boost/timer.hpp>
#include <boost/progress.hpp>
#include <fstream>
#include <iostream>
int main()
{
	std::ofstream ofile("test.csv");
	if (!ofile.is_open())
	{
		cerr << "OPEN OUTPUT FILE FAILED!" << endl;
		throw "OPEN OUTPUT FILE FAILED!";
	}
	boost::timer time_counter;
	Configuration_With_Sq::SET_GAP_LINE(3);
	Configuration_With_Sq::SET_HEAD_INFO_LINE(22);
	Configuration_With_Sq config_test("data.restart.shear.wi.1.4800");
	config_test.compute_rij();
	//cout << config_test.get_particle_num() << endl;
	vector<double> test = config_test.compute_struct_factor(2, 0.01, 2000, "x");
	//cout << test.size() << endl;
	//cout << "Sq time used: " << time_counter.elapsed() << "s" << endl;
	for (size_t i = 0; i < test.size(); i++)
	{
		ofile << test[i] << endl;
	}
	cout << "total time used: " << time_counter.elapsed() << "s" << endl;
	return 0;
}
