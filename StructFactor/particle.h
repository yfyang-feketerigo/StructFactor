//new struct for particle information
//yfyang@ciac.ac.cn
#pragma once
#include <vector>
#include <iostream>
#include <string>
using std::iostream;
using std::string;
using std::vector;

struct Particle
{
	size_t id = 0; //particle id
	int type = 0; //particle type

	double rx = 0;		//position in x direction
	double ry = 0;		//position in y direction
	double rz = 0;		//position in z direction

	int box_x = 0;		//box counter in x direction
	int box_y = 0;		//box counter in y direction
	int box_z = 0;		//box counter in z direction

	double vx = 0;		//velocity in x direction
	double vy = 0;		//velocity in y direction
	double vz = 0;		//velocity in z direction
};

Particle* seek_id(vector<Particle>& vec_p, size_t _id); //return the particle of given id;
