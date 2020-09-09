#include "configuration.h"

size_t Configuration::HEAD_INFO_LINE = 22;
size_t Configuration::GAP_LINE = 3;

Configuration::Configuration(std::string config_file)
{
	clog << "Reading configuration data file " << config_file << " ..." << endl;
	//timestep = _time;
	filename = config_file;
	string firstline;
	ifstream infile;
	infile.open(config_file);
	if (infile.is_open())
	{
		clog << "	Processing data file head information..." << endl;
		getline(infile, firstline); //1st line
		string str_timestep;
		for (size_t i = firstline.size() - 1; i > 0; i--)
		{
			if (' ' != firstline[i])
			{
				str_timestep.insert(str_timestep.begin(), firstline[i]);
			}
			else
			{
				break;
			}
		}
		timestep = stoull(str_timestep); //string to unsigned long long
		//cout << timestep << endl;

		infile.ignore(LINE_SKIP_MAX, '\n'); //2nd line

		infile >> particle_num;
		infile.ignore(LINE_SKIP_MAX, '\n'); //3rd line

		infile >> type_num;
		infile.ignore(LINE_SKIP_MAX, '\n'); //4th line

		infile.ignore(LINE_SKIP_MAX, '\n'); //5th line

		infile >> xlo >> xhi;
		infile.ignore(LINE_SKIP_MAX, '\n'); //6th line

		infile >> ylo >> yhi;
		infile.ignore(LINE_SKIP_MAX, '\n'); //7th line

		infile >> zlo >> zhi;
		infile.ignore(LINE_SKIP_MAX, '\n'); //8th line

		infile >> xy >> xz >> yz;
		infile.ignore(LINE_SKIP_MAX, '\n'); //9th line

		infile.ignore(LINE_SKIP_MAX, '\n'); //10th line

		string string_mass_info;
		for (size_t i = 0; i < 4; i++)
		{
			getline(infile, string_mass_info);
			strvec_mass_info.push_back(string_mass_info);
		} //10th - 13th line

		infile.ignore(LINE_SKIP_MAX, '\n'); //14th line

		string string_pair_info;
		for (size_t i = 0; i < 5; i++)
		{
			getline(infile, string_pair_info);
			strvec_pair_info.push_back(string_pair_info);
		} //15th - 19th line
		clog << "	head information has been processed" << endl;
	}
	else
	{
		cerr << "file " << config_file << " open failed!" << endl;
		throw config_file;
	}
	infile.close();

	clog << "	" << firstline << endl;;
	clog << "	configuration data file " << config_file << " has " << particle_num << " particles" << endl;
	clog << "	box parameters:" << endl;
	clog << "		xlo, xhi: " << xlo << ' ' << xhi << endl;
	clog << "		ylo, yhi: " << ylo << ' ' << yhi << endl;
	clog << "		zlo, zhi: " << zlo << ' ' << zhi << endl;
	clog << "		xy, xz, yz: " << xy << ' ' << xz << ' ' << yz << endl;
	//clog << "head line info read!" << endl;

	Input in_data(config_file, HEAD_INFO_LINE);
	clog << endl;
	clog << "	Reading coordinates..." << endl;
	in_data.open_file();
	in_data.skiphead();
	vec_particle.resize(particle_num);
	for (size_t i = 0; i < particle_num; i++)
	{
		in_data.read_line_data();
		vec_particle[i].id = (size_t)in_data.get_data()[0];
		vec_particle[i].type = (int)in_data.get_data()[1];
		vec_particle[i].rx = in_data.get_data()[2];
		vec_particle[i].ry = in_data.get_data()[3];
		vec_particle[i].rz = in_data.get_data()[4];
		vec_particle[i].box_x = (int)in_data.get_data()[5];
		vec_particle[i].box_y = (int)in_data.get_data()[6];
		vec_particle[i].box_z = (int)in_data.get_data()[7];
	}

	in_data.skip_line(GAP_LINE);
	clog << "	Coordinates have been read!" << endl;
	clog << "	Reading velocities..." << endl;
	Particle* p_particle = nullptr;
	for (size_t i = 0; i < particle_num; i++)
	{
		in_data.read_line_data();
		p_particle = seek_id(vec_particle, (size_t)in_data.get_data()[0]);
		(*p_particle).vx = in_data.get_data()[1];
		(*p_particle).vy = in_data.get_data()[2];
		(*p_particle).vz = in_data.get_data()[3];
	}
	p_particle = nullptr;
	clog << "	Velocities have been read!" << endl;
	clog << "Configuration data file " << config_file << " has been read!" << endl;
	clog << endl;
	infile.close();
}

const vector<double> Configuration::compute_RDF(size_t rdf_size, double r_cut)
{
	//double r_end = r_delta * rdf_size;
	//int rdf_size = r_end / r_delta;
	double r_delta = r_cut / rdf_size;
	double lx = xhi - xlo; //box length in x direction
	double ly = yhi - ylo;
	double lz = zhi - zlo;
	if (r_cut > (lx / 2) || r_cut > (ly / 2) || r_cut > (lz / 2))
	{
		cerr << "WARNING: rdf size too large, rdf can be incorrect" << endl;
	}

	vector<double> rdf(rdf_size, 0);
	for (size_t i = 0; i < particle_num - 1; i++)
	{
		for (size_t j = i + 1; j < particle_num; j++)
		{
			double rx_ij = vec_particle[i].rx - vec_particle[j].rx;
			double ry_ij = vec_particle[i].ry - vec_particle[j].ry;
			double rz_ij = vec_particle[i].rz - vec_particle[j].rz;

			rx_ij = convert_to_period_distance(rx_ij, lx);
			ry_ij = convert_to_period_distance(ry_ij, ly);
			rz_ij = convert_to_period_distance(rz_ij, lz);

			double rij_period = sqrt(rx_ij * rx_ij + ry_ij * ry_ij + rz_ij * rz_ij);
			if (rij_period < r_cut)
			{
				int r_index = (int)(rij_period / r_delta);
				rdf[r_index] ++;
			}
		}
	}


	double Vbox = lx * ly * lz; //volume of system
	double rho = particle_num / Vbox; //density of system 
	for (size_t i = 0; i < rdf.size(); i++)
	{
		double rhi = ((double)i + 1) * r_delta; //r
		double rlo = (double)i * r_delta; //r+¦¤r
		double dV = 4.0 / 3.0 * PI * (std::pow(rhi, 3.0) - std::pow(rlo, 3.0)); //volume of sphere shell
		rdf[i] = 2.0 * rdf[i] / dV / rho / particle_num; //only count half particles, so introduce factor 2
	}
	return rdf;
}




double convert_to_period_distance(double distance, double lbox)
{
	double half_lbox = lbox / 2;
	if (distance > half_lbox)
	{
		return distance - lbox;
	}
	if (distance < -half_lbox)
	{
		return distance + lbox;
	}
	return distance;
}

void Configuration_With_Sq::compute_rij()
{
	boost::progress_timer r_ij_timer;
	double lx = get_xhi() - get_xlo();
	double ly = get_yhi() - get_ylo();
	double lz = get_zhi() - get_zlo();
	vec_r_ij.resize(get_particle_num() * (get_particle_num() - 1) / 2);
	vec_rx_ij.resize(get_particle_num() * (get_particle_num() - 1) / 2);
	vec_ry_ij.resize(get_particle_num() * (get_particle_num() - 1) / 2);
	vec_rz_ij.resize(get_particle_num() * (get_particle_num() - 1) / 2);
	vector<double>::iterator iter_vec_r_ij = vec_r_ij.begin();
	vector<double>::iterator iter_vec_rx_ij = vec_rx_ij.begin();
	vector<double>::iterator iter_vec_ry_ij = vec_ry_ij.begin();
	vector<double>::iterator iter_vec_rz_ij = vec_rz_ij.begin();

	//size_t temp_counter = 0;
	cout << "computing r_ij..." << endl;
	boost::progress_display progress_displayer(get_particle_num() * (get_particle_num() - 1) / 2);

	for (size_t i = 0; i < get_particle().size() - 1; i++)
	{
		for (size_t j = i + 1; j < get_particle().size(); j++)
		{

			double rx_ij = get_particle()[i].rx - get_particle()[j].rx;
			double ry_ij = get_particle()[i].ry - get_particle()[j].ry;
			double rz_ij = get_particle()[i].rz - get_particle()[j].rz;
			*iter_vec_rx_ij = convert_to_period_distance(rx_ij, lx);
			*iter_vec_ry_ij = convert_to_period_distance(ry_ij, ly);
			*iter_vec_rz_ij = convert_to_period_distance(rz_ij, lz);
			double r_ij = sqrt(pow(*iter_vec_rx_ij, 2) + pow(*iter_vec_ry_ij, 2) + pow(*iter_vec_rz_ij, 2));
			*iter_vec_r_ij = r_ij;
			iter_vec_r_ij++;
			//*iter_vec_rx_ij = rx_ij;
			iter_vec_rx_ij++;
			//*iter_vec_ry_ij = ry_ij;
			iter_vec_ry_ij++;
			//*iter_vec_rz_ij = rz_ij;
			iter_vec_rz_ij++;
			//temp_counter++;
			//cout << temp_counter << endl;
			++progress_displayer;
		}
	}
	rij_update_flag = true;
	cout << "r_ij computing finished, time used: ";
}

const vector<double> Configuration_With_Sq::compute_struct_factor(double q_start, double q_delta, size_t q_size, string direction)
{
	boost::progress_timer sq_timer;
	if (direction != "x" && direction != "y" && direction != "z" && direction != "all")
	{
		cerr << "wrong dirction parameter when computing structure factor!" << endl;
		throw "wrong dirction parameter when computing structure factor!";
	}

	if (!rij_update_flag)
	{
		cerr << "WARNING: RIJ VECTOR MAY NOT BE UPDATED!" << endl;
	}

	for (size_t i = 0; i < q_size; i++)
	{
		vec_q.push_back(q_start + q_delta * i);
	}

	cout << "computing structure factor..." << endl;
	vector<double> vec_struct_factor(q_size);
	boost::progress_display progress_displayer(vec_struct_factor.size());
	if (direction == "all")
	{

		for (size_t q_index = 0; q_index < vec_struct_factor.size(); q_index++)
		{
			double struct_factor = 0;
			double q = q_start + q_delta * q_index;
			//double q = get_xhi() - get_xlo()
			//double qx, qy, qz;
			//qx = qy = qz = q / sqrt(3);
			for (size_t i = 0; i < vec_r_ij.size(); i++)
			{
				//struct_factor += cos(qx * vec_rx_ij[i] + qy * vec_ry_ij[i] + qz * vec_rz_ij[i]);
				struct_factor += std::sin(q * vec_r_ij[i]) / q / vec_r_ij[i];
			}
			struct_factor *= 2.0;
			struct_factor = struct_factor / get_particle_num() + 1;
			//cout << struct_factor_real << endl;
			vec_struct_factor[q_index] = struct_factor;
			++progress_displayer;
		}
	}
	if (direction == "x")
	{
		for (size_t q_index = 0; q_index < vec_struct_factor.size(); q_index++)
		{
			double struct_factor = 0;
			double q = q_start + q_delta * q_index;
			for (size_t i = 0; i < vec_rx_ij.size(); i++)
			{
				//struct_factor += pow(cos(q * vec_rx_ij[i]), 2);
				struct_factor += std::sin(q * vec_rx_ij[i]) / q / vec_rx_ij[i];
			}
			struct_factor *= 2.0;
			struct_factor = struct_factor / get_particle_num() + 1;
			vec_struct_factor[q_index] = struct_factor;
			++progress_displayer;
		}
	}

	if (direction == "y")
	{
		for (size_t q_index = 0; q_index < vec_struct_factor.size(); q_index++)
		{
			double struct_factor = 0;
			double q = q_start + q_delta * q_index;
			for (size_t i = 0; i < vec_ry_ij.size(); i++)
			{
				struct_factor += std::cos(q * vec_ry_ij[i]);
			}
			struct_factor *= 2.0;
			struct_factor = struct_factor / get_particle_num() + 1;
			vec_struct_factor[q_index] = struct_factor;
			++progress_displayer;
		}
	}

	if (direction == "z")
	{
		for (size_t q_index = 0; q_index < vec_struct_factor.size(); q_index++)
		{
			double struct_factor = 0;
			double q = q_start + q_delta * q_index;
			for (size_t i = 0; i < vec_rz_ij.size(); i++)
			{
				struct_factor += std::cos(q * vec_rz_ij[i]);
			}
			struct_factor *= 2.0;
			struct_factor = struct_factor / get_particle_num() + 1;
			vec_struct_factor[q_index] = struct_factor;
			++progress_displayer;
		}
	}
	cout << "struct factor computing finished, time used: ";
	return vec_struct_factor;
}
