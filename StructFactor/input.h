#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
//using namespace std;

using std::string;
using std::ifstream;
using std::vector;
using std::iostream;
using std::cerr;
using std::endl;
using std::istringstream;

std::string& trim(string& s);
class Input
{
private:
	static const int LINE_MAX = 2147483647;
	ifstream infile;
	string fname;
	size_t headline;
	size_t linePointer;
	vector<double> data;
	vector<string> str_data;
	size_t totalLine;
	size_t totalUnemptyLine;  //
public:
	Input(string _fname, size_t _headline);
	Input();

	inline const string& get_fname() { return fname; };
	inline const vector<double>& get_data() const {
		if (!data.empty())
		{
			return data;
		}
		else
		{
			cerr << "line " << linePointer << " is EMPTY" << endl;
			throw "trying get empty data";
		}

	}
	inline size_t get_linep() const { return linePointer; }
	inline size_t get_totall() const { return totalLine; }
	inline size_t get_unempty_line() const { return totalUnemptyLine; }
	inline const vector<string>& get_data_str() const { return str_data; }

	inline void reset_filename(string _fname) { fname = _fname; };
	bool open_file();
	void close_file();
	void skiphead();
	size_t move_to_line(size_t _line);
	size_t read_line_data();
	size_t read_line_str(size_t _num);
	size_t skip_line(size_t _num);

	bool check_EOF() const;
};

