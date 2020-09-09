#include "input.h"


Input::Input(string _fname, size_t _headline)
{
	fname = _fname;
	headline = _headline;
	linePointer = 0;
	infile.open(fname);
	totalLine = 0;
	string line;
	if (infile.is_open())
	{
		while (!infile.eof())
		{
			totalLine += 1;
			getline(infile, line);
			if (!line.empty())
			{
				totalUnemptyLine += 1;
			}
		}
		infile.close();
	}
	else
	{
		totalLine = 0;
		cerr << "open input file failed!" << endl;
		throw "open input file failed!";
	}
}

Input::Input()
{
	fname = "default";
	headline = 0;
	linePointer = 0;
	infile.open(fname);
	totalLine = 0;
	totalUnemptyLine = 0;
}

bool Input::open_file()
{
	if (!infile.is_open())
	{
		infile.open(fname);
	}
	else
	{
		cerr << "input file has already been opened!" << endl;
		throw "input file has already been opened!";
	}
	return infile.is_open();
}

void Input::close_file()
{
	if (infile.is_open())
	{
		infile.close();
		linePointer = 0;
		totalLine = 0;
		data.clear();
	}
	else
	{
		cerr << "input file is not opened!" << endl;
		throw "input file is not opened!";
	}
}

void Input::skiphead()
{
	if (infile.is_open())
	{
		if (linePointer == 0)
		{
			for (size_t i = 0; i < headline; i++)
			{
				infile.ignore(LINE_MAX, '\n');
				linePointer += 1;
			}
		}
		else
		{
			cerr << "head lines in file has already been skipped!" << endl;
			throw "head lines in file has already been skipped!";
		}
	}
	else
	{
		cerr << fname << " not open!" << endl;
		throw(fname + " not open!");
	}
}

size_t Input::move_to_line(size_t _line)
{
	if (_line > totalLine)
	{
		cerr << "out of file line range!" << endl;
		return linePointer;
	}
	infile.close();
	infile.open(fname);
	for (size_t i = 1; i < _line; i++)
	{
		infile.ignore(LINE_MAX, '\n');
	}
	return linePointer = _line;
}

size_t Input::read_line_data()
{
	data.clear();
	string line;
	if (infile.is_open())
	{
		if (infile.eof())
		{
			cerr << "reach eof!" << endl;
			cerr << "data remain unchanged!" << endl;
			return 0;
		}
		else
		{
			do
			{
				getline(infile, line);
				if (line.empty())
				{
					cerr << "skip reading EMPTY data line " << linePointer << endl;
					linePointer += 1;
				}
			} while (line.empty() && !infile.eof());
			istringstream ss(line);
			while (!ss.eof())
			{
				double tmp;
				ss >> tmp;
				data.push_back(tmp);
			}
		}
		linePointer += 1;
		return data.size();
	}

	else
	{
		cerr << "input file is not open!" << endl;
		throw "input file is not open!";
	}
}

size_t Input::read_line_str(size_t _num)
{
	string str_temp;
	if (!str_data.empty())
	{
		str_data.clear();
	}
	if (infile.is_open())
	{
		for (size_t i = 0; i < _num; i++)
		{
			getline(infile, str_temp);
			str_data.push_back(str_temp);
			linePointer++;
		}
	}
	return str_data.size();
}

size_t Input::skip_line(size_t _num)
{
	if (infile.is_open())
	{
		for (size_t i = 0; i < _num; i++)
		{
			infile.ignore(LINE_MAX, '\n');
			linePointer++;
		}
	}
	else
	{
		cerr << "file not open!" << endl;
	}
	return linePointer;
}

bool Input::check_EOF() const
{
	if (infile.eof())
	{
		return true;
	}
	else
	{
		return false;
	}
}

string& trim(string& s)
{
	if (s.empty())
		return s;
	s.erase(0, s.find_first_not_of(' '));
	s.erase(s.find_last_not_of(' ') + 1);
	return s;
}
