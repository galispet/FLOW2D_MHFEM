#pragma once


#include <iostream>
#include <string>

class GNUPLOT {

private:

	FILE * gnuplotpipe;

public:

	GNUPLOT();
	~GNUPLOT();

	void operator()(std::string const & command);

};

GNUPLOT::GNUPLOT() {

	gnuplotpipe = _popen("C:\\gnuplot\\bin\\gnuplot.exe", "w");

	if (!gnuplotpipe)
		std::cerr << ("Gnuplot was not found.");

};

GNUPLOT::~GNUPLOT() {

	fprintf(gnuplotpipe, "exit\n");
	_pclose(gnuplotpipe);

};

void GNUPLOT::operator()(std::string const & command) {

	fprintf(gnuplotpipe, "%s\n", command.c_str());
	fflush(gnuplotpipe);

};