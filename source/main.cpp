#include <string>
#include <iostream>
#include <algorithm>
#include <sstream>
#include "main.h"
#include "config.h"
#include "lodepng/lodepng.h"
#include <tclap/CmdLine.h>

int main(int argc, char* argv[])
{
	std::cout << PROJECT_NAME << std::endl
		<< PROJECT_DESCRIPTION << std::endl
		<< "uses LodePNG, Copyright (c) 2005-2010 Lode Vandevenne" << std::endl;

	try {
		TCLAP::CmdLine cmd("uses TCLAP, Copyright (c) 2017-2021 Google LLC, (c) 2012 - 2016 Daniel Aarno, (c) 2003 - 2012 Michael E.Smoot", ' ', PROJECT_VERSION);
		TCLAP::ValueArg<std::string> inpath_arg(
			"i",
			"input",
			"input path to folder of .png sequence",
			true,
			"",
			"path");
		TCLAP::ValueArg<std::string> outpath_arg(
			"o",
			"output",
			"output path to folder of .png sequence",
			true,
			"",
			"path");

		cmd.add(inpath_arg);
		cmd.add(outpath_arg);

		cmd.parse(argc, argv);

		std::string inpath = inpath_arg.getValue();
		std::string outpath = outpath_arg.getValue();

		std::cout << inpath << std::endl;
		std::cout << outpath << std::endl;

	}
	catch (TCLAP::ArgException& e)
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
	}

	// step 1: count input .png frames
	// step 2: process each frame through ProcessFrame()
	// step 3: reencode processed frames back into video
	// step 4: profit
	return 0;
}

void ProcessFrame()
{
	return;
}
