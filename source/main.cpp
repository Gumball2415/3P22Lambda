#include "main.h"
#include "config.h"
#include "lodepng/lodepng.h"

int main(int argc, char* argv[])
{
	std::cout << PROJECT_NAME <<" v."
		<< PROJECT_VERSION << std::endl
		<< PROJECT_DESCRIPTION << std::endl;
	std::cout << "uses LodePNG, Copyright (c) 2005-2010 Lode Vandevenne" << std::endl;

	// TODO: handle arguments

	// pause
	std::cout << "Press enter to continue . . . ";
	std::cin.get();
	return 0;

	// step 1: break down video into .png frames
	// step 2: process each frame through ProcessFrame()
	// step 3: reencode processed frames back into video
	// step 4: profit

}

void ProcessFrame()
{
	return;
}
