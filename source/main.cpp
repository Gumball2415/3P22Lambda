// main.cpp : Defines the entry point for the application.
//

#include "main.h"
#include "config.h"

int main()
{
	std::cout << "3P22Lambda v."
		<< PROJECT_VERSION << std::endl
		<< "CRT phosphor decay emulator" << std::endl;

	// pause
	std::cout << "Press any key to continue . . . ";
	std::cin.get();
	return 0;
}
