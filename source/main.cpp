#include "main.h"
#include "config.h"

int main()
{
	std::cout << PROJECT_NAME <<" v."
		<< PROJECT_VERSION << std::endl
		<< PROJECT_DESCRIPTION << std::endl;

	// pause
	std::cout << "Press enter to continue . . . ";
	std::cin.get();
	return 0;
}
