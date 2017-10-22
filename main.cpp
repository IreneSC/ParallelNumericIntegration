/**
 * @file main.cpp
 * @brief Performs all tests (in Print) and outputs the results
 * @author Irene Crowell
 */
#include "Print.h"
#include "Functions.h"

int main() {
	std::cout.precision(15);
	Functions functions;

	int threads = 4;
	double errorFast = 1e-3;
	double errorSlow = 1e-6;
	int subdivisionsFast = 1e2;
	int subdivisionsSlow = 1e5;
	int timeFast = 1;
	int timeSlow = 5;
	int pointsFast = 3;
	int pointsSlow = 6;
	int keyFast = GSL_INTEG_GAUSS15;
	int keySlow = GSL_INTEG_GAUSS61;

	std::cout << "running..." << std::endl;

	std::cout << std::endl << "NonAdaptiveNonParallel" << std::endl;
	printNonAdaptiveNonParallel(subdivisionsFast, subdivisionsSlow);

	std::cout << std::endl << "NonAdaptiveParallel" << std::endl;
	printNonAdaptiveParallel(subdivisionsFast, subdivisionsSlow, threads);

	std::cout << std::endl << "AdaptiveNonParallel" << std::endl;
	printAdaptiveNonParallel(subdivisionsSlow, timeFast, timeSlow, errorFast,
	      errorSlow);

	std::cout << std::endl << "AdaptiveParallel" << std::endl;
	printAdaptiveParallel(subdivisionsSlow, timeFast, timeSlow, errorFast,
	      errorSlow, threads);
	std::cout << std::endl << "All GSL" << std::endl;
	printAdvanced(subdivisionsSlow, errorFast, errorSlow, pointsFast, pointsSlow,
	      keyFast, keySlow, threads);

	std::cout << "done" << std::endl;
	return 0;
}
