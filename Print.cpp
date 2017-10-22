/**
 * @file Print.cpp
 * @brief %Functions to print test output to csv files
 * @author Irene Crowell
 */
#include "Print.h"

void printNonAdaptiveNonParallel(int subdivisionsFast, int subdivisionsSlow) {
	std::cout.precision(15);
	Functions functions;
	int subdivisions;

	std::fstream file;
	for (int i = 0; i < 2; i++) { //run two tests

		if (i == 0) {
			file.open("TestData/nonAdaptiveNonParallelFast.csv",
			      std::fstream::out);
			subdivisions = subdivisionsFast;
		} else {
			file.open("TestData/nonAdaptiveNonParallelAccurate.csv",
			      std::fstream::out);
			subdivisions = subdivisionsSlow;
		}

		file << "Subdivisions: " << subdivisions << "\n";
		file << ",Type,Integral,Result,Error,Time,\n";
		file << "Midpoint Rule:\n";
		for (Functions::integrableFunction &function : functions.functions) {
			std::cout << "Calculating " << function.name << "... " << std::flush;
			file << "," << std::defaultfloat << function.type << ","
			      << function.name << " from " << function.a << " to "
			      << function.b;
			std::clock_t start = std::clock();
			double value = MidpointRule::nonAdaptiveNonParallel(function.f,
			      function.a, function.b, subdivisions);
			double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
			file << "," << std::fixed << value;
			file << "," << fabs(value - function.value);
			file << "," << duration << std::endl;
			std::cout << "done." << std::endl;
		}

		file << "\nTrapezoid Rule:\n";
		for (Functions::integrableFunction &function : functions.functions) {
			std::cout << "Calculating " << function.name << "... " << std::flush;
			file << "," << std::defaultfloat << function.type << ","
			      << function.name << " from " << function.a << " to "
			      << function.b;
			std::clock_t start = std::clock();
			double value = TrapezoidRule::nonAdaptiveNonParallel(function.f,
			      function.a, function.b, subdivisions);
			double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
			file << "," << std::fixed << value;
			file << "," << fabs(value - function.value);
			file << "," << duration << std::endl;
			std::cout << "done." << std::endl;
		}

		file << "\nSimpson Rule:\n";
		for (Functions::integrableFunction &function : functions.functions) {
			std::cout << "Calculating " << function.name << "... " << std::flush;
			file << "," << std::defaultfloat << function.type << ","
			      << function.name << " from " << function.a << " to "
			      << function.b;
			std::clock_t start = std::clock();
			double value = SimpsonRule::nonAdaptiveNonParallel(function.f,
			      function.a, function.b, subdivisions);
			double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
			file << "," << std::fixed << value;
			file << "," << fabs(value - function.value);
			file << "," << duration << std::endl;
			std::cout << "done." << std::endl;
		}

		file << "\nSimpson 3/8 Rule:\n";
		for (Functions::integrableFunction &function : functions.functions) {
			std::cout << "Calculating " << function.name << "... " << std::flush;
			file << "," << std::defaultfloat << function.type << ","
			      << function.name << " from " << function.a << " to "
			      << function.b;
			std::clock_t start = std::clock();
			double value = Simpson38Rule::nonAdaptiveNonParallel(function.f,
			      function.a, function.b, subdivisions);
			double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
			file << "," << std::fixed << value;
			file << "," << fabs(value - function.value);
			file << "," << duration << std::endl;
			std::cout << "done." << std::endl;
		}

		file << "\nBoole's Rule:\n";
		for (Functions::integrableFunction &function : functions.functions) {
			std::cout << "Calculating " << function.name << "... " << std::flush;
			file << "," << std::defaultfloat << function.type << ","
			      << function.name << " from " << function.a << " to "
			      << function.b;
			std::clock_t start = std::clock();
			double value = BoolesRule::nonAdaptiveNonParallel(function.f,
			      function.a, function.b, subdivisions);
			double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
			file << "," << std::fixed << value;
			file << "," << fabs(value - function.value);
			file << "," << duration << std::endl;
			std::cout << "done." << std::endl;
		}
		file.close();
	}
}

void printNonAdaptiveParallel(int subdivisionsFast, int subdivisionsSlow,
      int threads) {
	std::cout.precision(15);
	Functions functions;
	int subdivisions;

	std::fstream file;
	for (int i = 0; i < 2; i++) { //run two tests

		if (i == 0) {
			file.open("TestData/nonAdaptiveParallelFast.csv", std::fstream::out);
			subdivisions = subdivisionsFast;
		} else {
			file.open("TestData/nonAdaptiveParallelAccurate.csv",
			      std::fstream::out);
			subdivisions = subdivisionsSlow;
		}

		file << "Subdivisions: " << subdivisions << "\n";
		file << ",Type,Integral,Result,Error,Time,\n";
		file << "Midpoint Rule:\n";
		for (Functions::integrableFunction &function : functions.functions) {
			std::cout << "Calculating " << function.name << "... " << std::flush;
			file << "," << std::defaultfloat << function.type << ","
			      << function.name << " from " << function.a << " to "
			      << function.b;
			std::clock_t start = std::clock();
			double value = MidpointRule::nonAdaptiveParallel(function.f,
			      function.a, function.b, subdivisions, threads);
			double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
			file << "," << std::fixed << value;
			file << "," << fabs(value - function.value);
			file << "," << duration << std::endl;
			std::cout << "done." << std::endl;
		}

		file << "\nTrapezoid Rule:\n";
		for (Functions::integrableFunction &function : functions.functions) {
			std::cout << "Calculating " << function.name << "... " << std::flush;
			file << "," << std::defaultfloat << function.type << ","
			      << function.name << " from " << function.a << " to "
			      << function.b;
			std::clock_t start = std::clock();
			double value = TrapezoidRule::nonAdaptiveParallel(function.f,
			      function.a, function.b, subdivisions, threads);
			double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
			file << "," << std::fixed << value;
			file << "," << fabs(value - function.value);
			file << "," << duration << std::endl;
			std::cout << "done." << std::endl;
		}

		file << "\nSimpson Rule:\n";
		for (Functions::integrableFunction &function : functions.functions) {
			std::cout << "Calculating " << function.name << "... " << std::flush;
			file << "," << std::defaultfloat << function.type << ","
			      << function.name << " from " << function.a << " to "
			      << function.b;
			std::clock_t start = std::clock();
			double value = SimpsonRule::nonAdaptiveParallel(function.f, function.a,
			      function.b, subdivisions, threads);
			double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
			file << "," << std::fixed << value;
			file << "," << fabs(value - function.value);
			file << "," << duration << std::endl;
			std::cout << "done." << std::endl;
		}

		file << "\nSimpson 3/8 Rule:\n";
		for (Functions::integrableFunction &function : functions.functions) {
			std::cout << "Calculating " << function.name << "... " << std::flush;
			file << "," << std::defaultfloat << function.type << ","
			      << function.name << " from " << function.a << " to "
			      << function.b;
			std::clock_t start = std::clock();
			double value = Simpson38Rule::nonAdaptiveParallel(function.f,
			      function.a, function.b, subdivisions, threads);
			double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
			file << "," << std::fixed << value;
			file << "," << fabs(value - function.value);
			file << "," << duration << std::endl;
			std::cout << "done." << std::endl;
		}

		file << "\nBoole's Rule:\n";
		for (Functions::integrableFunction &function : functions.functions) {
			std::cout << "Calculating " << function.name << "... " << std::flush;
			file << "," << std::defaultfloat << function.type << ","
			      << function.name << " from " << function.a << " to "
			      << function.b;
			std::clock_t start = std::clock();
			double value = BoolesRule::nonAdaptiveParallel(function.f, function.a,
			      function.b, subdivisions, threads);
			double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
			file << "," << std::fixed << value;
			file << "," << fabs(value - function.value);
			file << "," << duration << std::endl;
			std::cout << "done." << std::endl;
		}
		file.close();
	}
}

void printAdaptiveNonParallel(int max_subdivisions, int timeFast, int timeSlow,
      double errorFast, double errorSlow) {
	std::cout.precision(15);
	Functions functions;
	double error;
	int time;
	int subdivisions;

	std::fstream file;
	for (int i = 0; i < 2; i++) { //run two tests

		if (i == 0) {
			file.open("TestData/adaptiveNonParallelFast.csv", std::fstream::out);
			error = errorFast;
			time = timeFast;
		} else {
			file.open("TestData/adaptiveNonParallelAccurate.csv",
			      std::fstream::out);
			error = errorSlow;
			time = timeSlow;
		}
		file << "Error Goal: " << error;
		file << ",max_subdivisions: " << max_subdivisions << "\n";
		file << ",Type,Integral,Result,Error,Time,Subdivisions\n";
		file << "Midpoint Rule:\n";
		for (Functions::integrableFunction &function : functions.functions) {
			std::cout << "Calculating " << function.name << "... " << std::flush;
			file << "," << std::defaultfloat << function.type << ","
			      << function.name << " from " << function.a << " to "
			      << function.b;
			std::clock_t start = std::clock();
			double value = MidpointRule::adaptiveNonParallel(function.f,
			      function.a, function.b, error, max_subdivisions, time,
			      &subdivisions);
			double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
			file << "," << std::fixed << value;
			file << "," << fabs(value - function.value);
			file << "," << duration << "," << subdivisions << std::endl;
			std::cout << "done." << std::endl;
		}

		file << "\nTrapezoid Rule:\n";
		for (Functions::integrableFunction &function : functions.functions) {
			std::cout << "Calculating " << function.name << "... " << std::flush;
			file << "," << std::defaultfloat << function.type << ","
			      << function.name << " from " << function.a << " to "
			      << function.b;
			std::clock_t start = std::clock();
			double value = TrapezoidRule::adaptiveNonParallel(function.f,
			      function.a, function.b, error, max_subdivisions, time,
			      &subdivisions);
			double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
			file << "," << std::fixed << value;
			file << "," << fabs(value - function.value);
			file << "," << duration << "," << subdivisions << std::endl;
			std::cout << "done." << std::endl;
		}

		file << "\nSimpson Rule:\n";
		for (Functions::integrableFunction &function : functions.functions) {
			std::cout << "Calculating " << function.name << "... " << std::flush;
			file << "," << std::defaultfloat << function.type << ","
			      << function.name << " from " << function.a << " to "
			      << function.b;
			std::clock_t start = std::clock();
			double value = SimpsonRule::adaptiveNonParallel(function.f, function.a,
			      function.b, error, max_subdivisions, time, &subdivisions);
			double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
			file << "," << std::fixed << value;
			file << "," << fabs(value - function.value);
			file << "," << duration << "," << subdivisions << std::endl;
			std::cout << "done." << std::endl;
		}

		file << "\nSimpson 3/8 Rule:\n";
		for (Functions::integrableFunction &function : functions.functions) {
			std::cout << "Calculating " << function.name << "... " << std::flush;
			file << "," << std::defaultfloat << function.type << ","
			      << function.name << " from " << function.a << " to "
			      << function.b;
			std::clock_t start = std::clock();
			double value = Simpson38Rule::adaptiveNonParallel(function.f,
			      function.a, function.b, error, max_subdivisions, time,
			      &subdivisions);
			double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
			file << "," << std::fixed << value;
			file << "," << fabs(value - function.value);
			file << "," << duration << "," << subdivisions << std::endl;
			std::cout << "done." << std::endl;
		}

		file << "\nBoole's Rule:\n";
		for (Functions::integrableFunction &function : functions.functions) {
			std::cout << "Calculating " << function.name << "... " << std::flush;
			file << "," << std::defaultfloat << function.type << ","
			      << function.name << " from " << function.a << " to "
			      << function.b;
			std::clock_t start = std::clock();
			double value = BoolesRule::adaptiveNonParallel(function.f, function.a,
			      function.b, error, max_subdivisions, time, &subdivisions);
			double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
			file << "," << std::fixed << value;
			file << "," << fabs(value - function.value);
			file << "," << duration << "," << subdivisions << std::endl;
			std::cout << "done." << std::endl;
		}
		file.close();
	}
}

void printAdaptiveParallel(int max_subdivisions, int timeFast, int timeSlow,
      double errorFast, double errorSlow, int threads) {
	std::cout.precision(15);
	Functions functions;
	double error;
	int time;
	int subdivisions;

	std::fstream file;
	for (int i = 0; i < 2; i++) { //run two tests

		if (i == 0) {
			file.open("TestData/adaptiveParallelFast.csv", std::fstream::out);
			error = errorFast;
			time = timeFast;
		} else {
			file.open("TestData/adaptiveParallelAccurate.csv", std::fstream::out);
			error = errorSlow;
			time = timeSlow;
		}
		file << "Error Goal: " << error;
		file << ",max_subdivisions: " << max_subdivisions << "\n";
		file << ",Type,Integral,Result,Error,Time,Subdivisions\n";
		file << "Midpoint Rule:\n";
		for (Functions::integrableFunction &function : functions.functions) {
			std::cout << "Calculating " << function.name << "... " << std::flush;
			file << "," << std::defaultfloat << function.type << ","
			      << function.name << " from " << function.a << " to "
			      << function.b;
			std::clock_t start = std::clock();
			double value = MidpointRule::adaptiveParallel(function.f, function.a,
			      function.b, threads, error, max_subdivisions, time,
			      &subdivisions);
			double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
			file << "," << std::fixed << value;
			file << "," << fabs(value - function.value);
			file << "," << duration << "," << subdivisions << std::endl;
			std::cout << "done." << std::endl;
		}

		file << "\nTrapezoid Rule:\n";
		for (Functions::integrableFunction &function : functions.functions) {
			std::cout << "Calculating " << function.name << "... " << std::flush;
			file << "," << std::defaultfloat << function.type << ","
			      << function.name << " from " << function.a << " to "
			      << function.b;
			std::clock_t start = std::clock();
			double value = TrapezoidRule::adaptiveParallel(function.f, function.a,
			      function.b, threads, error, max_subdivisions, time,
			      &subdivisions);
			double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
			file << "," << std::fixed << value;
			file << "," << fabs(value - function.value);
			file << "," << duration << "," << subdivisions << std::endl;
			std::cout << "done." << std::endl;
		}

		file << "\nSimpson Rule:\n";
		for (Functions::integrableFunction &function : functions.functions) {
			std::cout << "Calculating " << function.name << "... " << std::flush;
			file << "," << std::defaultfloat << function.type << ","
			      << function.name << " from " << function.a << " to "
			      << function.b;
			std::clock_t start = std::clock();
			double value = SimpsonRule::adaptiveParallel(function.f, function.a,
			      function.b, threads, error, max_subdivisions, time,
			      &subdivisions);
			double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
			file << "," << std::fixed << value;
			file << "," << fabs(value - function.value);
			file << "," << duration << "," << subdivisions << std::endl;
			std::cout << "done." << std::endl;
		}

		file << "\nSimpson 3/8 Rule:\n";
		for (Functions::integrableFunction &function : functions.functions) {
			std::cout << "Calculating " << function.name << "... " << std::flush;
			file << "," << std::defaultfloat << function.type << ","
			      << function.name << " from " << function.a << " to "
			      << function.b;
			std::clock_t start = std::clock();
			double value = Simpson38Rule::adaptiveParallel(function.f, function.a,
			      function.b, threads, error, max_subdivisions, time,
			      &subdivisions);
			double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
			file << "," << std::fixed << value;
			file << "," << fabs(value - function.value);
			file << "," << duration << "," << subdivisions << std::endl;
			std::cout << "done." << std::endl;
		}

		file << "\nBoole's Rule:\n";
		for (Functions::integrableFunction &function : functions.functions) {
			std::cout << "Calculating " << function.name << "... " << std::flush;
			file << "," << std::defaultfloat << function.type << ","
			      << function.name << " from " << function.a << " to "
			      << function.b;
			std::clock_t start = std::clock();
			double value = BoolesRule::adaptiveParallel(function.f, function.a,
			      function.b, threads, error, max_subdivisions, time,
			      &subdivisions);
			double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
			file << "," << std::fixed << value;
			file << "," << fabs(value - function.value);
			file << "," << duration << "," << subdivisions << std::endl;
			std::cout << "done." << std::endl;
		}
		file.close();
	}
}
/**
 * Checks if the error code is actually an error (not empty etc.)
 * @param error_code the error code to check
 * @return true for error
 */
bool checkError(const char * error_code) {
	return ((error_code != NULL) && (error_code[0] != '\0'));
}

void printAdvanced(int max_subdivisions, double errorFast, double errorSlow,
      int pointsFast, int pointsSlow, int keyFast, int keySlow, int threads) {
	gsl_set_error_handler_off();
	const char * error_code = "";
	Functions functions;
	std::fstream file;
	double error;
	int points;
	int key;
	double abserror;
	std::vector<double> singularPoints;
	for (int i = 0; i < 2; i++) { //run two tests

		if (i == 0) {
			file.open("TestData/Advanced/Fast.csv", std::fstream::out);
			error = errorFast;
			points = pointsFast;
			key = keyFast;
		} else {
			file.open("TestData/Advanced/Slow.csv", std::fstream::out);
			error = errorSlow;
			points = pointsSlow;
			key = keySlow;
		}

		file << "Advanced\n";
		file << "max_subdivisions: " << max_subdivisions << "\n";
		file << ",Error Goal " << error << "\n";
		file << ",Type,Integral,Result,Error,Time, AbsError\n";

		for (int var1 = 0; var1 < AdvancedRules::algorithm::size; var1++) {
			std::clock_t start;
			double value;
			double duration;

			switch (var1) {
			case AdvancedRules::algorithm::a_gaussLegendreFixed:
				file << "gaussLegendreFixed:\n";
				for (Functions::integrableFunction &function : functions.functions) {
					std::cout << "Calculating " << function.name << "... "
					      << std::flush;
					start = std::clock();
					value = AdvancedRules::gaussLegendreFixed(function.f, function.a,
					      function.b, points);
					duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
					file << "," << std::defaultfloat << function.type << ","
					      << function.name << " from " << function.a << " to "
					      << function.b << std::fixed << value;
					file << "," << value;
					file << "," << fabs(value - function.value);
					file << "," << duration << "," << std::endl;
					std::cout << "done." << std::endl;
				}
				break;

			case AdvancedRules::algorithm::a_nonAdaptiveGaussKronrod:
				file << "nonAdaptiveGaussKronrod:\n";
				for (Functions::integrableFunction &function : functions.functions) {
					std::cout << "Calculating " << function.name << "... "
					      << std::flush;
					start = std::clock();
					value = AdvancedRules::nonAdaptiveGaussKronrod(error_code,
					      function.f, function.a, function.b, error, &abserror);
					duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
					file << "," << std::defaultfloat << function.type << ","
					      << function.name << " from " << function.a << " to "
					      << function.b << std::fixed << value;

					if (checkError(error_code)) {
						file << "," << error_code << "," << "," << duration
						      << std::endl;
						error_code = "";
					} else {
						file << "," << value;
						file << "," << fabs(value - function.value);
						file << "," << duration << "," << abserror << std::endl;
					}
					std::cout << "done." << std::endl;
				}
				break;

			case AdvancedRules::algorithm::a_gaussLegendreFixedParallel:
				for (Functions::integrableFunction &function : functions.functions) {
					std::cout << "Calculating " << function.name << "... "
					      << std::flush;
					start = std::clock();
					value = AdvancedRules::gaussLegendreFixedParallel(function.f,
					      function.a, function.b, points, threads);
					duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
					file << "," << std::defaultfloat << function.type << ","
					      << function.name << " from " << function.a << " to "
					      << function.b << std::fixed << value;
					file << "," << value;
					file << "," << fabs(value - function.value);
					file << "," << duration << std::endl;
					std::cout << "done." << std::endl;
				}
				break;

			case AdvancedRules::algorithm::a_nonAdaptiveGaussKronrodParallel:
				for (Functions::integrableFunction &function : functions.functions) {
					std::cout << "Calculating " << function.name << "... "
					      << std::flush;
					start = std::clock();
					value = AdvancedRules::nonAdaptiveGaussKronrodParallel(
					      error_code, function.f, function.a, function.b, error,
					      threads, &abserror);
					duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
					file << "," << std::defaultfloat << function.type << ","
					      << function.name << " from " << function.a << " to "
					      << function.b << std::fixed << value;
					if (checkError(error_code)) {
						file << "," << error_code << "," << "," << duration
						      << std::endl;
						error_code = "";
					} else {
						file << "," << value;
						file << "," << fabs(value - function.value);
						file << "," << duration << "," << abserror << std::endl;
					}
					std::cout << "done." << std::endl;
				}
				break;

			case AdvancedRules::algorithm::a_adaptiveGaussKronrod:
				gsl_integration_workspace * workspace;
				for (Functions::integrableFunction &function : functions.functions) {
					std::cout << "Calculating " << function.name << "... "
					      << std::flush;
					start = std::clock();
					value = AdvancedRules::adaptiveGaussKronrod(error_code,
					      function.f, function.a, function.b, error,
					      max_subdivisions, key, &abserror);
					duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
					file << "," << std::defaultfloat << function.type << ","
					      << function.name << " from " << function.a << " to "
					      << function.b << std::fixed << value;
					if (checkError(error_code)) {
						file << "," << error_code << "," << "," << duration
						      << std::endl;
						error_code = "";
					} else {
						file << "," << value;
						file << "," << fabs(value - function.value);
						file << "," << duration << "," << abserror << std::endl;
					}
					std::cout << "done." << std::endl;
				}

				break;

			case AdvancedRules::algorithm::a_adaptiveGaussKronrodSingular:
				for (Functions::integrableFunction &function : functions.functions) {
					std::cout << "Calculating " << function.name << "... "
					      << std::flush;
					start = std::clock();
					value = AdvancedRules::adaptiveGaussKronrodSingular(error_code,
					      function.f, function.a, function.b, error,
					      max_subdivisions, &abserror);
					duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
					file << "," << std::defaultfloat << function.type << ","
					      << function.name << " from " << function.a << " to "
					      << function.b << std::fixed << value;
					if (checkError(error_code)) {
						file << "," << error_code << "," << "," << duration
						      << std::endl;
						error_code = "";
					} else {
						file << "," << value;
						file << "," << fabs(value - function.value);
						file << "," << duration << "," << abserror << std::endl;
					}
					std::cout << "done." << std::endl;
				}
				break;

			case AdvancedRules::algorithm::a_adaptiveGaussKronrodKnownSingular:
				for (Functions::integrableFunction &function : functions.functions) {
					std::cout << "Calculating " << function.name << "... "
					      << std::flush;
					start = std::clock();
					value = AdvancedRules::adaptiveGaussKronrodKnownSingular(
					      error_code, function.f, function.a, function.b, error,
					      max_subdivisions, function.singularity, &abserror);
					duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
					file << "," << std::defaultfloat << function.type << ","
					      << function.name << " from " << function.a << " to "
					      << function.b << std::fixed << value;
					if (checkError(error_code)) {
						file << "," << error_code << "," << "," << duration
						      << std::endl;
						error_code = "";
					} else {
						file << "," << value;
						file << "," << fabs(value - function.value);
						file << "," << duration << "," << abserror << std::endl;
					}
					singularPoints.clear();
					std::cout << "done." << std::endl;
				}
				break;

			case AdvancedRules::algorithm::a_adaptiveGaussKronrodParallel:
				for (Functions::integrableFunction &function : functions.functions) {
					std::cout << "Calculating " << function.name << "... "
					      << std::flush;
					start = std::clock();
					value = AdvancedRules::adaptiveGaussKronrodParallel(error_code,
					      function.f, function.a, function.b, error,
					      max_subdivisions, key, threads, &abserror);
					duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
					file << "," << std::defaultfloat << function.type << ","
					      << function.name << " from " << function.a << " to "
					      << function.b << std::fixed << value;
					if (checkError(error_code)) {
						file << "," << error_code << "," << "," << duration
						      << std::endl;
						error_code = "";
					} else {
						file << "," << value;
						file << "," << fabs(value - function.value);
						file << "," << duration << "," << abserror << std::endl;
					}
					std::cout << "done." << std::endl;
				}
				break;

			case AdvancedRules::algorithm::a_adaptiveGaussKronrodSingularParallel:
				for (Functions::integrableFunction &function : functions.functions) {
					std::cout << "Calculating " << function.name << "... "
					      << std::flush;
					start = std::clock();
					value = AdvancedRules::adaptiveGaussKronrodSingularParallel(
					      error_code, function.f, function.a, function.b, error,
					      max_subdivisions, threads, &abserror);
					duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
					file << "," << std::defaultfloat << function.type << ","
					      << function.name << " from " << function.a << " to "
					      << function.b << std::fixed << value;
					if (checkError(error_code)) {
						file << "," << error_code << "," << "," << duration
						      << std::endl;
						error_code = "";
					} else {
						file << "," << value;
						file << "," << fabs(value - function.value);
						file << "," << duration << "," << abserror << std::endl;
					}
					std::cout << "done." << std::endl;
				}
				break;

			case AdvancedRules::algorithm::a_adaptiveGaussKronrodKnownSingularParallel:
				for (Functions::integrableFunction &function : functions.functions) {
					std::cout << "Calculating " << function.name << "... "
					      << std::flush;
					start = std::clock();
					value = AdvancedRules::adaptiveGaussKronrodKnownSingularParallel(
					      error_code, function.f, function.a, function.b, error,
					      max_subdivisions, function.singularity, threads,
					      &abserror);
					duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
					file << "," << std::defaultfloat << function.type << ","
					      << function.name << " from " << function.a << " to "
					      << function.b << std::fixed << value;
					if (checkError(error_code)) {
						file << "," << error_code << "," << "," << duration
						      << std::endl;
						error_code = "";
					} else {
						file << "," << value;
						file << "," << fabs(value - function.value);
						file << "," << duration << "," << abserror << std::endl;
					}
					singularPoints.clear();
					std::cout << "done." << std::endl;
				}
				break;

			}

		}
		file.close();
	}
}

