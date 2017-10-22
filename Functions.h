/**
 * @file Functions.h
 * @brief Contains mathematical functions to be integrated
 * @author Irene Crowell
 */
#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <math.h>
#include <array>
#include <string>
#include <vector>
/**
 * Provides an array of 24 integrableFunctions to be integrated
 */
class Functions {
public:
	Functions();
	/**
	 * Defines an area of a function to be integrated along with
	 * information about the function
	 */
	struct integrableFunction {
		gsl_function f; //!<the function to be integrated
		double a; //!<the left (starting) point
		double b; //!<the right (ending) point
		double value; //!<the symbolically calculated integral
		std::string name; //!<a description of the function, ex: sin(x)
		std::string type; //!<either simple, discontinuous, or badly behaved
		double singularity; //!<a singular point, if one exists
	};

	std::array<integrableFunction, 23> functions;
};

#endif /* FUNCTIONS_H_ */
