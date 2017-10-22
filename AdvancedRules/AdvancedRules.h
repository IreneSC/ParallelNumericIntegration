/**
 * @file AdvancedRules.h
 * @brief Contains function prototypes for the GSL implentations
 * @author Irene Crowell
 */

#ifndef ADVANCEDRULES_ADVANCEDRULES_H_
#define ADVANCEDRULES_ADVANCEDRULES_H_
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <thread>
#include <mutex>
#include <vector>
namespace AdvancedRules {
/**
 * Allows for specification of which algorithm to call
 */
enum algorithm {
	a_gaussLegendreFixed,
	a_nonAdaptiveGaussKronrod,

	a_gaussLegendreFixedParallel,
	a_nonAdaptiveGaussKronrodParallel,

	a_adaptiveGaussKronrod,
	a_adaptiveGaussKronrodSingular,
	a_adaptiveGaussKronrodKnownSingular,

	a_adaptiveGaussKronrodParallel,
	a_adaptiveGaussKronrodSingularParallel,
	a_adaptiveGaussKronrodKnownSingularParallel,

	size //!< provides the length of the enum (for looping)
};

double gaussLegendreFixed(gsl_function f, double a, double b, int points);

double gaussLegendreFixedParallel(gsl_function f, double a, double b,
      int points, int num_threads);

double nonAdaptiveGaussKronrod(const char * error_code, gsl_function f,
      double a, double b, double error, double *abserror);
double nonAdaptiveGaussKronrodParallel(const char * error_code, gsl_function f,
      double a, double b, double error, int num_threads, double *abserror);
double adaptiveGaussKronrod(const char * error_code, gsl_function f, double a,
      double b, double error, int max_subdivisions, int key, double *abserror);
double adaptiveGaussKronrodParallel(const char * error_code, gsl_function f,
      double a, double b, double error, int max_subdivisions, int key,
      int num_threads, double *abserror);
double adaptiveGaussKronrodSingular(const char * error_code, gsl_function f,
      double a, double b, double error, int max_subdivisions, double *abserror);
double adaptiveGaussKronrodSingularParallel(const char * error_code,
      gsl_function f, double a, double b, double error, int max_subdivisions,
      int num_threads, double *abserror);
double adaptiveGaussKronrodKnownSingular(const char * error_code,
      gsl_function f, double a, double b, double error, int max_subdivisions,
      double singularity, double *abserror);
double adaptiveGaussKronrodKnownSingularParallel(const char * error_code,
      gsl_function f, double a, double b, double error, int max_subdivisions,
      double singularity, int num_threads, double *abserror);
}

#endif /* ADVANCEDRULES_ADVANCEDRULES_H_ */
