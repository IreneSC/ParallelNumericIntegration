/**
 * @file AdvancedRules.cpp
 * @brief Contains functions to implement the GSL rules
 * @author Irene Crowell
 */
#include "AdvancedRules.h"
namespace AdvancedRules {

/**
 * Calculates the numerical integral using a Guass-Legendre Rule with a fixed number of points.
 * @ingroup AdvancedRules
 * @param f the function to integrate
 * @param a the left (starting) point of the integral
 * @param b the right (ending) point of the integral
 * @param points the number of points
 * @return the numerically integrated value
 */
double gaussLegendreFixed(gsl_function f, double a, double b, int points) {
	gsl_integration_glfixed_table table;
	table = *gsl_integration_glfixed_table_alloc(points);
	return gsl_integration_glfixed(&f, a, b, &table);
}

/**
 * For threading -- calculates a section of the integral using a Guass-Legendre
 * Rule with a fixed number of points.
 * @param f	the function to integrate
 * @param a the left (starting) point of the integral section
 * @param b the right (ending) point of the integral section
 * @param points the number of points
 * @param result_mutex pointer to a mutex for writing to the result
 * @param [out] result the sum of the numerical integral sections
 */
void gaussLegendreFixedThread(gsl_function f, double a, double b, int points,
      std::mutex *result_mutex, double *result) {
	gsl_integration_glfixed_table table;
	table = *gsl_integration_glfixed_table_alloc(points);
	double integrate = gsl_integration_glfixed(&f, a, b, &table);
	result_mutex->lock();
	(*result) += integrate;
	result_mutex->unlock();
}

/**
 * Calculates using parallel sections the numerical integral using a Guass-Legendre Rule with
 * a fixed number of points.
 * @param f the function to integrate
 * @param a the left (starting) point of the integral
 * @param b the right (ending) point of the integral
 * @param points the number of points
 * @param num_threads the number of parallel threads to run
 * @return the numerically integrated value
 */
double gaussLegendreFixedParallel(gsl_function f, double a, double b,
      int points, int num_threads) {
	std::thread threads[num_threads];
	std::mutex result_mutex;
	double result = 0;
	double width = b - a;

	for (int i = 0; i < num_threads; i++) {
		threads[i] = std::thread(gaussLegendreFixedThread, f,
		      a + (i * 0.25 * width), a + ((i + 1) * 0.25 * width), points,
		      &result_mutex, &result);
	}
	for (int i = 0; i < num_threads; i++) {
		threads[i].join();
	}
	return result;
}

/**
 * Calculates  the numerical integral using a Non-Adaptive Guass-Kronrod Rule with an error bound
 * @param [out] error_code pointer to store an error
 * @param f the function to integrate
 * @param a the left (starting) point of the integral
 * @param b the right (ending) point of the integral
 * @param error the error goal
 * @param [out] abserror the maximum error achieved
 * @return the numerically integrated value
 */
double nonAdaptiveGaussKronrod(const char * error_code, gsl_function f,
      double a, double b, double error, double *abserror) {
	double result;
	size_t evals;
	int status = gsl_integration_qng(&f, a, b, error, error, &result, abserror,
	      &evals);
	if (status)
		error_code = gsl_strerror(status);
	return result;
}

/**
 * For threading -- calculates a section of the integral using a Non-Adaptive Guass-Kronrod
 * Rule with an error bound
 * @param [out] error_code pointer to store an error
 * @param f the function to integrate
 * @param a the left (starting) point of the integral
 * @param b the right (ending) point of the integral
 * @param error the error goal
 * @param [out] abserror the maximum error achieved
 * @param result_mutex pointer to a mutex for writing to the result
 * @param [out] result the sum of the numerical integral sections
 */
void nonAdaptiveGaussKronrodThread(const char * error_code, gsl_function f,
      double a, double b, double error, std::mutex*result_mutex, double *result,
      double *abserror) {
	double integrate;
	double errorabs;
	size_t evals;
	int status = gsl_integration_qng(&f, a, b, error, error, &integrate,
	      &errorabs, &evals);
	result_mutex->lock();
	if (status)
		error_code = gsl_strerror(status);
	(*result) += integrate;
	(*abserror) += pow(errorabs, 2);
	result_mutex->unlock();
}

/**
 * Calculates using parallel sections the numerical integral using a Guass-Kronrod Rule
 * with an error bound
 * @param [out] error_code pointer to store an error
 * @param f the function to integrate
 * @param a the left (starting) point of the integral
 * @param b the right (ending) point of the integral
 * @param error the error goal
 * @param num_threads the number of parallel threads to run
 * @param [out] abserror the maximum error achieved
 * @return the numerically integrated value
 */
double nonAdaptiveGaussKronrodParallel(const char * error_code, gsl_function f,
      double a, double b, double error, int num_threads, double *abserror) {
	std::thread threads[num_threads];
	std::mutex result_mutex;
	double result = 0;
	double width = b - a;

	for (int i = 0; i < num_threads; i++) {
		threads[i] = std::thread(nonAdaptiveGaussKronrodThread, error_code, f,
		      a + (i * 0.25 * width), a + ((i + 1) * 0.25 * width), error,
		      &result_mutex, &result, abserror);
	}
	for (int i = 0; i < num_threads; i++) {
		threads[i].join();
	}
	(*abserror) = sqrt(*abserror);
	return result;
}
/**
 * Calculates the numerical integral using an adaptive Guass-Kronrod Rule
 * with an error bound and specified key
 * @param [out] error_code pointer to store an error
 * @param f the function to integrate
 * @param a the left (starting) point of the integral
 * @param b the right (ending) point of the integral
 * @param error the error goal
 * @param max_subdivisions the maximum subdivisions to use
 * @param key the "key" for the GSL Gauss-Kronrod Rule
 * @param [out] abserror the maximum error achieved
 * @return the numerically integrated value
 */
double adaptiveGaussKronrod(const char * error_code, gsl_function f, double a,
      double b, double error, int max_subdivisions, int key, double *abserror) {

	gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(
	      (unsigned long int) max_subdivisions);
	double result;
	int status = gsl_integration_qag(&f, a, b, error, error,
	      (size_t) max_subdivisions, key, workspace, &result, abserror);
	if (status)
		error_code = gsl_strerror(status);
	gsl_integration_workspace_free(workspace);
	return result;
}
/**
 * For threading -- calculates a section of the integral using an adaptive Guass-Kronrod Rule
 * with an error bound and specified key
 * @param [out] error_code pointer to store an error
 * @param f the function to integrate
 * @param a the left (starting) point of the integral
 * @param b the right (ending) point of the integral
 * @param error the error goal
 * @param max_subdivisions the maximum subdivisions to use
 * @param key the "key" for the GSL Gauss-Kronrod Rule
 * @param result_mutex pointer to a mutex for writing to the result
 * @param [out] result the sum of the numerical integral sections
 * @param [out] abserror the maximum error achieved
 */
void adaptiveGaussKronrodThread(const char * error_code, gsl_function f,
      double a, double b, double error, int max_subdivisions, int key,
      std::mutex*result_mutex, double *result, double *abserror) {
	gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(
	      (unsigned long int) max_subdivisions);
	double integrate;
	double errorabs;
	int status = gsl_integration_qag(&f, a, b, error, error,
	      (size_t) max_subdivisions, key, workspace, &integrate, &errorabs);
	result_mutex->lock();
	if (status)
		error_code = gsl_strerror(status);
	(*result) += integrate;
	(*abserror) += pow(errorabs, 2);
	result_mutex->unlock();
	gsl_integration_workspace_free(workspace);
}
/**
 * Calculates using parallel sections the numerical integral using an adaptive
 * Guass-Kronrod Rule with an error bound.
 * @param [out] error_code pointer to store an error
 * @param f the function to integrate
 * @param a the left (starting) point of the integral
 * @param b the right (ending) point of the integral
 * @param error the error goal
 * @param max_subdivisions the maximum subdivisions to use
 * @param key the "key" for the GSL Gauss-Kronrod Rule
 * @param num_threads the number of parallel threads to run
 * @param [out] abserror the maximum error achieved
 * @return the numerically integrated value
 */
double adaptiveGaussKronrodParallel(const char * error_code, gsl_function f,
      double a, double b, double error, int max_subdivisions, int key,
      int num_threads, double *abserror) {
	std::thread threads[num_threads];
	std::mutex result_mutex;
	double result = 0;
	double width = b - a;

	for (int i = 0; i < num_threads; i++) {
		threads[i] = std::thread(adaptiveGaussKronrodThread, error_code, f,
		      a + (i * 0.25 * width), a + ((i + 1) * 0.25 * width), error,
		      (max_subdivisions / num_threads), key, &result_mutex, &result,
		      abserror);
	}
	for (int i = 0; i < num_threads; i++) {
		threads[i].join();
	}
	(*abserror) = sqrt(*abserror);
	return result;
}
/**
 * Calculates nthe numerical integral using an adaptive Guass-Kronrod Rule
 * with an error bound, for functions with singularities.
 * @param [out] error_code pointer to store an error
 * @param f the function to integrate
 * @param a the left (starting) point of the integral
 * @param b the right (ending) point of the integral
 * @param error the error goal
 * @param max_subdivisions the maximum subdivisions to use
 * @param [out] abserror the maximum error achieved
 * @return the numerically integrated value
 */
double adaptiveGaussKronrodSingular(const char * error_code, gsl_function f,
      double a, double b, double error, int max_subdivisions,
      double *abserror) {
	gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(
	      (unsigned long int) max_subdivisions);
	double result;
	int status = gsl_integration_qags(&f, a, b, error, error,
	      (size_t) max_subdivisions, workspace, &result, abserror);
	if (status)
		error_code = gsl_strerror(status);
	gsl_integration_workspace_free(workspace);
	return result;
}
/**
 * For threading -- calculates a section of the integral using an adaptive
 * Guass-Kronrod Rule with an error bound, for functions with singularities.
 * @param [out] error_code pointer to store an error
 * @param f the function to integrate
 * @param a the left (starting) point of the integral
 * @param b the right (ending) point of the integral
 * @param error the error goal
 * @param max_subdivisions the maximum subdivisions to use
 * @param result_mutex pointer to a mutex for writing to the result
 * @param [out] result the sum of the numerical integral sections
 * @param [out] abserror the maximum error achieved
 */
void adaptiveGaussKronrodSingularThread(const char * error_code, gsl_function f,
      double a, double b, double error, int max_subdivisions,
      std::mutex*result_mutex, double *result, double *abserror) {
	gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(
	      (unsigned long int) max_subdivisions);
	double integrate;
	double errorabs;
	int status = gsl_integration_qags(&f, a, b, error, error,
	      (size_t) max_subdivisions, workspace, &integrate, &errorabs);
	result_mutex->lock();
	if (status)
		error_code = gsl_strerror(status);
	(*result) += integrate;
	(*abserror) += pow(errorabs, 2);
	result_mutex->unlock();
	gsl_integration_workspace_free(workspace);
}
/**
 * Calculates using parallel sections the numerical integral using using an adaptive
 * Guass-Kronrod Rule with an error bound, for functions with singularities.
 * @param [out] error_code pointer to store an error
 * @param f the function to integrate
 * @param a the left (starting) point of the integral
 * @param b the right (ending) point of the integral
 * @param error the error goal
 * @param max_subdivisions the maximum subdivisions to use
 * @param num_threads the number of parallel threads to run
 * @param [out] abserror the maximum error achieved
 * @return the numerically integrated value
 */
double adaptiveGaussKronrodSingularParallel(const char * error_code,
      gsl_function f, double a, double b, double error, int max_subdivisions,
      int num_threads, double *abserror) {
	std::thread threads[num_threads];
	std::mutex result_mutex;
	double result = 0;
	double width = b - a;

	for (int i = 0; i < num_threads; i++) {
		threads[i] = std::thread(adaptiveGaussKronrodSingularThread, error_code,
		      f, a + (i * 0.25 * width), a + ((i + 1) * 0.25 * width), error,
		      (max_subdivisions / num_threads), &result_mutex, &result, abserror);
	}
	for (int i = 0; i < num_threads; i++) {
		threads[i].join();
	}
	(*abserror) = sqrt(*abserror);
	return result;
}
/**
 * Calculates the numerical integral using an adaptive Guass-Kronrod Rule
 * with an error bound, for functions with a known singularity.
 * @param [out] error_code pointer to store an error
 * @param f the function to integrate
 * @param a the left (starting) point of the integral
 * @param b the right (ending) point of the integral
 * @param error the error goal
 * @param max_subdivisions the maximum subdivisions to use
 * @param singularity the singular point
 * @param [out] abserror the maximum error achieved
 * @return the numerically integrated value
 */
double adaptiveGaussKronrodKnownSingular(const char * error_code,
      gsl_function f, double a, double b, double error, int max_subdivisions,
      double singularity, double *abserror) {
	gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(
	      (unsigned long int) max_subdivisions);
	std::vector<double> points;
	points.push_back(a);
	if (singularity != 0)
		points.push_back(singularity);
	points.push_back(b);
	double result;
	int status = gsl_integration_qagp(&f, points.data(), points.size(), error,
	      error, (size_t) max_subdivisions, workspace, &result, abserror);
	if (status)
		error_code = gsl_strerror(status);
	gsl_integration_workspace_free(workspace);
	return result;
}
/**
 * For threading -- calculates a section of the integral using an adaptive
 * Guass-Kronrod Rule with an error bound, for functions with known singularities.
 * @param [out] error_code pointer to store an error
 * @param f the function to integrate
 * @param a the left (starting) point of the integral
 * @param b the right (ending) point of the integral
 * @param error the error goal
 * @param max_subdivisions the maximum subdivisions to use
 * @param points the end-points and singularities of the region
 * @param result_mutex pointer to a mutex for writing to the result
 * @param [out] result the sum of the numerical integral sections
 * @param [out] abserror the maximum error achieved
 */
void adaptiveGaussKronrodKnownSingularThread(const char * error_code,
      gsl_function f, double a, double b, double error, int max_subdivisions,
      std::vector<double> points, std::mutex*result_mutex, double *result,
      double *abserror) {
	gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(
	      (unsigned long int) max_subdivisions);
	double integrate;
	double errorabs;
	int status = gsl_integration_qagp(&f, points.data(), points.size(), error,
	      error, (size_t) max_subdivisions, workspace, &integrate, &errorabs);
	result_mutex->lock();
	if (status)
		error_code = gsl_strerror(status);
	(*result) += integrate;
	(*abserror) += pow(errorabs, 2);
	result_mutex->unlock();
	gsl_integration_workspace_free(workspace);
}
/**
 * Calculates using parallel sections the numerical integral using using an adaptive
 * Guass-Kronrod Rule with an error bound, for functions with a known singularity.
 * @param [out] error_code pointer to store an error
 * @param f the function to integrate
 * @param a the left (starting) point of the integral
 * @param b the right (ending) point of the integral
 * @param error the error goal
 * @param max_subdivisions the maximum subdivisions to use
 * @param singularity the singular point
 * @param num_threads the number of parallel threads to run
 * @param [out] abserror the maximum error achieved
 * @return the numerically integrated value
 */
double adaptiveGaussKronrodKnownSingularParallel(const char * error_code,
      gsl_function f, double a, double b, double error, int max_subdivisions,
      double singularity, int num_threads, double *abserror) {
	std::thread threads[num_threads];
	std::mutex result_mutex;
	double result = 0;
	double width = b - a;
	for (int i = 0; i < num_threads; i++) {
		double aNew = a + (i * 0.25 * width);
		double bNew = a + (i * 0.25 * width);
		std::vector<double> subPoints;
		subPoints.push_back(aNew);
		if (singularity != 0 && singularity > aNew && singularity < bNew)
			subPoints.push_back(singularity);
		subPoints.push_back(bNew);
		threads[i] = std::thread(adaptiveGaussKronrodKnownSingularThread,
		      error_code, f, aNew, bNew, error, (max_subdivisions / num_threads),
		      subPoints, &result_mutex, &result, abserror);
	}
	for (int i = 0; i < num_threads; i++) {
		threads[i].join();
	}
	(*abserror) = sqrt(*abserror);
	return result;
}
}
