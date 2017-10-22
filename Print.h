/**
 * @file Print.h
 * @brief Contains function prototypes for the test printer
 * @author Irene Crowell
 */

#ifndef PRINT_H_
#define PRINT_H_
#include <iostream>
#include <fstream>
#include <math.h>
#include <ctime>
#include <sstream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include "Functions.h"
#include "NewtonCotesRules/RuleHeaders.h"
#include "AdvancedRules/AdvancedRules.h"

/**
 * Prints the outputs of the Non Adaptive Non Parallel version of the Newton-Cotes rules to csv files.
 * Two tests are done with faster and slower parameters, printing to nonAdaptiveNonParallelFast.csv
 * and nonAdaptiveNonParallelAccurate.csv
 * @param subdivisionsFast The number of subdivisions for the faster test
 * @param subdivisionsSlow The number of subdivisions for the slower (more accurate) test
 */
void printNonAdaptiveNonParallel(int subdivisionsFast, int subdivisionsSlow);
/**
 * Prints the outputs of the Non Adaptive Parallel version of the Newton-Cotes rules to csv files.
 * Two tests are done with faster and slower parameters, printing to nonAdaptiveParallelFast.csv
 * and nonAdaptiveParallelAccurate.csv
 * @param subdivisionsFast The number of subdivisions for the faster test
 * @param subdivisionsSlow The number of subdivisions for the slower (more accurate) test
 * @param threads the number of threads to run in parallel
 */
void printNonAdaptiveParallel(int subdivisionsFast, int subdivisionsSlow,
      int threads);
/**
 * Prints the outputs of the AAdaptive Non Parallel version of the Newton-Cotes rules to csv files.
 * Two tests are done with faster and slower parameters, printing to adaptiveNonParallelFast.csv
 * and adaptiveNonParallelAccurate.csv
 * @param max_subdivisions the maximum subdivisions to be used for the test
 * @param timeFast the time limit for the faster test
 * @param timeSlow the time limit for the slower (more accurate) test
 * @param errorFast the error goal for the faster test
 * @param errorSlow the error goal for the slower (more accurate) test
 */
void printAdaptiveNonParallel(int max_subdivisions, int timeFast, int timeSlow,
      double errorFast, double errorSlow);
/**
 * Prints the outputs of the Adaptive Parallel version of the Newton-Cotes rules to csv files.
 * Two tests are done with faster and slower parameters, printing to adaptiveParallelFast.csv
 * and adaptiveParallelAccurate.csv
 * @param max_subdivisions the maximum subdivisions to be used for the test
 * @param timeFast the time limit for the faster test
 * @param timeSlow the time limit for the slower (more accurate) test
 * @param errorFast the error goal for the faster test
 * @param errorSlow the error goal for the slower (more accurate) test
 * @param threads the number of threads to run in parallel
 */
void printAdaptiveParallel(int max_subdivisions, int timeFast, int timeSlow,
      double errorFast, double errorSlow, int threads);
/**
 * Prints the outputs of the Adaptive Parallel version of the Newton-Cotes rules to csv files.
 * Two tests are done with faster and slower parameters, printing to adaptiveParallelFast.csv
 * and adaptiveParallelAccurate.csv
 * @param max_subdivisions
 * @param timeFast the time limit for the faster test
 * @param timeSlow the time limit for the slower (more accurate) test
 * @param errorFast the error goal for the faster test
 * @param errorSlow the error goal for the slower (more accurate) test
 * @param threads the number of threads to run in parallel
 * @param max_subdivisions the maximum subdivisions to be used for the test
 * @param errorFast the error goal for the faster test
 * @param errorSlow the error goal for the slower (more accurate) test
 * @param pointsFast the number of points (Gauss-Legendre) for the faster test
 * @param pointsSlow the number of points (Gauss-Legendre) for the slower (more accurate) test
 * @param keyFast the key (Gauss-Kronrod) for the faster test
 * @param keySlow the key (Gauss-Kronrod) for the slower (more accurate) test
 * @param threads the number of threads to use for the parallel tests
 */
void printAdvanced(int max_subdivisions, double errorFast, double errorSlow,
      int pointsFast, int pointsSlow, int keyFast, int keySlow, int threads);

#endif /* PRINT_H_ */
