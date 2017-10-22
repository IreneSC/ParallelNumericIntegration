/**
 * @file RuleHeaders.h
 * @brief Contains function prototypes for the Newton Cotes Rules
 * @author Irene Crowell
 */
#ifndef RULEHEADERS_H_
#define RULEHEADERS_H_

namespace MidpointRule {
double nonAdaptiveNonParallel(gsl_function f, double a, double b,
      int subdivisions);
double nonAdaptiveParallel(gsl_function f, double a, double b, int subdivisions,
      int num_threads);
double adaptiveNonParallel(gsl_function f, double a, double b, double error,
      int max_subdivisions, int max_time, int *subdivisions);
double adaptiveParallel(gsl_function f, double a, double b, int num_threads,
      double error, int max_subdivisions, int max_time, int *subdivisions);

}
namespace TrapezoidRule {
double nonAdaptiveNonParallel(gsl_function f, double a, double b,
      int subdivisions);
double nonAdaptiveParallel(gsl_function f, double a, double b, int subdivisions,
      int num_threads);
double adaptiveNonParallel(gsl_function f, double a, double b, double error,
      int max_subdivisions, int max_time, int *subdivisions);
double adaptiveParallel(gsl_function f, double a, double b, int num_threads,
      double error, int max_subdivisions, int max_time, int *subdivisions);
}

namespace SimpsonRule {
double nonAdaptiveNonParallel(gsl_function f, double a, double b,
      int subdivisions);
double nonAdaptiveParallel(gsl_function f, double a, double b, int subdivisions,
      int num_threads);
double adaptiveNonParallel(gsl_function f, double a, double b, double error,
      int max_subdivisions, int max_time, int *subdivisions);
double adaptiveParallel(gsl_function f, double a, double b, int num_threads,
      double error, int max_subdivisions, int max_time, int *subdivisions);
}
namespace Simpson38Rule {
double nonAdaptiveNonParallel(gsl_function f, double a, double b,
      int subdivisions);
double nonAdaptiveParallel(gsl_function f, double a, double b, int subdivisions,
      int num_threads);
double adaptiveNonParallel(gsl_function f, double a, double b, double error,
      int max_subdivisions, int max_time, int *subdivisions);
double adaptiveParallel(gsl_function f, double a, double b, int num_threads,
      double error, int max_subdivisions, int max_time, int *subdivisions);
}
namespace BoolesRule {
double nonAdaptiveNonParallel(gsl_function f, double a, double b,
      int subdivisions);
double nonAdaptiveParallel(gsl_function f, double a, double b, int subdivisions,
      int num_threads);
double adaptiveNonParallel(gsl_function f, double a, double b, double error,
      int max_subdivisions, int max_time, int *subdivisions);
double adaptiveParallel(gsl_function f, double a, double b, int num_threads,
      double error, int max_subdivisions, int max_time, int *subdivisions);
}

#endif /* RULEHEADERS_H_ */
