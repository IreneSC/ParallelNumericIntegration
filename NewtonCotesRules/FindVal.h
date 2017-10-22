/**
 * @file FindVal.h
 * @brief Allows for calculation of a gsl_function at a point
 * @author Irene Crowell
 */
#ifndef FINDVAL_H_
#define FINDVAL_H_
#include <gsl/gsl_math.h>
#include <iostream>
/**
 * Calculates a function's value at a point.
 * If the point is a singularity, it averages points nearby (relative to the width).
 * @param f the function to find the value of
 * @param x the point to find the value at
 * @param width the "width" to consider as close to the point.
 * @return the value of f at x
 */
inline double findVal(gsl_function f, double x, double width) {
	double* q = 0;
	double val = f.function(x, q);
	if (isnan(val) || isinf(val)) {
		double fa = f.function(x - 0.001 * width, q);
		double fb = f.function(x + 0.001 * width, q);
		if (isnan(fa) || isinf(fa)) {
			double fm = f.function(x + 0.0005 * width, q);
			val = 2 * fm - fb; //same as fm-(fb-fm) -> linear guess at value of f(x)
		} else if (isnan(fb) || isinf(fb)) {
			double fm = f.function(x - 0.0005 * width, q);
			val = 2 * fm - fa; //same as fm-(fa-fm) -> linear guess at value of f(x)
		} else {
			val = (fa + fb) / 2.0;
		}
	}
	return val;
}

#endif /* FINDVAL_H_ */
