/**
 * @file BoolesRule.cpp
 * @brief contains functions for calculations with Boole's Rule.
 * Boole's Rule approximates the function interval as a series of
 * 4th degree polynomials, and their area is summed.
 * @author Irene Crowell
 */
#include <math.h>
#include <thread>
#include <queue>
#include <mutex>
#include "FindVal.h"

namespace BoolesRule {
/**
 * A simple section to integrate
 */
struct interval {
	double a; //!<the left (starting) point
	double m; //!<the midpoint
	double b; //!<the right (ending) point
	double fa; //!<the value of the function at m
	double fm; //!<the value of the function at a
	double fb; //!<the value of the function at b
};
/**
 * An already integrated section divided in two
 */
struct leftRightInterval {
	interval left; //!<the left half of the interval
	interval right; //!<the right half of the interval
	double integrated; //!<the calculated integral of the interval, for estimating error in adaptive rules
};
/**
 * Integrates and divides an interval into two, resulting in a leftRightInterval (for adaptive rules)
 * @param currentInterval the interval to integrate and divide
 * @param f the function to integrate
 * @return the divided and integrated interval
 */
leftRightInterval getLeftRight(interval currentInterval, gsl_function f) {
	double width = currentInterval.b - currentInterval.a;
	double m1 = currentInterval.a + width / 4;
	double m2 = currentInterval.b - width / 4;
	double fm1 = findVal(f, m1, width);
	double fm2 = findVal(f, m2, width);
	double integrate = (width / 90)
	      * (7 * currentInterval.fa + 32 * fm1 + 12 * currentInterval.fm
	            + 32 * fm2 + 7 * currentInterval.fb);
	interval leftInterval = { currentInterval.a, m1, currentInterval.m,
	      currentInterval.fa, fm1, currentInterval.fm };
	interval rightInterval = { currentInterval.m, m2, currentInterval.b,
	      currentInterval.fm, fm2, currentInterval.fb };

	return {leftInterval, rightInterval, integrate};
}
/**
 * Calculates the numerical integral using Simpson's 3/8 rule
 * @param f	the function to integrate
 * @param a the left (starting) point of the integral section
 * @param b the right (ending) point of the integral section
 * @param subdivisions the number of subdivisions to use
 * @return the numerically integrated value
 */
double nonAdaptiveNonParallel(gsl_function f, double a, double b,
      int subdivisions) {
	double result = 0;
	double width = (b - a) / subdivisions;

	result += 7 * findVal(f, a, width);

	double x = a + width / 4;
	while (x < b - width) {
		result += 32 * findVal(f, x, width);
		x += width * 0.25;
		result += 12 * findVal(f, x, width);
		x += width * 0.25;
		result += 32 * findVal(f, x, width);
		x += width * 0.25;
		result += 14 * findVal(f, x, width);
		x += width * 0.25;
	}

	result += 32 * findVal(f, b - width * 0.75, width);
	result += 12 * findVal(f, b - width * 0.5, width);
	result += 32 * findVal(f, b - width * 0.25, width);
	result += 7 * findVal(f, b, width);

	return result * width / 90;
}
/**
 * For threading -- calculates a section of the integral using Boole's rule
 * @param f	the function to integrate
 * @param a the left (starting) point of the integral section
 * @param b the right (ending) point of the integral section
 * @param width the size of the interval (b-a)
 * @param threads the total number of threads being run
 * @param threadNum the identifying member of this thread (0 to threads-1)
 * @param result_mutex pointer to a mutex for writing to the result
 * @param [out] result the sum of the numerical integral sections
 */
void nonAdaptiveThread(gsl_function f, double a, double b, double width,
      int threads, int threadNum, std::mutex *result_mutex, double *result) {
	double x = a + width * (threadNum + 0.25);
	double integrate;
	while (x < b - width) {
		integrate += 32 * findVal(f, x, width);
		x += width * 0.25;
		integrate += 12 * findVal(f, x, width);
		x += width * 0.25;
		integrate += 32 * findVal(f, x, width);
		x += width * 0.25;
		integrate += 14 * findVal(f, x, width);
		x += width * (threads - 0.75);
	}
	result_mutex->lock();
	(*result) += integrate;
	result_mutex->unlock();
}
/**
 * Calculates using parallel sections the integral using Boole's Rule
 * @param f	the function to integrate
 * @param a the left (starting) point of the integral section
 * @param b the right (ending) point of the integral section
 * @param subdivisions the number of subdivisions to use
 * @param num_threads the number of parallel threads to run
 * @return the numerically integrated value
 */
double nonAdaptiveParallel(gsl_function f, double a, double b, int subdivisions,
      int num_threads) {
	std::thread threads[num_threads];
	std::mutex result_mutex;

	double result = 0;
	double width = (b - a) / subdivisions;

	result += findVal(f, a, width);

	for (int i = 0; i < num_threads; i++) {
		threads[i] = std::thread(nonAdaptiveThread, f, a, b, width, num_threads,
		      i, &result_mutex, &result);
	}
	for (int i = 0; i < num_threads; i++) {
		threads[i].join();
	}

	result += 32 * findVal(f, b - width * 0.75, width);
	result += 12 * findVal(f, b - width * 0.5, width);
	result += 32 * findVal(f, b - width * 0.25, width);
	result += 7 * findVal(f, b, width);

	return result * width / 90;
}

/**
 * Calculates the numerical integral using an adaptive Boole's rule.
 * The adaptive rule divides each section in two until the error goal is met.
 * @param f	the function to integrate
 * @param a the left (starting) point of the integral section
 * @param b the right (ending) point of the integral section
 * @param error the error goal
 * @param max_subdivisions the maximum number of subdivisions to use
 * @param max_time the time limit for the calculation
 * @param [out] subdivisions the number of subdivisions used
 * @return the numerically integrated value
 */
double adaptiveNonParallel(gsl_function f, double a, double b, double error,
      int max_subdivisions, int max_time, int *subdivisions) {
	(*subdivisions) = 1;
	std::clock_t start = std::clock();
	std::queue<leftRightInterval> intervals;
	interval whole = { a, (a + b) / 2, b, findVal(f, a, b - a), findVal(f,
	      (a + b) / 2, b - a), findVal(f, b, b - a) };
	intervals.push(getLeftRight(whole, f));
	double result = 0;
	bool subdivisions_exceeded = false;
	bool time_exceeded = false;
	while (!intervals.empty()) {
		if (((std::clock() - start) / (double) CLOCKS_PER_SEC) > max_time)
			time_exceeded = true;
		leftRightInterval currentInterval = intervals.front();
		if (intervals.size() > max_subdivisions)
			subdivisions_exceeded = true;
		intervals.pop();
		double width = currentInterval.right.b - currentInterval.left.a;
		leftRightInterval left = getLeftRight(currentInterval.left, f);
		leftRightInterval right = getLeftRight(currentInterval.right, f);

		if (fabs(left.integrated + right.integrated - currentInterval.integrated)
		      < (255 * width * error) //Richardson Extrapolation for error ((4^4)-1)
		|| subdivisions_exceeded || time_exceeded) {
			result += left.integrated + right.integrated;
		} else {
			intervals.push(left);
			intervals.push(right);
			(*subdivisions)++;
		}
	}
	return result;

}
/**
 * For threading -- Calculates a section of the  numerical integral using an
 * adaptive Boole's rule.
 * The adaptive rule divides each section in two until the error goal is met.
 * @param f	the function to integrate
 * @param error the error goal
 * @param max_subdivisions the maximum number of subdivisions to use
 * @param max_time the time limit for the calculation
 * @param intervals_mutex pointer to a mutex for changing the interval queue
 * @param intervals pointer to the queue of intervals to be integrated
 * @param result_mutex pointer to a mutex for writing to the result
 * @param [out] result the sum of the integrals of the sections
 * @param [out] subdivisions the number of subdivisions used
 */
void adaptiveThread(gsl_function f, double error, int max_subdivisions,
      int max_time, std::mutex *intervals_mutex,
      std::queue<leftRightInterval> *intervals, std::mutex *result_mutex,
      double *result, int *subdivisions) {

	std::clock_t start = std::clock();
	bool subdivisions_exceeded = false;
	bool time_exceeded = false;
	intervals_mutex->lock();
	while (!intervals->empty()) {
		if (((std::clock() - start) / (double) CLOCKS_PER_SEC) > max_time)
			time_exceeded = true;
		leftRightInterval currentInterval = intervals->front();
		if (intervals->size() > max_subdivisions)
			subdivisions_exceeded = true;
		intervals->pop();
		intervals_mutex->unlock();
		double width = currentInterval.right.b - currentInterval.left.a;
		leftRightInterval left = getLeftRight(currentInterval.left, f);
		leftRightInterval right = getLeftRight(currentInterval.right, f);
		if (fabs(left.integrated + right.integrated - currentInterval.integrated)
		      < (255 * width * error) //Richardson Extrapolation for error ((4^4)-1)
		|| subdivisions_exceeded || time_exceeded) {
			result_mutex->lock();
			(*result) += left.integrated + right.integrated;
			result_mutex->unlock();
		} else { //divide interval in 2 to be integrated again
			intervals_mutex->lock();
			intervals->push(left);
			intervals->push(right);
			(*subdivisions)++;
			intervals_mutex->unlock();
		}
		intervals_mutex->lock();
	}
	intervals_mutex->unlock();
}
/**
 * Calculates using parallel threads the numerical integral using an adaptive Boole's rule.
 * The adaptive rule divides each section in two until the error goal is met.
 * @param f	the function to integrate
 * @param a the left (starting) point of the integral section
 * @param b the right (ending) point of the integral section
 * @param num_threads the number of parallel threads to run
 * @param error the error goal
 * @param max_subdivisions the maximum number of subdivisions to use
 * @param max_time the time limit for the calculation
 * @param [out] subdivisions the number of subdivisions used
 * @return the numerically integrated value
 */
double adaptiveParallel(gsl_function f, double a, double b, int num_threads,
      double error, int max_subdivisions, int max_time, int *subdivisions) {
	(*subdivisions) = 1;
	std::queue<leftRightInterval> intervals;
	interval whole = { a, (a + b) / 2, b, findVal(f, a, b - a), findVal(f,
	      (a + b) / 2, b - a), findVal(f, b, b - a) };
	intervals.push(getLeftRight(whole, f));
	double result = 0;
	std::thread threads[num_threads];
	std::mutex result_mutex;
	std::mutex intervals_mutex;

	for (int i = 0; i < num_threads; i++) {
		threads[i] = std::thread(adaptiveThread, f, error, max_subdivisions,
		      max_time, &intervals_mutex, &intervals, &result_mutex, &result,
		      subdivisions);
	}
	for (int i = 0; i < num_threads; i++) {
		threads[i].join();
	}
	return result;

}
}
