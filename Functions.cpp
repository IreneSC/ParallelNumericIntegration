/**
 * @file Functions.cpp
 * @brief Contains mathematical functions to be integrated
 * @author Irene Crowell
 */

#include "Functions.h"

double f_sin(double x, void * params) {
	return sin(x);
}
double f_square(double x, void * params) {
	return pow(x, 2);
}
double f_cube(double x, void * params) {
	return pow(x, 3);
}
double f_10000power(double x, void * params) {
	return pow(x, 10000);
}
double f_polynomial4(double x, void * params) {
	return 3 * pow(x, 4) + 4 * pow(x, 3) + 76 * pow(x, 2) + 58 * x + 4;
}
double f_inverse(double x, void * params) {
	return 1 / x;
}
double f_exp(double x, void * params) {
	return exp(x);
}
double f_sqrt(double x, void * params) {
	return sqrt(x);
}
//discontinuous
double f_ex(double x, void * params) {
	return (exp(x) - 1) / x;
}
double f_sqrt_abs_inv(double x, void * params) {
	return 1.0 / (sqrt(fabs(x)));
}
double f_sin_sqrt(double x, void * params) {
	return sin(x) * sqrt(1 - pow(x, 2));
}
double f_sin_x(double x, void * params) {
	return sin(x) / x;
}
double f_sqrt_abs(double x, void * params) {
	return sqrt(fabs(x - 0.7));
}
double f_rational(double x, void * params) {
	return (x - 2) * (x + 2) / (x - 2);
}
double f_floor(double x, void * params) {
	return floor(x);
}
double f_piecewise(double x, void * params) {
	if (x < 1)
		return 0.5 * x;
	if (x > 1)
		return 1.5 * x;
	return nan("");
}

//bad behavior
double f_ex_lnsin(double x, void * params) {
	return exp(x) * log(sin(x));
}
double f_ln(double x, void * params) {
	return log(x);
}
double f_ln_squared(double x, void * params) {
	return log(pow(x, 2));
}
double f_ln_inv(double x, void * params) {
	return log(1 / x);
}
double f_squared_sin_inv(double x, void * params) {
	return pow(x, 2) * sin(1 / x);
}
double f_tan(double x, void * params) {
	return tan(x);
}
double f_why(double x, void * params) {
	return sqrt(1 - pow(x, 4)) / pow(x, (1.0 - 1.0 / M_PI));
}
double f_why2(double x, void * params) {
	return exp(x) / pow(x, 1 / M_PI);
}
/**
 * Initializes the integrableFunctions
 */
Functions::Functions() {
	double alpha = 1.0;
	functions = {};
	functions = {};
	functions = {};

	functions[0].f.function = &f_sin;
	functions[0].f.params = &alpha;
	functions[0].a = 0.0;
	functions[0].b = 1.0;
	functions[0].name = "sin(x)";
	functions[0].type = "simple";
	functions[0].value = 1.0 - cos(1);

	functions[1].f.function = &f_square;
	functions[1].f.params = &alpha;
	functions[1].a = 0.0;
	functions[1].b = 1.0;
	functions[1].name = "x^2";
	functions[1].type = "simple";
	functions[1].value = 1 / 3.0;

	functions[2].f.function = &f_cube;
	functions[2].f.params = &alpha;
	functions[2].a = 0.0;
	functions[2].b = 1.0;
	functions[2].name = "x^3";
	functions[2].type = "simple";
	functions[2].value = 0.25;

	functions[3].f.function = &f_10000power;
	functions[3].f.params = &alpha;
	functions[3].a = 0.0;
	functions[3].b = 1.0;
	functions[3].name = "x^10000";
	functions[3].type = "simple";
	functions[3].value = 1.0 / 10001;

	functions[4].f.function = &f_polynomial4;
	functions[4].f.params = &alpha;
	functions[4].a = -20.0;
	functions[4].b = 0.0;
	functions[4].name = "3*x^4+4*x^3+76*x^2+58*x+4";
	functions[4].type = "simple";
	functions[4].value = 5853440.0 / 3;

	functions[5].f.function = &f_inverse;
	functions[5].f.params = &alpha;
	functions[5].a = 1.0;
	functions[5].b = 2.0;
	functions[5].name = "x^-1";
	functions[5].type = "simple";
	functions[5].value = log(2);

	functions[6].f.function = &f_exp;
	functions[6].f.params = &alpha;
	functions[6].a = 0.0;
	functions[6].b = 1.0;
	functions[6].name = "e^x";
	functions[6].type = "simple";
	functions[6].value = exp(1) - 1.0;

	functions[7].f.function = &f_sqrt;
	functions[7].f.params = &alpha;
	functions[7].a = 0.0;
	functions[7].b = 1.0;
	functions[7].name = "x^1/2";
	functions[7].type = "simple";
	functions[7].value = 2.0 / 3.0;

	functions[8].f.function = &f_ex;
	functions[8].f.params = &alpha;
	functions[8].a = 0.0;
	functions[8].b = 1.0;
	functions[8].name = "(e^x-1)/x (edge discontinuity)";
	functions[8].type = "discontinuous";
	functions[8].value = gsl_sf_expint_Ei(1) - M_EULER;

	functions[9].f.function = &f_sqrt_abs_inv;
	functions[9].f.params = &alpha;
	functions[9].a = -9.0;
	functions[9].b = 1.0;
	functions[9].name = "1/sqrt(abs(x)) (edge discontinuity)";
	functions[9].type = "discontinuous";
	functions[9].value = 8.0;

	functions[10].f.function = &f_sin_sqrt;
	functions[10].f.params = &alpha;
	functions[10].a = 0.0;
	functions[10].b = 1.0;
	functions[10].name =
	      "sin(x)*sqrt(1-x^2) (edge discontinuity) *note: answer only accurate to e-6";
	functions[10].type = "discontinuous";
	functions[10].value = 0.311736;

	functions[11].f.function = &f_sin_x;
	functions[11].f.params = &alpha;
	functions[11].a = 0.0;
	functions[11].b = 3.0;
	functions[11].name = "sin(x)/x (edge discontinuity)";
	functions[11].type = "discontinuous";
	functions[11].value = gsl_sf_Si(3.0);

	functions[12].f.function = &f_sqrt_abs;
	functions[12].f.params = &alpha;
	functions[12].a = 0.0;
	functions[12].b = 1.0;
	functions[12].name =
	      "sqrt(abs(x-0.7)) (middle discontinuity) *note: answer only accurate to e-6";
	functions[12].type = "discontinuous";
	functions[12].value = 0.499986;
	functions[12].singularity = 0.7;

	functions[13].f.function = &f_rational;
	functions[13].f.params = &alpha;
	functions[13].a = 1.0;
	functions[13].b = 3.0;
	functions[13].name = "(x+2)(x-2)/(x-2) (middle discontinuity)";
	functions[13].type = "discontinuous";
	functions[13].value = 8.0;
	functions[13].singularity = 2;

	functions[14].f.function = &f_floor;
	functions[14].f.params = &alpha;
	functions[14].a = 0.0;
	functions[14].b = 2.0;
	functions[14].name = "floor(x) (middle discontinuity)";
	functions[14].type = "discontinuous";
	functions[14].value = 1.0;
	functions[14].singularity = 1;

	functions[15].f.function = &f_piecewise;
	functions[15].f.params = &alpha;
	functions[15].a = 0.0;
	functions[15].b = 2.0;
	functions[15].name = "piecewise (middle discontinuity)";
	functions[15].type = "discontinuous";
	functions[15].value = 2.5;
	functions[15].singularity = 1;

	functions[16].f.function = &f_ex_lnsin;
	functions[16].f.params = &alpha;
	functions[16].a = 0.0;
	functions[16].b = M_PI;
	functions[16].name = "e^x*ln(sin(x)) *note: answer only accurate to e-4";
	functions[16].type = "badly behaved";
	functions[16].value = -20.8449;

	functions[17].f.function = &f_ln;
	functions[17].f.params = &alpha;
	functions[17].a = 0.0;
	functions[17].b = 1.0;
	functions[17].name = "ln(x)";
	functions[17].type = "badly behaved";
	functions[17].value = -1;

	functions[18].f.function = &f_ln_squared;
	functions[18].f.params = &alpha;
	functions[18].a = 8.0;
	functions[18].b = 9.0;
	functions[18].name = "ln(x^2)";
	functions[18].type = "badly behaved";
	functions[18].value = -2.0 * (1.0 + 8.0 * log(8) - 9.0 * log(9));

	functions[19].f.function = &f_ln_inv;
	functions[19].f.params = &alpha;
	functions[19].a = 0.0;
	functions[19].b = 1.0;
	functions[19].name = "ln(1/x)";
	functions[19].type = "badly behaved";
	functions[19].value = 1;

	functions[20].f.function = &f_squared_sin_inv;
	functions[20].f.params = &alpha;
	functions[20].a = 0.0;
	functions[20].b = 2.0 / M_PI;
	functions[20].name = "x^2*sin(1/x) (oscillatory)";
	functions[20].type = "badly behaved";
	functions[20].value = (16.0
	      + pow(M_PI, 2) * (-2.0 + M_PI * gsl_sf_Ci(M_PI / 2)))
	      / (6.0 * pow(M_PI, 3));

	functions[21].f.function = &f_why;
	functions[21].f.params = &alpha;
	functions[21].a = 0.0;
	functions[21].b = 1.0;
	functions[21].name = "sqrt(1-x^4)/x^(1-1/pi)";
	functions[21].type = "badly behaved";
	functions[21].value = (sqrt(M_PI) * gsl_sf_gamma(1.0 / (4.0 * M_PI)))
	      / (8.0 * gsl_sf_gamma((6.0 + 1.0 / M_PI) / 4.0));

	functions[22].f.function = &f_why2;
	functions[22].f.params = &alpha;
	functions[22].a = 0.0;
	functions[22].b = 1.0;
	functions[22].name = "e^x/(x^(1/pi)";
	functions[22].type = "badly behaved";
	functions[22].value = 2.303904211820843;

}

