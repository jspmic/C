#include <stdio.h>
#include <stdlib.h>

#define MAX_ITERATIONS 100

/*

 * @description: Newton-Cotes formulas are a group of formulas for numerical integration based on evaluating the integrand at equally spaced points.

 * @conventions: 
	 - In this file, `a` and `b` represents the lower and upper bound, respectively.
	 - The integrand(function to integrate) is called 'f'
	 - All integrands take a double as an argument, even if no decimal numbers are provided or needed(for generalization purposes and precision)

 * @reference: For more information, consult https://en.wikipedia.org/wiki/Newton%E2%80%93Cotes_formulas

*/


/*
 * @brief: Function that represents the function to integrate
 * @args: Takes a double as an argument to generalize all values of the input
 * @examples: f(x) = y, x is the input, y is the output
*/
typedef double (*function)(double);


/*

 * @brief: Function that represents a method to solve an integral
 * @args:
 - function: Function to integrate
 - An integer that represents the number of iterations to make when estimating the ssolution
 - A double that represents the lower bound
 - Another double that represents the upper bound

*/

typedef double (*method)(function, int, double, double);

// Functions representing all the solving methods implemented

double trapezoid(function, int, double, double);
double simpson1_3(function, int, double, double);
double simpson3_8(function, int, double, double);
double mid_point_rule(function, int, double, double);
double boole(function, int, double, double);

// Example functions to integrate

double f(double x){
	return x;
}

double f1(double x){
	return x*x;
}

double f2(double x){
	return x*x*x;
}

/*
 * @brief: Function that generate an example for each solving method function
 * @args:
	 - method: The function representing the method used
	 - method_name: The name of the method used
	 - f: function to integrate
	 - a: the lower bound
	 - b: the upper bound
*/

void example(method _method, const char* method_name, function f, double a, double b){
	printf("\nIntegral of the given function between %.3lf and %.3lf(%s method)\n", a, b, method_name);
	printf("---------------\n");
	printf("With %d iterations: %lf\n", MAX_ITERATIONS, _method(f, MAX_ITERATIONS, a, b));
	printf("With %d iterations: %lf\n", MAX_ITERATIONS*2, _method(f, MAX_ITERATIONS*2, a, b));
	printf("---------------\n");
}

/*
 * @brief: The main function
 * @args: Takes no argument
*/

int main(void){
	int a, b;
	printf("Integration bounds(separated by a space): ");
	scanf("%d %d", &a, &b);
	example(trapezoid, "trapezoid", f1, a, b);
	example(simpson1_3, "simpson 1/3", f1, a, b);
	example(simpson1_3, "simpson 3/8", f1, a, b);
	example(mid_point_rule, "mid-point", f1, a, b);
	example(boole, "boole", f1, a, b);
	return EXIT_SUCCESS;
}

/*
 The section below contains different methods to solve definite integrals
*/

// Closed methods section

double trapezoid(function f, int n, double a, double b){
	double area=0.5*f(a)+0.5*f(b), h = (b-a)/n;

	for (int i=1; i<n; i++){
		area += f(a + i*h);
	}
	return area*h;
}

// Simpson 1/3
double simpson1_3(function f, int n, double a, double b){
	double area=f(a)+f(b), h = (b-a)/n;

	for (int i=1; i<n; i++){
		area += i%2 == 0 ? 2*f(a + i*h) : 4*f(a + i*h);
	}
	return (area*h)/3;
}

double simpson3_8(function f, int n, double a, double b){
	double area=f(a)+f(b), h = (b-a)/n;

	for (int i=1; i<n; i++){
		area += i%3 == 0 ? 2*f(a + i*h) : 3*f(a + i*h);
	}
	return 3*(area*h)/8;
}

double boole(function f, int n, double a, double b){
	double area=7*(f(a)+f(b)), h = (b-a)/n;

	for (int i=1; i<n; i++){
		if (i%2 != 0)
			area += 32*f(a + i*h);
		else{
			area += i%4 == 0 ? 14*f(a + i*h) : 12*f(a + i*h);
		}
	}
	return 2*(area*h)/45;
}

// Open methods section

double mid_point_rule(function f, int n, double a, double b){
	double h = (b-a)/n;
	double area = 0.0;
	for (int i=0; i<n; i++){
		area += f((a+h/2.0)+i*h);
	}
	return area*h;
}
