
#ifndef MAIN_H
#define MAIN_H

#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <exception>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#define DIM 3

// Position vector
using Real3 = std::array<double,DIM>;

// Contains the indices of two atoms which constitute a pair
using Pair = std::array<int,2>;

// Apply the minimum image convention
void minImage(const Real3& x1, const Real3& x2, const Real3& boxL, Real3& x12, double& distSq);

// Keep the atom in the simulaton box
void keepInBox(const Real3& boxL, Real3& x);

// Checks whether each oxygen has 2 hydrogens
bool areHydrogensCorrectlyPlaced(const std::vector<int>& hydrogenCounts);

// Vector norm
double norm(const Real3& x);

// Use the dot product to get the angle between the vectors
double angleBetweenVectors(const Real3& a, const Real3& b);

// Read a single frame from a .gro file 
void readFrameFromGroFile(
	const std::string& gro_file,
	// Output
	std::vector<Real3>& coords, 
	std::vector<std::string>& atomTypes,
	std::vector<int>& atomSerials, 
	Real3& boxL
);

// Check that the given structure is acceptable
void check_structure(
	const double tol,  // tolerance (fractional error permitted, percent/100)
	// Frame
	const std::vector<Real3>& x_all,
	const Real3& boxL,
	// Model parameters
	const std::vector<double>& charges,  // charges for each atom in a water molecule
	const double l_OH, 
	const double l_HH, 
	const double l_OM,
	// Structure 
	const double l_HB,   // hydrogen bond length
	const double r_nn,   // distance between oxygen nearest neighbors
	const double theta,  // tetrahedral angle
	// Output
	Real3& dipole,  // net dipole
	double& mu      // norm of net dipole
);


#endif // MAIN_H
