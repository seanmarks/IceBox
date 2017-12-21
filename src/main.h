
#ifndef MAIN_H
#define MAIN_H

#include <array>
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
typedef int Pair[2];

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

#endif // MAIN_H
