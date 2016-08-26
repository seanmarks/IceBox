
#ifndef MAIN_H
#define MAIN_H

// Position vector
typedef double rvec[3];

// Contains the indices of two atoms which constitute a pair
typedef int Pair[2];

// Apply the minimum image convention
void minImage(const rvec x1, const rvec x2, const rvec boxL, rvec x12, double& distSq);

// Keep the atom in the simulaton box
void keepInBox(const rvec boxL, rvec x);

// Checks whether each oxygen has 2 hydrogens
bool areHydrogensCorrectlyPlaced(std::vector<int> hydrogenCounts);

// Vector norm
double norm(const rvec x);

// Use the dot product to get the angle between the vectors
double angleBetweenVectors(const rvec a, const rvec b);

#endif // MAIN_H
