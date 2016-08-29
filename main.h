
#ifndef MAIN_H
#define MAIN_H

// Position vector
typedef double Rvec[3];

// Contains the indices of two atoms which constitute a pair
typedef int Pair[2];

// Apply the minimum image convention
void minImage(const Rvec x1, const Rvec x2, const Rvec boxL, Rvec x12, double& distSq);

// Keep the atom in the simulaton box
void keepInBox(const Rvec boxL, Rvec x);

// Checks whether each oxygen has 2 hydrogens
bool areHydrogensCorrectlyPlaced(std::vector<int> hydrogenCounts);

// Vector norm
double norm(const Rvec x);

// Use the dot product to get the angle between the vectors
double angleBetweenVectors(const Rvec a, const Rvec b);

#endif // MAIN_H
