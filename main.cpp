
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

// Project headers
#include "main.h"

int main(int argc, char** argv)
{
	//--- Input ---//

	// FIXME Input number of waters

	//---- Geometric constants -----//

	// Water geometry: TIP4P/Ice
	double PI = 3.14159265358979323846;
	double theta = 109.47 * 2.0*PI/360.0; 	  // tetrahedral bond angle [radians]
	double theta_HOH = 104.52 * 2.0*PI/360.0; // H-O-H bond angle [radians]
	double l = 0.275;      		// distance between nearest-neighbor oxygens [nm]
	double l_OH = 0.09572; 		// length of O-H bonds [nm]
	double l_OM = 0.01577;		// O-M distance [nm]

	// Derived geometric constants
	double l_xy = l*sqrt(2.0/3.0*(1.0 - cos(theta)));
	double dx = sqrt(3.0)/2.0*l_xy;
	double dy = l_xy/2.0;
	double dz = 0.5*sqrt(l*l - l_xy*l_xy);

	// - Distance between H's [nm]
	double l_HH = l_OH*sqrt(2.0*(1.0 - cos(theta_HOH))); // ~0.15139 nm

	//--- Construct the HCP unit cell for the oxygens ---//
	// - Start with perfectly tetrahedral geometry
	std::cout << "  IceBox: Constructing HCP unit cell for oxygens." << std::endl;
	rvec unitCell[8];

	// Unit cell box lengths
	rvec unitCellBoxL;
	unitCellBoxL[0] = 2.0*dx;
	unitCellBoxL[1] = 2.0*(l_xy + dy);
	unitCellBoxL[2] = 2.0*(l + 2.0*dz);

	// Lower z-plane
	// 
	//          d
	//          |
	//          u
	//        /
	//      d
	//      |
	//      u
	//
	unitCell[0][0] = 0.0;	
	unitCell[0][1] = 0.0;	
	unitCell[0][2] = 0.0;

	unitCell[1][0] = unitCell[0][0];	
	unitCell[1][1] = unitCell[0][1] + l_xy;
	unitCell[1][2] = unitCell[0][2];
	
	unitCell[2][0] = unitCell[1][0] + dx;
	unitCell[2][1] = unitCell[1][1] + dy;
	unitCell[2][2] = unitCell[1][2];

	unitCell[3][0] = unitCell[2][0];
	unitCell[3][1] = unitCell[2][1] + l_xy;
	unitCell[3][2] = unitCell[2][2];

	// Upper z-plane: directly above the first one
	// 
	//          u
	//          |
	//          d
	//        /
	//      u
	//      |
	//      d
	//
	for ( int i=0; i<4; ++i )
	{
		for ( int j=0; j<3; ++j )
		{
			unitCell[i + 4][j] = unitCell[i][j];
		}
		unitCell[i + 4][2] += l + 2.0*dz; // spacing between basal planes
	}

	// Shift all atoms up/down relative to the basal plane (see diagrams above)
	for ( int i=0; i<8; ++i )
	{
		// Shift up these atoms
		if ( (i == 0) || (i == 2) || (i == 5) || (i == 7 ) )
		{
			unitCell[i][2] += dz;
		}
		// Shift down the rest
		else
		{
			unitCell[i][2] -= dz;
		}
	}

	// Offset all atoms so that none lie on the edges of the unit cell box
	rvec x_shift;
	x_shift[0] = dx/2.0;
	x_shift[1] = dy/2.0;
	x_shift[2] = l/2.0;

	for ( int i=0; i<8; ++i )
	{
		for ( int j=0; j<3; ++j )
		{
			unitCell[i][j] += x_shift[j];
		}
	}

	//--- Replicate the unit cell to produce the oxygen HCP lattice ---//
	// - To get about 4,000 ice molecules in a roughly cubic volume, use:
	//   			int unitCellGrid[3] = { 12, 6, 7 };
	int unitCellGrid[3] = { 13, 6, 7 };
	//int unitCellGrid[3] = { 2, 2, 2 };

	std::cout << "  IceBox: Constructing full HCP lattice of oxygens." << std::endl;

	// Eight waters per unit cell
	int numWaters = 8*unitCellGrid[0]*unitCellGrid[1]*unitCellGrid[2];

	// Allocate memory
	rvec* x_HCP = (rvec*) malloc( numWaters*sizeof(rvec) );
	int atomCounter = 0;
	for ( int i=0; i<unitCellGrid[0]; ++i ) // x cells
	{
		for ( int j=0; j<unitCellGrid[1]; ++j ) // y cells
		{
			for ( int k=0; k<unitCellGrid[2]; ++k ) // z cells
			{
				for ( int m=0; m<8; ++m ) // unit cell atoms
				{
					// x, y, z dimensions
					x_HCP[atomCounter][0] = unitCell[m][0] + i*unitCellBoxL[0];
					x_HCP[atomCounter][1] = unitCell[m][1] + j*unitCellBoxL[1];
					x_HCP[atomCounter][2] = unitCell[m][2] + k*unitCellBoxL[2];

					++atomCounter;
				}
			}
		}
	}

	// Simulation box lengths
	rvec boxL;
	for ( int j=0; j<3; ++j )
	{
		boxL[j] = unitCellGrid[j]*unitCellBoxL[j];
	}

	//----- Place hydrogens using a Monte Carlo routine -----//
	// Method of Buch et al., J. Phys. Chem. B 102.44 (1998)

	// Set up RNG engine
	std::seed_seq seedSequence = { 749725171 };
	std::mt19937  rng(seedSequence);	// Mersenne Twister

	// Find all O-O nearest-neighbor pairs, and randomly assign an H to one atom in each pair
	int   numH = 2*numWaters;
	int   numPairs = numH;
	Pair* pairs = (Pair*) malloc(numPairs*sizeof(Pair)); // NOTE: index over pairs = index over H's

	std::vector<int>  hydrogenCounts(numWaters, 0); 		// how many bound H's each oxygen has
	std::vector<int>  hydrogenOwners(numH, 0);	   			// which O has the H
	std::vector<bool> doesFirstAtomHaveHydrogen(numPairs);	// Does pair[][0] have the H?

	std::uniform_int_distribution<int> coin_flip(0, 1);

	int    pairIndex = 0;
	double dist, distSq, pairDist = 1.01*l;
	rvec   x_i_j; // Direction: i --> j
	for ( int i=0; i<numWaters; ++i )
	{
		for ( int j=i+1; j<numWaters; ++j )
		{
			// Minimum image vector: i --> j
			minImage(x_HCP[i], x_HCP[j], boxL, x_i_j, distSq);
			dist = sqrt(distSq); // = || x_i_j ||_2

			if ( dist <= pairDist )
			{
				// Record the pair
				pairs[pairIndex][0] = i;
				pairs[pairIndex][1] = j;

				// Randomly assign the H to one of the paired atoms
				if ( coin_flip(rng) == 0 )
				{
					// Give to 'i'
					++(hydrogenCounts[i]);
					hydrogenOwners[pairIndex] = i;
					doesFirstAtomHaveHydrogen[pairIndex] = true;
				}
				else
				{
					// Give to 'j'
					++(hydrogenCounts[j]);
					hydrogenOwners[pairIndex] = j;
					doesFirstAtomHaveHydrogen[pairIndex] = false;
				}

				++pairIndex;
			}
		}

		// numPairs == numHydrogens
		if ( pairIndex == numPairs )
		{
			break;
		}
	}

	// Stochastically reassign hydrogens until all O's have two H's each
	// - For each pair, the owner is the O closest to the H; the other O is its partner
	std::cout << "  IceBox: Stochastically assigning hydrogen atoms." << std::endl;
	int delta, delta_trial, owner, partner;
	bool isSwapAccepted;
	std::uniform_int_distribution<int> randomPair(0, numPairs - 1);

	while ( areHydrogensCorrectlyPlaced(hydrogenCounts) == false ) //TODO
	{
		pairIndex = randomPair(rng);
		owner = hydrogenOwners[pairIndex];

		// Where are the owner and partner in the pairs array?
		if ( pairs[pairIndex][0] == owner )
		{
			partner = pairs[pairIndex][1];
		}
		else
		{
			partner = pairs[pairIndex][0];
		}

		// Current difference 
		delta = abs(hydrogenCounts[owner] - hydrogenCounts[partner]);

		// Result of proposed swap
		delta_trial = abs( (hydrogenCounts[owner] - 1) - (hydrogenCounts[partner] + 1) );

		// Determine whether to accept the proposed swap
		if ( delta_trial < delta )
		{
			isSwapAccepted = true;
		}
		else if ( delta_trial == delta )
		{
			// 50/50 chance to accept anyway
			isSwapAccepted = static_cast<bool>( coin_flip(rng) );
		}
		else
		{
			isSwapAccepted = false;
		}

		// Make any necessary changes
		if ( isSwapAccepted )
		{
			// NOTE: owner/partner here refer to the roles in the "old" arrangement

			hydrogenOwners[pairIndex] = partner;

			--(hydrogenCounts[owner]);
			++(hydrogenCounts[partner]);

			if ( doesFirstAtomHaveHydrogen[pairIndex] == true )
			{
				doesFirstAtomHaveHydrogen[pairIndex] = false;
			}
			else
			{
				doesFirstAtomHaveHydrogen[pairIndex] = true;
			}
		}
	}

	//----- Construct the full coordinates array -----//
	std::cout << "  IceBox: Constructing full array of atomic positions." << std::endl;
	int    numAtoms = 4*numWaters, count, atomIndex;
	double norm_x;
	rvec   x_OW, x_MW, x_HW1, x_HW2;
	rvec   x_OO_1, x_OO_2, unitBisector, unitNormal; // working vectors
	rvec*  x_all = (rvec*) malloc(numAtoms*sizeof(rvec));

	// Length of the projection of the O-H bond onto the H-O-H unit bisector
	double l_proj = sqrt(l_OH*l_OH - l_HH*l_HH/4.0);
	
	for ( int i=0; i<numWaters; ++i )
	{
		// Oxygen
		for ( int d=0; d<3; ++d )
		{
			x_OW[d] = x_HCP[i][d];
		}

		// For the rest, first get vectors connecting the associated oxygens
		count = 0;
		for ( int j=0; j<numH; ++j )
		{
			
			if ( hydrogenOwners[j] == i ) // If the owner of hydrogen j is oxygen i...
			{
				owner = i;

				// Infer the index of the partner
				if ( doesFirstAtomHaveHydrogen[j] == true )
				{
					partner = pairs[j][1];
				}
				else
				{
					partner = pairs[j][0];
				}

				// Connecting vector: owner --> partner
				if ( count == 0 )
				{
					minImage(x_HCP[owner], x_HCP[partner], boxL, x_OO_1, distSq);
					++count;
				}
				else if ( count == 1 )
				{
					minImage(x_HCP[owner], x_HCP[partner], boxL, x_OO_2, distSq);
					++count;
				}
				else
				{
					// Found both H's
					break;
				}
			}
		}

		// Unit vector along the H-O-H bisector
		for ( int d=0; d<3; ++d ) { unitBisector[d] = x_OO_1[d] + x_OO_2[d]; }
		norm_x = norm(unitBisector);
		for ( int d=0; d<3; ++d ) { unitBisector[d] /= norm_x; }

		// Unit vector normal to the H-O-H bisector, in the H-O-H plane
		for ( int d=0; d<3; ++d ) { unitNormal[d] = x_OO_2[d] - x_OO_1[d]; }
		norm_x = norm(unitNormal);
		for ( int d=0; d<3; ++d ) { unitNormal[d] /= norm_x; }

		// Use these unit vectors to construct HW1, HW2, and M
		for ( int d=0; d<3; ++d )
		{
			x_HW1[d] = x_OW[d] + l_proj*unitBisector[d] + (l_HH/2.0)*unitNormal[d];
			x_HW2[d] = x_OW[d] + l_proj*unitBisector[d] - (l_HH/2.0)*unitNormal[d];
			x_MW[d]  = x_OW[d] + l_OM*unitBisector[d];
		}

		// DON'T correct for PBCs in this way! GROMACS likes to be given whole molecules
		// - mdrun "sees" how to handle PBCs based on the settings in the *.mdp file
		//
		//keepInBox(boxL, x_HW1);
		//keepInBox(boxL, x_HW2);
		//keepInBox(boxL, x_MW);
		//

		// Put results into composite array
		atomIndex = 4*i;
		for ( int d=0; d<3; ++d ) { x_all[atomIndex][d] = x_OW[d]; }

		++atomIndex;
		for ( int d=0; d<3; ++d ) { x_all[atomIndex][d] = x_HW1[d]; }

		++atomIndex;
		for ( int d=0; d<3; ++d ) { x_all[atomIndex][d] = x_HW2[d]; }

		++atomIndex;
		for ( int d=0; d<3; ++d ) { x_all[atomIndex][d] = x_MW[d]; }
	}

	// Check results (use 1% margin of error)
	std::cout << "  IceBox: Checking configuration." << std::endl;
	for ( int i=0; i<numWaters; ++i )
	{
		atomIndex = 4*i;

		// O-H1
		minImage(x_all[atomIndex], x_all[atomIndex+1], boxL, x_i_j, distSq);
		dist = sqrt(distSq);
		if ( dist > 1.01*l_OH )
		{
			std::cerr << "  Molecule " << i << ": OW-HW1 bond is " << dist << " nm " 
					  << "(should be " << l_OH << " nm)." << std::endl;
		}

		// O-H2
		minImage(x_all[atomIndex], x_all[atomIndex+2], boxL, x_i_j, distSq);
		dist = sqrt(distSq);
		if ( dist > 1.01*l_OH )
		{
			std::cerr << "  Molecule " << i << ": OW-HW2 bond is " << dist << " nm "
				      << "(should be " << l_OH << " nm)." << std::endl;
		}

		// H1-H2
		minImage(x_all[atomIndex+1], x_all[atomIndex+2], boxL, x_i_j, distSq);
		dist = sqrt(distSq);
		if ( dist > 1.01*l_HH )
		{
			std::cerr << "  Molecule " << i << ": HW1-HW2 bond is " << dist << " nm "
					  << "(should be " << l_HH << " nm)." << std::endl;
		}

		// O-M
		minImage(x_all[atomIndex], x_all[atomIndex+3], boxL, x_i_j, distSq);
		dist = sqrt(distSq);
		if ( dist > 1.01*l_OM ) 
		{
			std::cerr << "  Molecule " << i << ": OW-MW bond is " << dist << " nm "
					  << "(should be " << l_OM << " nm)." << std::endl;
		}
	}

	//----- Write .gro file -----//
	// - Note: .gro file indexing starts with 1!
	std::cout << "  IceBox: Writing .gro file." << std::endl;

	// File name
	std::ostringstream ossGroFileName;
	ossGroFileName << "initconf_ice_" << numWaters << ".gro";
	std::string groFileName( ossGroFileName.str() );

	FILE* pGroFile;
	pGroFile = fopen(groFileName.c_str(), "w");
	fprintf( pGroFile, "TIP4P/Ice in a box\n" );
	fprintf( pGroFile, "    %d\n", numAtoms );

	for ( int i=0; i<numWaters; ++i )
	{
		// OW
		atomIndex = 4*i;	   // Index in the coordinates array "x_all"

		fprintf( pGroFile, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
							i+1, "ICE", "OW", atomIndex+1,
						    x_all[atomIndex][0], x_all[atomIndex][1], x_all[atomIndex][2] );

		// HW1
		++atomIndex;

		fprintf( pGroFile, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
							i+1, "ICE", "HW1", atomIndex+1,
						    x_all[atomIndex][0], x_all[atomIndex][1], x_all[atomIndex][2] );

		// HW2
		++atomIndex;

		fprintf( pGroFile, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
							i+1, "ICE", "HW2", atomIndex+1,
							x_all[atomIndex][0], x_all[atomIndex][1], x_all[atomIndex][2] );

		// MW
		++atomIndex;

		fprintf( pGroFile, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
							i+1, "ICE", "MW", atomIndex+1,
							x_all[atomIndex][0], x_all[atomIndex][1], x_all[atomIndex][2] );

	}
	fprintf( pGroFile, "   %2.5f   %2.5f   %2.5f\n", 
					   boxL[0], boxL[1], boxL[2] );

	fclose(pGroFile);

	//----- Cleanup -----//
	free(x_HCP);
	free(x_all);
	free(pairs);

}

// Applies the minimum image convention
void minImage(const rvec x1, const rvec x2, const rvec boxL, rvec x12, double& distSq)
{
	distSq = 0.0;

	for ( int d=0; d<3; d++ )
	{
		x12[d] = x2[d] - x1[d];
		// Apply minimum image convention
		if      ( x12[d] >  0.5*boxL[d] ) { x12[d] -= boxL[d]; }
		else if ( x12[d] < -0.5*boxL[d] ) { x12[d] += boxL[d]; }

		distSq += x12[d]*x12[d];
	}

	return;
}

// Keeps the atom in the simulaton box by applying PBCs
void keepInBox(const rvec boxL, rvec x)
{
	for ( int d=0; d<3; d++ )
	{
		// Apply minimum image convention
		if      ( x[d] > boxL[d] ) { x[d] -= boxL[d]; }
		else if ( x[d] < 0.0 )     { x[d] += boxL[d]; }
	}

	return;

}

// Checks whether each oxygen has 2 hydrogens
bool areHydrogensCorrectlyPlaced(std::vector<int> hydrogenCounts)
{
	int numWaters = hydrogenCounts.size();
	for ( int i=0; i<numWaters; i++ )
	{
		if ( hydrogenCounts[i] != 2 )
		{
			return false;
		}
	}
	return true;
}

// Vector norm
double norm(const rvec x)
{
	return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}
