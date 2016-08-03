
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

	/*
	int waterGrid[3] = { 16, 16, 16 };

	//--- Inpuut checks ---// (FIXME add error messages)
	if ( waterGrid[0] % 2 != 0 )
	{
		exit(1);
	}
	else if ( waterGrid[1] % 4 != 0 )
	{
		exit(1);
	}
	else if ( waterGrid[2] % 2 != 0 )
	{
		exit(1);
	}
	*/

	//---- Geometric constants -----//

	// Water geometry
	double PI = 3.14159265358979323846;
	double theta = 109.47 * 2.0*PI/360.0; // [radians]
	double l = 0.34;  // distance between nearest-neighbor oxygens [nm]
	double l_OH = 0.1; // length of O-H bonds [nm]

	// Derived geometric constants
	double l_xy = l*sqrt(2.0/3.0*(1.0 - cos(theta)));
	double dx = sqrt(3.0)/2.0*l_xy;
	double dy = l_xy/2.0;
	double dz = 0.5*sqrt(l*l - l_xy*l_xy);

	//--- Construct the HCP unit cell for the oxygens ---//
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
	int unitCellGrid[3] = { 13, 6, 7 };

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
	double dist, distSq, pairDist = 1.01*l, scale_factor;
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

		// FIXME
		/*
		if ( (hydrogenCounts[owner] > hydrogenCounts[partner]) && (!isSwapAccepted) )
		{
			std::cerr << "Swap NOT accepted when it should have been!" << std::endl;
		}
		else if ( (hydrogenCounts[owner] < hydrogenCounts[partner]) && (isSwapAccepted) )
		{
			std::cerr << "Swap accepted when it should NOT have been!" << std::endl;
		}
		*/
		/*
		std::cout << "  ATTEMPTED SWAP" << std::endl;
		std::cout << "    Owner:   " << owner << "(has " << hydrogenCounts[owner] << " H's)" << std::endl;
		std::cout << "    Partner: " << partner << "(has " << hydrogenCounts[partner] << " H's)" << std::endl;
		std::cout << "    SwapAccepted: " << isSwapAccepted << std::endl;
		*/
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
	int    numAtoms = 4*numWaters, count, atomIndex;
	rvec   x_OW, x_MW, x_HW1, x_HW2;
	rvec*  x_all = (rvec*) malloc(numAtoms*sizeof(rvec));
	double a_M = 0.1345833509;

	for ( int i=0; i<numWaters; ++i )
	{
		// Oxygen
		for ( int d=0; d<3; ++d )
		{
			x_OW[d] = x_HCP[i][d];
		}

		// Set the positions of this O's hydrogens
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
				minImage(x_HCP[owner], x_HCP[partner], boxL, x_i_j, distSq);
				dist = sqrt(distSq);
				scale_factor = l_OH/dist;

				// HW1
				if ( count == 0 )
				{
					for ( int d=0; d<3; ++d ) 
					{
						x_HW1[d] = x_HCP[owner][d] + scale_factor*x_i_j[d];
					}
					keepInBox(boxL, x_HW1);

					++count;
				}
				// HW2
				else if ( count == 1 )
				{
					for ( int d=0; d<3; ++d ) 
					{
						x_HW2[d] = x_HCP[owner][d] + scale_factor*x_i_j[d];
					}
					keepInBox(boxL, x_HW2);

					++count;
				}
				// Done
				else
				{
					break;
				}
			}
		}

		// Construct virtual site "M"
		for ( int d=0; d<3; ++d )
		{
			x_MW[d] = (1.0 - 2.0*a_M)*x_OW[d] + a_M*(x_HW1[d] + x_HW2[d]);
		}
		keepInBox(boxL, x_MW);

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

	//----- Write .gro file -----//
	// - Note: .gro file indexing starts with 1!

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

	/*
	for ( int i=0; i<numH; ++i )
	{
		free(pairs[i]);
	}
	free(pairs);
	*/

	/*
	if ( argc < 2 )
	{
		std::cerr << "Error: Must submit pdb file." << std::endl;
		exit(1);
	}

	// Input file
	std::string pdbFile(argv[1]);
	int numAtomsPdb = 1296;
	int numWaters   = 432;

	// Allocate memory for input array
	rvec* x_pdb = (rvec *) malloc(numAtomsPdb*sizeof(rvec));



	//----- Read .pdb file -----//
	std::ifstream 	   ifs(pdbFile);
	std::istringstream iss;
	std::string   	   line, atomType, atomName, atomClass;
	int			  	   atomIndex, moleculeIndex;

	int atomCounter = 0;
	while ( getline(ifs, line) && (atomCounter < numAtomsPdb) )
	{
		iss.str(line);

		iss >> atomClass; iss >> atomIndex; 
		iss >> atomName;  iss >> atomType;
		iss >> moleculeIndex;

		// {x,y,z} coordinates [Angstroms]
		iss >> x_pdb[atomCounter][0];
		iss >> x_pdb[atomCounter][1];
		iss >> x_pdb[atomCounter][2];

		atomCounter++;
	}
	ifs.close();

	// Convert from [Angstroms] --> [nm]
	for ( int i=0; i<numAtomsPdb; ++i )
	{
		for ( int j=0; j<3; ++j )
		{
			x_pdb[i][j] /= 10.0;
		}
	}

	// Find min and max x,y,z --> get box lengths
	float ranges[3][2]; // x,y,z and min, max
	for ( int i=0; i<3; ++i )
	{
		for ( int j=0; j<2; ++j )
		{
			ranges[i][j] = 0.0;
		}
	}

	for ( int i=0; i<numAtomsPdb; ++i )
	{
		for ( int d=0; d<3; ++d )
		{
			// Min
			if ( x_pdb[i][d] < ranges[d][0] )
			{
				ranges[d][0] = x_pdb[i][d];
			}

			// Max
			if ( x_pdb[i][d] > ranges[d][1] )
			{
				ranges[d][1] = x_pdb[i][d];
			}

		}
	}

	// Print results
	for ( int d=0; d<3; ++d )
	{
		std::cout << "\td=" << d << ": min=" << ranges[d][0]
				  << ", max=" << ranges[d][1] << std::endl;
	}

	// Shift the entire box so that the origin is at a corner
	rvec boxL, x_shift;
	for ( int d=0; d<3; ++d )
	{
		boxL[d]    = 1.04*(ranges[d][1] - ranges[d][0]);
		x_shift[d] = boxL[d]/2.0; // old box has origin at center
	}

	for ( int i=0; i<numAtomsPdb; ++i )
	{
		for ( int d=0; d<3; ++d )
		{
			x_pdb[i][d] += x_shift[d];
		}
	}

	//----- Write .gro file -----//
	// - Note: .gro file indexing starts with 1!

	// Include virtual sites
	int    numAtomsGro = 4*numWaters;	// Includes M sites
	int	   oldIndex;
	float  a = 0.1345833509;
	float* x_O, * x_HW1, * x_HW2;
	float  x_M[3];

	FILE* pGroFile;
	pGroFile = fopen("outconf.gro", "w");
	fprintf( pGroFile, "Box of TIP4P-type ice derived from LSBU site's pdb file \n");
	fprintf( pGroFile, "    %d\n", numAtomsGro );

	atomCounter = 0;
	for ( int i=0; i<numWaters; ++i )
	{
		// Index of molecule in new .gro file
		moleculeIndex = i+1;

		// OW
		oldIndex = 3*i;			// Index in the old SPC/E-type array
		atomIndex = 4*i + 1;	// Index in the new .gro file
		x_O = x_pdb[oldIndex];

		fprintf( pGroFile, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
							moleculeIndex, "ICE", "OW", atomIndex,
						    x_O[0], x_O[1], x_O[2] );

		// HW1
		++oldIndex;
		++atomIndex;
		x_HW1 = x_pdb[oldIndex];

		fprintf( pGroFile, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
							moleculeIndex, "ICE", "HW1", atomIndex,
						    x_HW1[0], x_HW1[1], x_HW1[2] );

		// HW2
		++oldIndex;
		++atomIndex;
		x_HW2 = x_pdb[oldIndex];

		fprintf( pGroFile, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
							moleculeIndex, "ICE", "HW2", atomIndex,
							x_HW2[0], x_HW2[1], x_HW2[2] );

		// MW
		++atomIndex;
		++atomIndex;
		for ( int d=0; d<3; ++d )
		{
			x_M[d] = (1.0 - 2.0*a)*x_O[d] + a*x_HW1[d] + a*x_HW2[d];
		}
		fprintf( pGroFile, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
							moleculeIndex, "ICE", "MW", atomIndex,
							x_M[0], x_M[1], x_M[2] );

	}
	fprintf( pGroFile, "   %2.5f   %2.5f   %2.5f\n", 
					   boxL[0], boxL[1], boxL[2] );

	fclose(pGroFile);

	// Cleanup
	free(x_pdb);
	*/
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
