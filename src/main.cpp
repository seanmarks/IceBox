// Icebox.exe
//
// INPUT:
// - Phase
//   - I_hex/hexagonal/I_h
//   - cubic/I_c
// - Number of unit cells along each axis (Nx Ny Nz)
//   - To get about 4,000 ice molecules in a roughly cubic volume (L ~ 5 nm) ...
//     - hex:    13  6  7
//     - cubic:   8  8  8
//
// TODO:
// - Fix issues placing H's with unitCellGrid = {1, 1, 1} (possible with this algorithm?)
// - Fix fixed-format output when #atoms or #waters is >5 char long in decimal format
// - Change atom names to "ice" (e.g. OI instead of OW)

#include "main.h"

int main(int argc, char* argv[])
{
	//----- Input -----//

	// FIXME Input number of waters or something like that so the user can specify a kind of scale
	if ( argc < 2 ) {
		std::cout << "  Error: Must submit an ice phase type (\"I_h\" or \"I_c\")" << "\n";
		exit(1);
	}
	else if ( argc < 5 ) {
		std::cout << "  Error: Must provide number of unit cells along each axis (Nx Ny Nz)\n";
		exit(1);
	}

	// Crystal polymorph
	std::string phase(argv[1]);

	// Unit cell mesh
	int unitCellGrid[DIM];
	unitCellGrid[0] = std::atoi( argv[2] );
	unitCellGrid[1] = std::atoi( argv[3] );
	unitCellGrid[2] = std::atoi( argv[4] );

	//----- Constants -----//

	// Water geometry: TIP4P/Ice
	const double PI = 3.14159265358979323846;
	double theta = 109.47 * 2.0*PI/360.0; 	  // tetrahedral bond angle [radians]
	double theta_HOH = 104.52 * 2.0*PI/360.0; // H-O-H bond angle [radians]
	double l = 0.275;         // distance between nearest-neighbor oxygens [nm]
	                          // - About the same for I_h and I_c
	double l_OH = 0.09572; // length of O-H bonds [nm]
	double l_OM = 0.01577; // O-M distance [nm]

	// Distance between H's [nm]
	double l_HH = l_OH*sqrt(2.0*(1.0 - cos(theta_HOH))); // ~0.15139 nm

	// Partial charges
	double q_H = 0.520;
	double q_M = -2.0*q_H;

	//----- Construct the unit cell for the oxygens -----//

	Real3* unitCell;      // Unit cell oxygen positions
	Real3  unitCellBoxL;	// Unit cell box size lengths
	int numOxygensPerUnitCell;

	if ( phase == "I_h" or phase == "hex" or phase == "hexagonal" ) { // Hexagonal ice
		// Credit for approach:
		//	Johannes Zierenberg
		//  Master's Thesis: "Tip4p Water Model in the Ice Ih Configuration"
		//  Supervisor: Dr. Wolfhard Janke
		//	Institut für Theoretische Physik
		//	Fakultät für Physik und Geowissenschaften, Universität Leipzig

		// Derived constants for this method
		double l_xy = l*sqrt(2.0/3.0*(1.0 - cos(theta)));
		double dx = sqrt(3.0)/2.0*l_xy;
		double dy = l_xy/2.0;
		double dz = 0.5*sqrt(l*l - l_xy*l_xy);

		// Start with perfectly tetrahedral geometry
		std::cout << "  IceBox: Constructing HCP unit cell for oxygens." << std::endl;
		numOxygensPerUnitCell = 8;
		unitCell = new Real3[numOxygensPerUnitCell];

		// Unit cell box lengths
		unitCellBoxL[0] = 2.0*dx;
		unitCellBoxL[1] = 2.0*(l_xy + dy);
		unitCellBoxL[2] = 2.0*(l + 2.0*dz);

		// Lower z-plane
		// - u/d => atom is shifted up/down relative the "midpoint plane"
		//   
		// 
		//              d               3
		//              |               |
		//              u               2
		//            /               /
		//          d               1
		//    y     |               |
		//    ^     u               0
		//    |
		//    ---> x
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
		for ( int i=0; i<4; ++i ) {
			for ( int j=0; j<DIM; ++j ) {
				unitCell[i + 4][j] = unitCell[i][j];
			}
			unitCell[i + 4][2] += l + 2.0*dz; // spacing between planes
		}

		// Shift all atoms up/down relative to the plane (see diagrams above)
		for ( int i=0; i<numOxygensPerUnitCell; ++i ) {
			// Shift up these atoms
			if ( (i == 0) || (i == 2) || (i == 5) || (i == 7 ) ) {
				unitCell[i][2] += dz;
			}
			// Shift down the rest
			else {
				unitCell[i][2] -= dz;
			}
		}

		// Offset all atoms so that none lie on the edges of the unit cell box
		Real3 x_shift = {{ dx/2.0, dy/2.0, l/2.0 }};
		for ( int i=0; i<numOxygensPerUnitCell; ++i ) {
			for ( int j=0; j<DIM; ++j ) {
				unitCell[i][j] += x_shift[j];
			}
		}
	}
	else if ( phase == "I_c" or phase == "cubic" ) { // Cubic ice
		// Derived constants
		double a = 4.0/sqrt(3.0)*l; // side length of FCC unit cell [nm]
		double b = sqrt(8.0/3.0)*l; // side length of "internal" tetrahedra [nm]
		double d = b/sqrt(3.0);     // distance from center of cell to oxygens off of the FCC lattice

		// Start with perfectly tetrahedral geometry
		std::cout << "  IceBox: Constructing diamond unit cell for oxygens." << std::endl;
		numOxygensPerUnitCell = 8;
		unitCell = new Real3[numOxygensPerUnitCell];

		// Unit cell box lengths
		unitCellBoxL = {{ a, a, a }};

		// First, place the atoms which fall on a regular FCC lattice (easy)
		unitCell[0] = {{ 0.0,   0.0,   0.0   }};
		unitCell[1] = {{ a/2.0, 0.0,   a/2.0 }};
		unitCell[2] = {{ 0.0,   a/2.0, a/2.0 }};
		unitCell[3] = {{ a/2.0, a/2.0, a     }};

		//----- Now place the 4 non-FCC oxygens -----//

		Real3 x_center = {{ a/2.0, a/2.0, a/2.0 }}; // center of unit cell
		Real3 x_corner; // opposing corner
		int   index;    // index in unit cell

		// Top-front
		x_corner = {{ 0.0, 0.0, a }};
		index = 4;
		for ( int j=0; j<DIM; ++j ) {
			unitCell[index][j] = x_center[j] + 0.5*(x_corner[j] - x_center[d]);
		}

		// Top-back
		x_corner = {{ a, a, a }};
		index = 5;
		for ( int j=0; j<DIM; ++j ) {
			unitCell[index][j] = x_center[j] + 0.5*(x_corner[j] - x_center[d]);
		}

		// Bottom-right
		x_corner = {{ a, 0.0, 0.0 }};
		index = 6;
		for ( int j=0; j<DIM; ++j ) {
			unitCell[index][j] = x_center[j] + 0.5*(x_corner[j] - x_center[d]);
		}

		// Bottom-left
		x_corner = {{ 0.0, a, 0.0 }};
		index = 7;
		for ( int j=0; j<DIM; ++j ) {
			unitCell[index][j] = x_center[j] + 0.5*(x_corner[j] - x_center[d]);
		}
	}
	else {
		// Error TODO
		std::cerr << "Error - polymorph \"" << phase << "\" not recognized.\n";
		exit(1);
	}


	//----- Replicate the unit cell to produce the oxygen lattice -----//

	std::cout << "  IceBox: Constructing full lattice of oxygens." << std::endl;

	// Eight waters per unit cell (for both I_h and I_c)
	int numUnitCells = unitCellGrid[0]*unitCellGrid[1]*unitCellGrid[2];
	int numWaters = numOxygensPerUnitCell*numUnitCells;

	std::cout << "  IceBox: The lattice will have " << numWaters << " waters." << std::endl;

	// Allocate memory
	Real3* x_O_lattice;
	try {
		x_O_lattice = new Real3[numWaters];
	}
	catch (std::bad_alloc& ba) {
		std::cerr << "  IceBox: Unable to allocate memory for oxygen lattice "
		          << "(exception: " << ba.what() << ")." << "\n";
		exit(1);
	}

	// Replicate unit cell
	int atomCounter = 0;
	for ( int i=0; i<unitCellGrid[0]; ++i ) { // x cells
		for ( int j=0; j<unitCellGrid[1]; ++j ) { // y cells
			for ( int k=0; k<unitCellGrid[2]; ++k ) { // z cells
				for ( int m=0; m<numOxygensPerUnitCell; ++m ) { // unit cell atoms
					// x, y, z dimensions
					x_O_lattice[atomCounter][0] = unitCell[m][0] + i*unitCellBoxL[0];
					x_O_lattice[atomCounter][1] = unitCell[m][1] + j*unitCellBoxL[1];
					x_O_lattice[atomCounter][2] = unitCell[m][2] + k*unitCellBoxL[2];

					++atomCounter;
				}
			}
		}
	}

	// Simulation box lengths
	Real3 boxL;
	for ( int j=0; j<DIM; ++j ) {
		boxL[j] = unitCellGrid[j]*unitCellBoxL[j];
	}

	//----- Place hydrogens using a Monte Carlo routine -----//
	// Method: Buch et al., J. Phys. Chem. B 102.44 (1998)

	std::cout << "  Icebox: Finding all O-O nearest-neighbor pairs" << "\n";

	// Set up RNG engine
	//std::seed_seq seed_sequence = { 749725171 };
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::seed_seq seed_sequence = { seed };
	std::mt19937 rng(seed_sequence);	// Mersenne Twister

	// Find all O-O nearest-neighbor pairs, and randomly assign an H to one atom in each pair
	int   numH = 2*numWaters;
	int   numPairs = numH;
	Pair* pairs = new Pair[numPairs];  // NOTE: index over pairs = index over H's

	std::vector<int>  hydrogenCounts(numWaters, 0); 		// how many bound H's each oxygen has
	std::vector<int>  hydrogenOwners(numH, 0);	   			// which O has the H
	std::vector<bool> doesFirstAtomHaveHydrogen(numPairs);	// Does pair[][0] have the H?

	std::uniform_int_distribution<int> coin_flip(0, 1);

	// Working variables
	int    pairIndex = 0;
	double dist, distSq;
	double pairDist = 1.001*l; // Distance between nearest-neighbor oxygens
	Real3   x_i_j;              // Direction: i --> j

	for ( int i=0; i<numWaters; ++i ) {
		for ( int j=i+1; j<numWaters; ++j ) {
			// Minimum image vector: i --> j
			minImage(x_O_lattice[i], x_O_lattice[j], boxL, x_i_j, distSq);
			dist = sqrt(distSq); // = || x_i_j ||_2

			if ( dist <= pairDist ) {
				// Record the pair
				pairs[pairIndex][0] = i;
				pairs[pairIndex][1] = j;

				// Randomly assign the H to one of the paired atoms
				if ( coin_flip(rng) == 0 ) {
					// Give to 'i'
					++(hydrogenCounts[i]);
					hydrogenOwners[pairIndex] = i;
					doesFirstAtomHaveHydrogen[pairIndex] = true;
				}
				else {
					// Give to 'j'
					++(hydrogenCounts[j]);
					hydrogenOwners[pairIndex] = j;
					doesFirstAtomHaveHydrogen[pairIndex] = false;
				}

				++pairIndex;
			}
		}

		// numPairs == numHydrogens
		if ( pairIndex == numPairs ) {
			break;
		}
	}

	// Stochastically reassign hydrogens until all O's have two H's each
	// - For each pair, the owner is the O closest to the H; the other O is its partner
	std::cout << "  IceBox: Stochastically assigning hydrogen atoms." << std::endl;
	int delta, delta_trial, owner, partner;
	bool isSwapAccepted;
	std::uniform_int_distribution<int> randomPair(0, numPairs - 1);
	while ( areHydrogensCorrectlyPlaced(hydrogenCounts) == false ) {
		pairIndex = randomPair(rng);
		owner = hydrogenOwners[pairIndex];

		// Where are the owner and partner in the pairs array?
		if ( pairs[pairIndex][0] == owner ) {
			partner = pairs[pairIndex][1];
		}
		else {
			partner = pairs[pairIndex][0];
		}

		// Current difference 
		delta = abs(hydrogenCounts[owner] - hydrogenCounts[partner]);

		// Result of proposed swap
		delta_trial = abs( (hydrogenCounts[owner] - 1) - (hydrogenCounts[partner] + 1) );

		// Determine whether to accept the proposed swap
		if ( delta_trial < delta ) {
			isSwapAccepted = true;
		}
		else if ( delta_trial == delta ) {
			// 50/50 chance to accept anyway
			isSwapAccepted = static_cast<bool>( coin_flip(rng) );
		}
		else {
			isSwapAccepted = false;
		}

		// Make any necessary changes
		if ( isSwapAccepted ) {
			// NOTE: owner/partner here refer to the roles in the "old" arrangement

			hydrogenOwners[pairIndex] = partner;

			--(hydrogenCounts[owner]);
			++(hydrogenCounts[partner]);

			if ( doesFirstAtomHaveHydrogen[pairIndex] == true ) {
				doesFirstAtomHaveHydrogen[pairIndex] = false;
			}
			else {
				doesFirstAtomHaveHydrogen[pairIndex] = true;
			}
		}
	}

	//----- Construct the full coordinates array -----//

	std::cout << "  IceBox: Constructing full array of atomic positions." << std::endl;
	int    numAtoms = 4*numWaters, count, atomIndex;
	double norm_x;
	Real3   zero3 = {{ 0.0, 0.0, 0.0 }};
	Real3   x_OW, x_MW, x_HW1, x_HW2;
	Real3   x_OO_1 = zero3, x_OO_2 = zero3, unitBisector, unitNormal; // working vectors
	Real3*  x_all = new Real3[numAtoms];

	// Length of the projection of the O-H bond onto the H-O-H unit bisector
	double l_proj = sqrt(l_OH*l_OH - l_HH*l_HH/4.0);
	
	for ( int i=0; i<numWaters; ++i ) {
		// Oxygen
		x_OW = x_O_lattice[i];

		// For the rest, first get vectors connecting the associated oxygens
		count = 0;
		for ( int j=0; j<numH; ++j ) {
			
			if ( hydrogenOwners[j] == i ) { // If the owner of hydrogen j is oxygen i...
				owner = i;

				// Infer the index of the partner
				if ( doesFirstAtomHaveHydrogen[j] == true ) {
					partner = pairs[j][1];
				}
				else {
					partner = pairs[j][0];
				}

				// Connecting vector: owner --> partner
				if ( count == 0 ) {
					minImage(x_O_lattice[owner], x_O_lattice[partner], boxL, x_OO_1, distSq);
					++count;
				}
				else if ( count == 1 ) {
					minImage(x_O_lattice[owner], x_O_lattice[partner], boxL, x_OO_2, distSq);
					++count;
				}
				else {
					// Found both H's
					break;
				}
			}
		}

		// Unit vector along the H-O-H bisector
		for ( int d=0; d<DIM; ++d ) { unitBisector[d] = x_OO_1[d] + x_OO_2[d]; }
		norm_x = norm(unitBisector);
		for ( int d=0; d<DIM; ++d ) { unitBisector[d] /= norm_x; }

		// Unit vector normal to the H-O-H bisector, in the H-O-H plane
		for ( int d=0; d<DIM; ++d ) { unitNormal[d] = x_OO_2[d] - x_OO_1[d]; }
		norm_x = norm(unitNormal);
		for ( int d=0; d<DIM; ++d ) { unitNormal[d] /= norm_x; }

		// Use these unit vectors to construct HW1, HW2, and M
		for ( int d=0; d<DIM; ++d ) {
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
		for ( int d=0; d<DIM; ++d ) { x_all[atomIndex][d] = x_OW[d]; }

		++atomIndex;
		for ( int d=0; d<DIM; ++d ) { x_all[atomIndex][d] = x_HW1[d]; }

		++atomIndex;
		for ( int d=0; d<DIM; ++d ) { x_all[atomIndex][d] = x_HW2[d]; }

		++atomIndex;
		for ( int d=0; d<DIM; ++d ) { x_all[atomIndex][d] = x_MW[d]; }
	}

	//----- Check final results -----//

	std::cout << "  IceBox: Checking final configuration." << std::endl;

	// Check the internal structure of each molecule (use 1% margin of error)
	for ( int i=0; i<numWaters; ++i ) {
		// Index of the oxygen
		atomIndex = 4*i;

		// O-H1 distance
		minImage(x_all[atomIndex], x_all[atomIndex+1], boxL, x_i_j, distSq);
		dist = sqrt(distSq);
		if ( dist > 1.01*l_OH ) {
			std::cerr << "  Molecule " << i << ": OW-HW1 bond is " << dist << " nm " 
					  << "(should be " << l_OH << " nm)." << std::endl;
		}

		// O-H2 distance
		minImage(x_all[atomIndex], x_all[atomIndex+2], boxL, x_i_j, distSq);
		dist = sqrt(distSq);
		if ( dist > 1.01*l_OH ) {
			std::cerr << "  Molecule " << i << ": OW-HW2 bond is " << dist << " nm "
				      << "(should be " << l_OH << " nm)." << std::endl;
		}

		// H1-H2 distance
		minImage(x_all[atomIndex+1], x_all[atomIndex+2], boxL, x_i_j, distSq);
		dist = sqrt(distSq);
		if ( dist > 1.01*l_HH ) {
			std::cerr << "  Molecule " << i << ": HW1-HW2 bond is " << dist << " nm "
					  << "(should be " << l_HH << " nm)." << std::endl;
		}

		// O-M distance
		minImage(x_all[atomIndex], x_all[atomIndex+3], boxL, x_i_j, distSq);
		dist = sqrt(distSq);
		if ( dist > 1.01*l_OM ) {
			std::cerr << "  Molecule " << i << ": OW-MW bond is " << dist << " nm "
					  << "(should be " << l_OM << " nm)." << std::endl;
		}
	}

	// Neighbor list for oxygen atoms
	std::vector<std::vector<int>> neighborList(numWaters);

	for ( int i=0; i<numWaters; ++i ) {
		for ( int j=i+1; j<numWaters; ++j ) {
			minImage(x_all[4*i], x_all[4*j], boxL, x_i_j, distSq);
			dist = sqrt(distSq);
			if ( dist <= pairDist ) {
				neighborList[i].push_back(j);
				neighborList[j].push_back(i);
			}
		}
	}

	// Check that all oxygens are tetrahedrally coordinated
	std::cout << "  IceBox: Checking oxygen lattice for tetrahedral geometry "
	          << "and number of nearest-neighbor oxygens." << "\n";

	Real3   x_i_k;
	int    numNeighbors;
	double theta_high = 1.01*theta;
	double theta_low  = 0.99*theta;
	double theta_ijk;
	for ( int i=0; i<numWaters; ++i )
	{
		numNeighbors = neighborList[i].size();
		if ( numNeighbors != 4 )
		{
			std::cerr << "  IceBox: Atom " << i << " in the oxygen lattice has "
			          << numNeighbors << " neighbors (should have 4)!" << "\n";
			exit(1);
		}

		// Choose 'j' as the first neighbor of 'i' in the list
		int j = neighborList[i][0];
		minImage(x_all[4*i], x_all[4*j], boxL, x_i_j, distSq);

		for ( int n=1; n<numNeighbors; ++n ) {
			int k = neighborList[i][n];
			minImage(x_all[4*i], x_all[4*k], boxL, x_i_k, distSq);

			theta_ijk = angleBetweenVectors(x_i_j, x_i_k);
			if ( (theta_ijk < theta_low) || (theta_high < theta_ijk) ) {
				std::cerr << "  IceBox: Angle " << 4*i << "-" << 4*j << "-" << 4*k
                          << " (between neighboring oxygen atoms) is " << theta_ijk
				          << " (should be " << theta << ")!" << "\n";
				exit(1);
			}
		}
	}

	// Check the hydrogens that are coordinated with (but not bonded to) each oxygen 
	// - Hydrogen Bond length (from H to coordinated O) [nm]
	double l_HB = sqrt( (l - l_proj)*(l - l_proj) + l_HH*l_HH/4.0 );

	pairDist = 1.001*l_HB;
	int numCoordinatedHydrogens;
	for ( int i=0; i<numWaters; ++i ) {
		numCoordinatedHydrogens = 0;

		for ( int j=0; j<numWaters; ++j ) {
			if ( i != j ) {
				// O of water 'i' and H_k of water 'j'
				for ( int k=1; k<=2; ++k ) {
					minImage(x_all[4*i], x_all[4*j + k], boxL, x_i_k, distSq);
					dist = sqrt(distSq);

					if ( (dist <= pairDist) ) { //&& (dist >= 0.999*l_OH) )
						++numCoordinatedHydrogens;
					}
				}
			}
		}

		// Check
		if ( numCoordinatedHydrogens != 2 ) {
			std::cerr << "  IceBox: Water " << i << " is coordinated with "
			          << numCoordinatedHydrogens << " hydrogens that are not "
			          << "bonded to it (should be 2!)." << "\n";
			exit(1);
		}
	}

	//----- Calculate net dipole -----//

	std::vector<double> dipole(DIM, 0.0);
	std::vector<double> charges = { 0.0, q_H, q_H, q_M };

	for ( int i=0; i<numWaters; ++i ) {
		for ( int j=0; j<4; ++j ) {
			for ( int k=0; k<DIM; ++k ) {
				dipole[k] += charges[j]*x_all[4*i + j][k];
			}
		}
	}

	// Convert units from e*nm to Debyes
	const double E_NM_PER_DEBYE = 0.0208194; // [(e*nm)/Debye]
	for ( int k=0; k<DIM; ++k ) { dipole[k] /= E_NM_PER_DEBYE; }

	// Norm
	double mu = 0.0;
	for ( int k=0; k<DIM; ++k ) { mu += dipole[k]*dipole[k]; }
	mu = sqrt(mu);

	// Print
	std::cout << "  IceBox: Dipole is ( "
              << dipole[0] << ", " << dipole[1] << ", " << dipole[2] << " ) D "
 	          << " (norm: " << mu << " D total, " 
	          << mu/numUnitCells << " D per unit cell)" << "\n";


	//----- Write .gro file -----//
	// - Note: .gro file indexing starts with 1!

	std::cout << "  IceBox: Writing .gro file." << std::endl;

	// File name
	std::ostringstream ossGroFileName;
	ossGroFileName << "initconf_ice_" << phase << "_" << numWaters << ".gro";
	std::string groFileName( ossGroFileName.str() );

	FILE* pGroFile;
	pGroFile = fopen(groFileName.c_str(), "w");
	fprintf( pGroFile, "TIP4P/Ice in a box\n" );
	fprintf( pGroFile, "    %d\n", numAtoms );

	for ( int i=0; i<numWaters; ++i ) {
		// OW
		atomIndex = 4*i;	   // Index in the coordinates array "x_all"

		fprintf( pGroFile, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
							i+1, "ICE", "OI", atomIndex+1,
						    x_all[atomIndex][0], x_all[atomIndex][1], x_all[atomIndex][2] );

		// HW1
		++atomIndex;

		fprintf( pGroFile, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
							i+1, "ICE", "HI1", atomIndex+1,
						    x_all[atomIndex][0], x_all[atomIndex][1], x_all[atomIndex][2] );

		// HW2
		++atomIndex;

		fprintf( pGroFile, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
							i+1, "ICE", "HI2", atomIndex+1,
							x_all[atomIndex][0], x_all[atomIndex][1], x_all[atomIndex][2] );

		// MW
		++atomIndex;

		fprintf( pGroFile, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
							i+1, "ICE", "MI", atomIndex+1,
							x_all[atomIndex][0], x_all[atomIndex][1], x_all[atomIndex][2] );

	}
	fprintf( pGroFile, "   %2.5f   %2.5f   %2.5f\n", 
					   boxL[0], boxL[1], boxL[2] );

	fclose(pGroFile);

	std::cout << "  IceBox: Done.\n";

	//----- Cleanup -----//
	delete[] x_O_lattice;
	delete[] x_all;
	delete[] pairs;
	delete[] unitCell;
}



// Applies the minimum image convention
void minImage(const Real3& x1, const Real3& x2, const Real3& boxL, Real3& x12, double& distSq)
{
	distSq = 0.0;

	for ( int d=0; d<DIM; d++ ) {
		x12[d] = x2[d] - x1[d];
		// Apply minimum image convention
		if      ( x12[d] >  0.5*boxL[d] ) { x12[d] -= boxL[d]; }
		else if ( x12[d] < -0.5*boxL[d] ) { x12[d] += boxL[d]; }

		distSq += x12[d]*x12[d];
	}

	return;
}



// Keeps the atom in the simulaton box by applying PBCs
void keepInBox(const Real3& boxL, Real3& x)
{
	for ( int d=0; d<DIM; d++ ) {
		// Apply minimum image convention
		if      ( x[d] > boxL[d] ) { x[d] -= boxL[d]; }
		else if ( x[d] < 0.0 )     { x[d] += boxL[d]; }
	}

	return;

}



// Checks whether each oxygen has 2 hydrogens
bool areHydrogensCorrectlyPlaced(const std::vector<int>& hydrogenCounts)
{
	int numWaters = hydrogenCounts.size();
	for ( int i=0; i<numWaters; i++ ) {
		if ( hydrogenCounts[i] != 2 ) {
			return false;
		}
	}
	return true;
}



// Vector norm
double norm(const Real3& x)
{
	return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}



// Use the dot product to get the angle between the vectors
double angleBetweenVectors(const Real3& a, const Real3& b)
{
	double norm_a = 0.0, norm_b = 0.0, a_dot_b = 0.0;
	for ( int i=0; i<DIM; ++i ) {
		norm_a  += a[i]*a[i];
		norm_b  += b[i]*b[i];
		a_dot_b += a[i]*b[i];
	}
	norm_a = sqrt(norm_a);
	norm_b = sqrt(norm_b);

	double theta = std::acos(a_dot_b/(norm_a*norm_b));

	return theta;
}

