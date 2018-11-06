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
// - Fix issues placing H's with unit_cell_grid = {1, 1, 1} (possible with this algorithm?)
// - Fix fixed-format output when #atoms or #waters is >5 char long in decimal format
// - Change atom names to "ice" (e.g. OI instead of OW)

#include "main.h"

int main(int argc, char* argv[])
{
	//----- Constants -----//

	// Water geometry: TIP4P/Ice
	const double PI = 3.14159265358979323846;
	const double theta = 109.47 * 2.0*PI/360.0; 	  // tetrahedral bond angle [radians]
	const double theta_HOH = 104.52 * 2.0*PI/360.0; // H-O-H bond angle [radians]
	const double l = 0.275;         // distance between nearest-neighbor oxygens [nm]
	                          // - About the same for I_h and I_c
	const double l_OH = 0.09572; // length of O-H bonds [nm]
	const double l_OM = 0.01577; // O-M distance [nm]

	// Distance between H's [nm]
	const double l_HH = l_OH*sqrt(2.0*(1.0 - cos(theta_HOH))); // ~0.15139 nm

	// Partial charges
	const double q_H = 0.520;
	const double q_M = -2.0*q_H;
	const std::vector<double> charges = { 0.0, q_H, q_H, q_M };
	int   num_atoms_per_water = charges.size();

	// Length of the projection of the O-H bond onto the H-O-H unit bisector
	const double l_proj = sqrt(l_OH*l_OH - l_HH*l_HH/4.0);
	
	// Hydrogen Bond length (from H to coordinated O) [nm]
	const double l_HB = sqrt( (l - l_proj)*(l - l_proj) + l_HH*l_HH/4.0 );


	//----- Input -----//

	// Check a configuration stored as a gro file
	if ( argc >= 3 and std::string(argv[1]) == "check" ) {
		// Check structure
		std::cout << "  IceBox: Checking configuration." << std::endl;

		// Read frame
		std::string gro_file( argv[2] );
		std::vector<Real3> coords;
		std::vector<std::string> atomTypes;
		std::vector<int> atomSerials;
		Real3 boxL;
		readFrameFromGroFile(gro_file, coords, atomTypes, atomSerials, boxL);

		Real3 dipole;
		double mu;
		double tol = 0.05;
		if ( argc == 4 ) {
			tol = std::stod( argv[3] );
		}
		check_structure( tol, coords, boxL, charges, l_OH, l_HH, l_OM, l_HB, l, theta,
										 dipole, mu );

		// Print dipole
		std::cout << "  IceBox: Dipole is ( "
								<< dipole[0] << ", " << dipole[1] << ", " << dipole[2] << " ) D "
							<< " (norm: " << mu << " D total\n";
		return 0;
	}

	// Usual operation: make a box of ice
	// TODO Input number of waters or something like that so the user can specify a kind of scale?
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
	int unit_cell_grid[DIM];
	unit_cell_grid[0] = std::atoi( argv[2] );
	unit_cell_grid[1] = std::atoi( argv[3] );
	unit_cell_grid[2] = std::atoi( argv[4] );



	//----- Construct the unit cell for the oxygens -----//

	std::vector<Real3> unit_cell;      // Unit cell oxygen positions
	Real3              unit_cellBoxL;  // Unit cell box size lengths
	int num_oxygens_per_unit_cell;

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
		num_oxygens_per_unit_cell = 8;
		unit_cell.resize(num_oxygens_per_unit_cell);

		// Unit cell box lengths
		unit_cellBoxL[0] = 2.0*dx;
		unit_cellBoxL[1] = 2.0*(l_xy + dy);
		unit_cellBoxL[2] = 2.0*(l + 2.0*dz);

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
		unit_cell[0][0] = 0.0;	
		unit_cell[0][1] = 0.0;	
		unit_cell[0][2] = 0.0;

		unit_cell[1][0] = unit_cell[0][0];	
		unit_cell[1][1] = unit_cell[0][1] + l_xy;
		unit_cell[1][2] = unit_cell[0][2];
	
		unit_cell[2][0] = unit_cell[1][0] + dx;
		unit_cell[2][1] = unit_cell[1][1] + dy;
		unit_cell[2][2] = unit_cell[1][2];

		unit_cell[3][0] = unit_cell[2][0];
		unit_cell[3][1] = unit_cell[2][1] + l_xy;
		unit_cell[3][2] = unit_cell[2][2];

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
				unit_cell[i + 4][j] = unit_cell[i][j];
			}
			unit_cell[i + 4][2] += l + 2.0*dz; // spacing between planes
		}

		// Shift all atoms up/down relative to the plane (see diagrams above)
		for ( int i=0; i<num_oxygens_per_unit_cell; ++i ) {
			// Shift up these atoms
			if ( (i == 0) || (i == 2) || (i == 5) || (i == 7 ) ) {
				unit_cell[i][2] += dz;
			}
			// Shift down the rest
			else {
				unit_cell[i][2] -= dz;
			}
		}

		// Offset all atoms so that none lie on the edges of the unit cell box
		Real3 x_shift = {{ dx/2.0, dy/2.0, l/2.0 }};
		for ( int i=0; i<num_oxygens_per_unit_cell; ++i ) {
			for ( int j=0; j<DIM; ++j ) {
				unit_cell[i][j] += x_shift[j];
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
		num_oxygens_per_unit_cell = 8;
		unit_cell.resize(num_oxygens_per_unit_cell);

		// Unit cell box lengths
		unit_cellBoxL = {{ a, a, a }};

		// First, place the atoms which fall on a regular FCC lattice (easy)
		unit_cell[0] = {{ 0.0,   0.0,   0.0   }};
		unit_cell[1] = {{ a/2.0, 0.0,   a/2.0 }};
		unit_cell[2] = {{ 0.0,   a/2.0, a/2.0 }};
		unit_cell[3] = {{ a/2.0, a/2.0, a     }};

		//----- Now place the 4 non-FCC oxygens -----//

		Real3 x_center = {{ a/2.0, a/2.0, a/2.0 }}; // center of unit cell
		Real3 x_corner; // opposing corner
		int   index;    // index in unit cell

		// Top-front
		x_corner = {{ 0.0, 0.0, a }};
		index = 4;
		for ( int j=0; j<DIM; ++j ) {
			unit_cell[index][j] = x_center[j] + 0.5*(x_corner[j] - x_center[d]);
		}

		// Top-back
		x_corner = {{ a, a, a }};
		index = 5;
		for ( int j=0; j<DIM; ++j ) {
			unit_cell[index][j] = x_center[j] + 0.5*(x_corner[j] - x_center[d]);
		}

		// Bottom-right
		x_corner = {{ a, 0.0, 0.0 }};
		index = 6;
		for ( int j=0; j<DIM; ++j ) {
			unit_cell[index][j] = x_center[j] + 0.5*(x_corner[j] - x_center[d]);
		}

		// Bottom-left
		x_corner = {{ 0.0, a, 0.0 }};
		index = 7;
		for ( int j=0; j<DIM; ++j ) {
			unit_cell[index][j] = x_center[j] + 0.5*(x_corner[j] - x_center[d]);
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
	int numUnitCells = unit_cell_grid[0]*unit_cell_grid[1]*unit_cell_grid[2];
	int num_waters = num_oxygens_per_unit_cell*numUnitCells;

	std::cout << "  IceBox: The lattice will have " << num_waters << " waters." << std::endl;

	// Replicate unit cell
	std::vector<Real3> x_O_lattice(num_waters);
	int atom_counter = 0;
	for ( int i=0; i<unit_cell_grid[0]; ++i ) { // x cells
		for ( int j=0; j<unit_cell_grid[1]; ++j ) { // y cells
			for ( int k=0; k<unit_cell_grid[2]; ++k ) { // z cells
				for ( int m=0; m<num_oxygens_per_unit_cell; ++m ) { // unit cell atoms
					// x, y, z dimensions
					x_O_lattice[atom_counter][0] = unit_cell[m][0] + i*unit_cellBoxL[0];
					x_O_lattice[atom_counter][1] = unit_cell[m][1] + j*unit_cellBoxL[1];
					x_O_lattice[atom_counter][2] = unit_cell[m][2] + k*unit_cellBoxL[2];

					++atom_counter;
				}
			}
		}
	}

	// Simulation box lengths
	Real3 boxL;
	for ( int j=0; j<DIM; ++j ) {
		boxL[j] = unit_cell_grid[j]*unit_cellBoxL[j];
	}

	// Density
	double molar_mass = 18.02;        // g/mol
	const double N_Av = 6.022e23;     // Avogadro's number (#/mol)
	const double cm_per_nm = 1.0e-7;  // conversion factor
	double box_volume = 1.0;          // in cm&3
	for ( int d=0; d<DIM; ++d ) {
		box_volume *= boxL[d]*cm_per_nm;
	}
	//              (molecules/cm^3)         g/(molecule)
	double rho = (num_waters/box_volume)*(molar_mass/N_Av);

	std::cout << "  Icebox: Density is " << rho << " g/cm^3.\n";


	//----- Place hydrogens using a Monte Carlo routine -----//
	// Method: Buch et al., J. Phys. Chem. B 102.44 (1998)

	std::cout << "  Icebox: Finding all O-O nearest-neighbor pairs" << "\n";

	// Set up RNG engine
	std::seed_seq seed_sequence = { 749725171 };
	/*
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::seed_seq seed_sequence = { seed };
	*/
	std::mt19937 rng(seed_sequence);	// Mersenne Twister

	// Find all O-O nearest-neighbor pairs, and randomly assign an H to one atom in each pair
	int   numH = 2*num_waters;
	int   numPairs = numH;
	std::vector<Pair> pairs(numPairs);  // NOTE: index over pairs = index over H's

	std::vector<int>  hydrogen_counts(num_waters, 0); 		// how many bound H's each oxygen has
	std::vector<int>  hydrogenOwners(numH, 0);	   			// which O has the H
	std::vector<bool> doesFirstAtomHaveHydrogen(numPairs);	// Does pair[][0] have the H?

	std::uniform_int_distribution<int> coin_flip(0, 1);

	// Working variables
	int    pair_index = 0;
	double dist, dist_sq;
	double pair_dist = 1.001*l; // Distance between nearest-neighbor oxygens
	Real3   x_i_j;              // Direction: i --> j

	for ( int i=0; i<num_waters; ++i ) {
		for ( int j=i+1; j<num_waters; ++j ) {
			// Minimum image vector: i --> j
			minImage(x_O_lattice[i], x_O_lattice[j], boxL, x_i_j, dist_sq);
			dist = sqrt(dist_sq); // = || x_i_j ||_2

			if ( dist <= pair_dist ) {
				// Record the pair
				pairs[pair_index][0] = i;
				pairs[pair_index][1] = j;

				// Randomly assign the H to one of the paired atoms
				if ( coin_flip(rng) == 0 ) {
					// Give to 'i'
					++(hydrogen_counts[i]);
					hydrogenOwners[pair_index] = i;
					doesFirstAtomHaveHydrogen[pair_index] = true;
				}
				else {
					// Give to 'j'
					++(hydrogen_counts[j]);
					hydrogenOwners[pair_index] = j;
					doesFirstAtomHaveHydrogen[pair_index] = false;
				}

				++pair_index;
			}
		}

		// numPairs == numHydrogens
		if ( pair_index == numPairs ) {
			break;
		}
	}

	// Stochastically reassign hydrogens until all O's have two H's each
	// - For each pair, the owner is the O closest to the H; the other O is its partner
	std::cout << "  IceBox: Stochastically assigning hydrogen atoms." << std::endl;
	int delta, delta_trial, owner, partner;
	bool is_swap_accepted;
	std::uniform_int_distribution<int> random_pair_generator(0, numPairs - 1);
	while ( areHydrogensCorrectlyPlaced(hydrogen_counts) == false ) {
		pair_index = random_pair_generator(rng);
		owner = hydrogenOwners[pair_index];

		// Where are the owner and partner in the pairs array?
		if ( pairs[pair_index][0] == owner ) {
			partner = pairs[pair_index][1];
		}
		else {
			partner = pairs[pair_index][0];
		}

		// Current difference 
		delta = abs(hydrogen_counts[owner] - hydrogen_counts[partner]);

		// Result of proposed swap
		delta_trial = abs( (hydrogen_counts[owner] - 1) - (hydrogen_counts[partner] + 1) );

		// Determine whether to accept the proposed swap
		if ( delta_trial < delta ) {
			is_swap_accepted = true;
		}
		else if ( delta_trial == delta ) {
			// 50/50 chance to accept anyway
			is_swap_accepted = static_cast<bool>( coin_flip(rng) );
		}
		else {
			is_swap_accepted = false;
		}

		// Make any necessary changes
		if ( is_swap_accepted ) {
			// NOTE: owner/partner here refer to the roles in the "old" arrangement

			hydrogenOwners[pair_index] = partner;

			--(hydrogen_counts[owner]);
			++(hydrogen_counts[partner]);

			if ( doesFirstAtomHaveHydrogen[pair_index] == true ) {
				doesFirstAtomHaveHydrogen[pair_index] = false;
			}
			else {
				doesFirstAtomHaveHydrogen[pair_index] = true;
			}
		}
	}

	//----- Construct the full coordinates array -----//

	std::cout << "  IceBox: Constructing full array of atomic positions." << std::endl;
	int    numAtoms = 4*num_waters, count, atomIndex;
	double norm_x;
	Real3   zero3 = {{ 0.0, 0.0, 0.0 }};
	Real3   x_OW, x_MW, x_HW1, x_HW2;
	Real3   x_OO_1 = zero3, x_OO_2 = zero3, unitBisector, unitNormal; // working vectors
	std::vector<Real3> x_all(numAtoms);

	for ( int i=0; i<num_waters; ++i ) {
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
					minImage(x_O_lattice[owner], x_O_lattice[partner], boxL, x_OO_1, dist_sq);
					++count;
				}
				else if ( count == 1 ) {
					minImage(x_O_lattice[owner], x_O_lattice[partner], boxL, x_OO_2, dist_sq);
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

	Real3 dipole;
	double mu;
	double tol = 0.01;
	check_structure( tol, x_all, boxL, charges, l_OH, l_HH, l_OM, l_HB, l, theta,
	                 dipole, mu );

	// Print dipole
	std::cout << "  IceBox: Dipole is ( "
              << dipole[0] << ", " << dipole[1] << ", " << dipole[2] << " ) D "
 	          << " (norm: " << mu << " D total, " 
	          << mu/numUnitCells << " D per unit cell)" << "\n";


	//----- Write .gro file -----//
	// - Note: .gro file indexing starts with 1!

	std::cout << "  IceBox: Writing .gro file." << std::endl;

	// File name
	std::ostringstream ossGroFileName;
	ossGroFileName << "initconf_ice_" << phase << "_" << num_waters << ".gro";
	std::string gro_fileName( ossGroFileName.str() );

	FILE* gro_file_ptr;
	gro_file_ptr = fopen(gro_fileName.c_str(), "w");
	fprintf( gro_file_ptr, "TIP4P/Ice in a box\n" );
	fprintf( gro_file_ptr, "    %d\n", numAtoms );

	for ( int i=0; i<num_waters; ++i ) {
		// OW
		atomIndex = 4*i;	   // Index in the coordinates array "x_all"

		fprintf( gro_file_ptr, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
							i+1, "ICE", "OI", atomIndex+1,
						    x_all[atomIndex][0], x_all[atomIndex][1], x_all[atomIndex][2] );

		// HW1
		++atomIndex;

		fprintf( gro_file_ptr, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
							i+1, "ICE", "HI1", atomIndex+1,
						    x_all[atomIndex][0], x_all[atomIndex][1], x_all[atomIndex][2] );

		// HW2
		++atomIndex;

		fprintf( gro_file_ptr, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
							i+1, "ICE", "HI2", atomIndex+1,
							x_all[atomIndex][0], x_all[atomIndex][1], x_all[atomIndex][2] );

		// MW
		++atomIndex;

		fprintf( gro_file_ptr, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
							i+1, "ICE", "MI", atomIndex+1,
							x_all[atomIndex][0], x_all[atomIndex][1], x_all[atomIndex][2] );

	}
	fprintf( gro_file_ptr, "   %2.5f   %2.5f   %2.5f\n", 
					   boxL[0], boxL[1], boxL[2] );

	fclose(gro_file_ptr);

	std::cout << "  IceBox: Done.\n";
}



// Applies the minimum image convention
void minImage(const Real3& x1, const Real3& x2, const Real3& boxL, Real3& x12, double& dist_sq)
{
	dist_sq = 0.0;

	for ( int d=0; d<DIM; d++ ) {
		x12[d] = x2[d] - x1[d];
		// Apply minimum image convention
		if      ( x12[d] >  0.5*boxL[d] ) { x12[d] -= boxL[d]; }
		else if ( x12[d] < -0.5*boxL[d] ) { x12[d] += boxL[d]; }

		dist_sq += x12[d]*x12[d];
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
bool areHydrogensCorrectlyPlaced(const std::vector<int>& hydrogen_counts)
{
	int num_waters = hydrogen_counts.size();
	for ( int i=0; i<num_waters; i++ ) {
		if ( hydrogen_counts[i] != 2 ) {
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


// Read a single frame from a .gro file 
void readFrameFromGroFile(
		const std::string& gro_file,
		// Output
		std::vector<Real3>& coords, std::vector<std::string>& atomTypes,
		std::vector<int>& atomSerials, Real3& boxL)
{
	// Get the number of atoms for a consistency check vs. the xtc file
	std::ifstream ifs(gro_file);
	if ( ! ifs.is_open() ) {
		std::cerr << "readFrameFromGroFile - Unable to open .gro file.\n"
							<< "(input: " << gro_file << ").\n";
		exit(1);
	}

	std::string        line;
	std::istringstream ss;

	int numAtoms = -1;
	getline(ifs, line);		// Ignore first line
	getline(ifs, line);		// Second line contains the number of atoms 
	ss.str(line);
	ss >> numAtoms;

	// Checks
	if ( numAtoms < 0 ) {
		std::cout << "  XdrFileTools::readFrameFromGroFile - Invalid number of atoms "
                  << "read from .gro file." << "\n";
		exit(1);
	}

	// Allocate memory
	coords.resize(numAtoms);
	atomTypes.resize(numAtoms);
	atomSerials.resize(numAtoms);

	// Continue reading
	int               atom_counter = 0, moleculeIndex;
	std::string       moleculeName, atomName;
	std::stringstream parsingBuffer;

	while ( (atom_counter<numAtoms) && getline(ifs, line) ) {
		// Parse line
		ss.str(line);

		// GROMACS FORMAT: "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
		// - Number of characters per field: 5, 5, 5, 5, 8, 8, 8, 8, 8, 8
		// - Use a stringstream to neatly handle the different variable types

		// Clear the parsing buffer
		parsingBuffer.str( std::string() );
		parsingBuffer.clear();

		// Field 1: molecule index (0-4, 5 char)
		parsingBuffer << line.substr(0, 5) << "\t";

		// Field 2: molecule name (5-9, 5 char)
		parsingBuffer << line.substr(5, 5) << "\t";

		// Field 3: atom name (10-14, 5 char)
		parsingBuffer << line.substr(10, 5) << "\t";

		// Field 4: atom index (15-19, 5 char)
		parsingBuffer << line.substr(15, 5) << "\t";

		// Field 5: x-position (20-27, 8 char)
		parsingBuffer << line.substr(20, 8) << "\t";

		// Field 6: y-position (28-35, 8 char)
		parsingBuffer << line.substr(28, 8) << "\t";

		// Field 7: z-position (36-43, 8 char)
		parsingBuffer << line.substr(36, 8) << "\t";

		// Store variables
		parsingBuffer >> moleculeIndex;
		parsingBuffer >> moleculeName;
		parsingBuffer >> atomName;
		parsingBuffer >> atomSerials[atom_counter];
		for ( int d=0; d<DIM; ++d ) {
			parsingBuffer >> coords[atom_counter][d];
		}
		atomTypes[atom_counter] = atomName;

		atom_counter++;
	}

	//----- Last line has the box lengths -----//
	getline(ifs, line);

	// Clear flags before using with a new line
	ss.clear();
	ss.str(line);

	// Clear the parsing buffer
	parsingBuffer.str("");
	parsingBuffer.clear();

	for ( int i=0; i<3; ++i ) {
		ss >> boxL[i];
	}

	// Close .gro file
	ifs.close();

	// Checks 
	for ( int i=0; i<3; ++i ) {
		if ( boxL[i] <= 0.0 ) {
			std::cerr << "readFrameFromGroFile - Invalid box length for dimension " 
                      << i+1 << " of 3 (boxL = " << boxL[i] << ")." << "\n";
			exit(1);
		}
	}

	return;
}


void check_structure(
	const double tol,  // tolerance (%err/100)
	// Frame
	const std::vector<Real3>& x_all,
	const Real3& boxL,
	// Model parameters
	const std::vector<double>& charges,
	const double l_OH, const double l_HH, const double l_OM,
	// Packing structure
	const double l_HB, const double r_nn, const double theta,
	// Output
	Real3& dipole, double& mu
)
{
	int num_atoms_per_water = charges.size();
	int num_waters = static_cast<int>(x_all.size())/num_atoms_per_water;

	// Check the internal structure of each molecule (use 1% margin of error)
	int atom_index;
	Real3 x_i_j;
	double dist, dist_sq;
	for ( int i=0; i<num_waters; ++i ) {
		// Index of the oxygen
		atom_index = num_atoms_per_water*i;

		// O-H1 distance
		minImage(x_all[atom_index], x_all[atom_index+1], boxL, x_i_j, dist_sq);
		dist = sqrt(dist_sq);
		if ( dist > (1.0 + tol)*l_OH ) {
			std::cerr << "  Molecule " << i << ": OW-HW1 bond is " << dist << " nm " 
			          << "(should be " << l_OH << " nm)." << std::endl;
		}

		// O-H2 distance
		minImage(x_all[atom_index], x_all[atom_index+2], boxL, x_i_j, dist_sq);
		dist = sqrt(dist_sq);
		if ( dist > (1.0 + tol)*l_OH ) {
			std::cerr << "  Molecule " << i << ": OW-HW2 bond is " << dist << " nm "
			          << "(should be " << l_OH << " nm)." << std::endl;
		}

		// H1-H2 distance
		minImage(x_all[atom_index+1], x_all[atom_index+2], boxL, x_i_j, dist_sq);
		dist = sqrt(dist_sq);
		if ( dist > (1.0 + tol)*l_HH ) {
			std::cerr << "  Molecule " << i << ": HW1-HW2 bond is " << dist << " nm "
			          << "(should be " << l_HH << " nm)." << std::endl;
		}

		// O-M distance
		/* FIXME tolerance too small for reading at gro file precision?
		minImage(x_all[atom_index], x_all[atom_index+3], boxL, x_i_j, dist_sq);
		dist = sqrt(dist_sq);
		if ( dist > (1.0 + tol)*l_OM ) {
			std::cerr << "  Molecule " << i << ": OW-MW bond is " << dist << " nm "
			          << "(should be " << l_OM << " nm)." << std::endl;
		}
		*/
	}

	// Neighbor list for oxygen atoms
	double pair_dist = (1.0 + tol)*r_nn;
	std::vector<std::vector<int>> neighbor_list(num_waters);
	for ( int i=0; i<num_waters; ++i ) {
		for ( int j=i+1; j<num_waters; ++j ) {
			minImage(x_all[num_atoms_per_water*i], x_all[num_atoms_per_water*j], boxL, x_i_j, dist_sq);
			dist = sqrt(dist_sq);
			if ( dist <= pair_dist ) {
				neighbor_list[i].push_back(j);
				neighbor_list[j].push_back(i);
			}
		}
	}

	// Check that all oxygens are tetrahedrally coordinated
	std::cout << "  IceBox: Checking oxygen lattice for tetrahedral geometry "
	          << "and number of nearest-neighbor oxygens." << "\n";

	Real3  x_i_k;
	int    num_neighbors;
	double theta_high = (1.0 + tol)*theta;
	double theta_low  = (1.0 - tol)*theta;
	double theta_ijk;
	const double PI = 3.14159265358979323846;
	const double deg_per_rad = 180.0/PI;
	for ( int i=0; i<num_waters; ++i ) {
		num_neighbors = neighbor_list[i].size();
		if ( num_neighbors != 4 ) {
			std::cerr << "  IceBox: Atom " << i << " in the oxygen lattice has "
			          << num_neighbors << " neighbors (should have 4)!" << "\n";
			exit(1);
		}

		// Choose 'j' as the first neighbor of 'i' in the list
		int j = neighbor_list[i][0];
		minImage(x_all[4*i], x_all[4*j], boxL, x_i_j, dist_sq);

		for ( int n=1; n<num_neighbors; ++n ) {
			int k = neighbor_list[i][n];
			minImage(x_all[4*i], x_all[4*k], boxL, x_i_k, dist_sq);

			theta_ijk = angleBetweenVectors(x_i_j, x_i_k);
			if ( (theta_ijk < theta_low) || (theta_high < theta_ijk) ) {
				std::cerr << "  IceBox: Angle " << 4*i << "-" << 4*j << "-" << 4*k
                          << " (between neighboring oxygen atoms) is " << theta_ijk*deg_per_rad
				          << " degrees (should be " << theta*deg_per_rad << " degrees)!" << "\n";
				//exit(1);
			}
		}
	}

	// Check the hydrogens that are coordinated with (but not bonded to) each oxygen 
	int num_coordinated_hydrogens;
	pair_dist = (1.0 + tol)*l_HB;
	for ( int i=0; i<num_waters; ++i ) {
		num_coordinated_hydrogens = 0;

		for ( int j=0; j<num_waters; ++j ) {
			if ( i != j ) {
				// O of water 'i' and H_k of water 'j'
				for ( int k=1; k<=2; ++k ) {
					minImage(x_all[4*i], x_all[4*j + k], boxL, x_i_k, dist_sq);
					dist = sqrt(dist_sq);

					if ( (dist <= pair_dist) ) { //&& (dist >= 0.999*l_OH) )
						++num_coordinated_hydrogens;
					}
				}
			}
		}

		// Check
		if ( num_coordinated_hydrogens != 2 ) {
			std::cerr << "  IceBox: Water " << i << " is coordinated with "
			          << num_coordinated_hydrogens << " hydrogens that are not "
			          << "bonded to it (should be 2!)." << "\n";
			exit(1);
		}
	}

	//-----  Compute the net dipole -----//

	dipole.fill(0.0);
	for ( int i=0; i<num_waters; ++i ) {
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
	mu = 0.0;
	for ( int k=0; k<DIM; ++k ) { mu += dipole[k]*dipole[k]; }
	mu = sqrt(mu);
}
