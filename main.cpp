
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <exception>
#include <fstream>
#include <iostream>
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
	double l = 1.0; // [nm]

	// Derived geometric constants
	double l_xy = l*sqrt(2.0/3.0*(1.0 - cos(theta)));
	double dx = sqrt(3.0)/2.0*l_xy;
	double dy = l_xy/2.0;
	double dz = 0.5*sqrt(l*l - l_xy*l_xy);

	//--- Construct the HCP unit cell for the oxygens ---//
	typedef double rvec[3];
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
//	int numWaters = waterGrid[0]*waterGrid[1]*waterGrid[2];
//	int unitCellGrid[3] = { waterGrid[0]/2, waterGrid[1]/4, waterGrid[2]/2 };

	int unitCellGrid[3] = { 2, 1, 1 };
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

	//----- Write .gro file -----//
	// - Note: .gro file indexing starts with 1!

	// Include virtual sites
	FILE* pGroFile;
	pGroFile = fopen("outconf.gro", "w");
	fprintf( pGroFile, "Oxygen HCP lattice \n");
	fprintf( pGroFile, "    %d\n", numWaters );

	atomCounter = 0;
	for ( int i=0; i<numWaters; ++i )
	{
		fprintf( pGroFile, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
							i+1, "ICE", "OW", i+1,
						    x_HCP[i][0], x_HCP[i][1], x_HCP[i][2] );
	}
	fprintf( pGroFile, "   %2.5f   %2.5f   %2.5f\n", 
					   boxL[0], boxL[1], boxL[2] );

	fclose(pGroFile);

	//----- Cleanup -----//
	free(x_HCP);

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
