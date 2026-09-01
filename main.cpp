#include <iostream>
#include <fstream>
#include <chrono>
#include "delaunay.h"

typedef TetMesh_t<basicVec3d> TetMesh;

bool loadPoints(const char* filename, std::vector<double>& coords) {
	std::ifstream file(filename);

	if (!file.is_open()) {
		std::cerr << "Error: cannot open file '" << filename << "'." << std::endl;
		return false;
	}

	double x, y, z;
	while (file >> x >> y >> z) {
		coords.push_back(x);
		coords.push_back(y);
		coords.push_back(z);
	}

	if (!file.eof() && file.fail()) {
		std::cerr << "Error: invalid file format." << std::endl;
		return false;
	}

	return true;
}

// saveOutputFile
// 
// 'tin' is a tet mesh.
// 'filename' is the name of the output file without extension.
// The file produced will be called 'filename.tet' and/or
// 'filename.off' (if 's' option is used).
// 'options' is a (possibly empty) string of characters, each controlling
// one option as follows:
// n: binary output
// s: saves skin to an ASCII OFF file (outer convex hull)
// m: saves mesh to MEDIT format instead of TET

bool saveOutputFile(TetMesh& tin, const char* filename, const char* options) {
	bool binary = false, skin = false, medit = false;
	for (int i = 0; i < strlen(options); i++) switch (options[i]) {
	case 'n':
		binary = true; break;
	case 's':
		skin = true; break;
	case 'm':
		medit = true; break;
	}

	char tetfilename[2048], offfilename[2048];

	bool ret = true;

	if (medit) {
		sprintf(tetfilename, "%s.mesh", filename);
		ret &= tin.saveMEDIT(tetfilename, true);
	}
	else {
		sprintf(tetfilename, "%s.tet", filename);
		if (binary) ret &= tin.saveBinaryTET(tetfilename);
		else ret &= tin.saveTET(tetfilename, true);
	}

	if (skin) {
		sprintf(offfilename, "%s.off", filename);
		ret &= tin.saveBoundaryToOFF(offfilename);
	}

	return ret;
}

#ifdef _MSC_VER
#ifndef NDEBUG
#define DEBUG
#endif
#endif

int main(int argc, char* argv[])
{
#ifndef DEBUG
	if (argc < 2) {
		std::cout << "Delaunay3d - Create a Delaunay tetrahedrization out of a set of 3D points.\n";
		std::cout << "USAGE: Delaunay3d [-lbvfqnrs] filename.txt\n";
		std::cout << "Example: Delaunay3d -v test.txt\n";
		std::cout << "OPTIONS:\n";
		std::cout << "-v: verbose mode\n";
		std::cout << "-n: binary output\n";
		std::cout << "-m: use MEDIT format instead of TET\n";
		std::cout << "-s: saves skin to an ASCII OFF file (outer convex hull)\n";
		std::cout << "OUTPUT:\n";
		std::cout << "Output has same name (and path) as input with an extension appended.\n";
		std::cout << "E.g. Delaunay3d my_dir/test.txt produces my_dir/test.txt.tet\n";
		std::cout << "E.g. Delaunay3d -s my_dir/test.txt produces my_dir/test.txt.tet and my_dir/test.txt.off\n";
		return 0;
	}
#endif

	char filename[2048] = "..\\Input_file\\boeing_part.txt";

	std::string options = "";
#ifdef DEBUG
	options += "v";
#endif

	for (int i = 1; i < argc; i++)
		if (argv[i][0] == '-') {
			for (int j = 1; j < strlen(argv[i]); j++) options += argv[i][j];
		}
		else memcpy(filename, argv[i], strlen(argv[i])+1);

	// Load input file
	std::vector<double> coordinates;
	if (loadPoints(filename, coordinates)) {
		TetMesh* tin = new TetMesh;

		std::chrono::steady_clock::time_point time_point = std::chrono::steady_clock::now();
		tin->init_vertices(coordinates.data(), (uint32_t)coordinates.size() / 3);
		tin->tetrahedrize();
		tin->markInnerTets();
		std::chrono::steady_clock::time_point time_point2 = std::chrono::steady_clock::now();
		uint64_t ms = std::chrono::duration_cast<std::chrono::milliseconds>(time_point2 - time_point).count();
		std::cout << "Meshing terminated in " << ms << " milliseconds\n";

		if (saveOutputFile(*tin, filename, options.c_str()))
			printf("Finished\n");
	}

	return 0;
}
