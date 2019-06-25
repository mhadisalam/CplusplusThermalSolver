//LLPD Computational Two-Fluid Plasma Dynamics Thesis
//Reading OpenFOAM polyMesh + Two-Volume Thermal Conduction Solver
//M. Hadi Salam, American University of Beirut, Mechanical Engineering Department
//Sunday, February 17th, 2019

//Compiling: g++ TwoThermal.cpp -o TwoThermal.bin -pthread -std=c++11

//pre-processor directives

#include <iostream>
#include <fstream>
#include <string>
#include <new>
#include <thread>
#include <cstdlib>

//using statement

using namespace std;

//Function declarations

void Create0Dir()
{
	const int dir_err = system("mkdir 0");

	if(dir_err==-1)
	{
		cout << "\nError creating directory.\n" << endl;
		exit(1);
	}
}

void CreateDir(const char* cmd)
{
	const int dir_err = system(cmd);

	if(dir_err==-1)
	{
		cout << "\nError creating directory.\n" << endl;
		exit(1);
	}
}

//Function Main

int main()
{
	unsigned long int nPoints = 0;

	string dummyStr;

	char dummyChar;

	//double test = 0;

	double * pointsX;

	double * pointsY;

	double * pointsZ;

	unsigned long int index = 0;

	cout.precision(15);

	//Read Points

	cout << "\n" << "Reading polyMesh points file..." << endl;

	ifstream pointsFile ("./constant/polyMesh/points");

	//skip file header to line 19 containing number of points

	for(int counter1=1; counter1<=18; counter1++)
	{
		getline(pointsFile,dummyStr);
	}

	pointsFile >> nPoints;

	cout << "\nThere are " << nPoints << " points, which are the cell vertices in this OpenFOAM polyMesh.\n" << endl;

	pointsX = new double [nPoints];

	pointsY = new double [nPoints];

	pointsZ = new double [nPoints];

	getline(pointsFile,dummyStr);
	getline(pointsFile,dummyStr);

	//For Loop to read the Points' Coordinates

	for(unsigned long int pointsCounter=1; pointsCounter<=nPoints; pointsCounter++)
	{
		index = pointsCounter-1;
		pointsFile >> dummyChar;
		pointsFile >> pointsX[index];
		//pointsFile >> dummyChar;
		pointsFile >> pointsY[index];
		//pointsFile >> dummyChar;
		pointsFile >> pointsZ[index];
		pointsFile >> dummyChar;
		//pointsFile >> dummyChar;
	}

	//Close Points File

	pointsFile.close();

	//Show Read Points data

	cout << "Displaying Read Points Data as Point# X Y Z:\n" << endl;

	for(unsigned long int displayCounter=1; displayCounter<=nPoints; displayCounter++)
	{
		index = displayCounter-1;

		cout << "Point" << index << " " << pointsX[index] << " " << pointsY[index] << " " << pointsZ[index] << endl;
	}

	//Read Faces

	cout << "\n" << "Reading polyMesh faces file..." << endl;

	ifstream facesFile ("./constant/polyMesh/faces");

	//skip file header to line 19 containing number of faces

	for(int counter2=1; counter2<=18; counter2++)
	{
		getline(facesFile,dummyStr);
	}

	unsigned long int nFaces=0;

	facesFile >> nFaces;

	cout << "\nThere are " << nFaces << " faces in this OpenFOAM polyMesh.\n" << endl;

	int PointsPerFace=0;

	getline(facesFile,dummyStr);
	getline(facesFile,dummyStr);

	facesFile >> PointsPerFace;

	cout << "\nThere are " << PointsPerFace << " points per face in this OpenFOAM polyMesh faces file.\n" << endl;

	unsigned long int ** FacePoints = new unsigned long int *[nFaces];

	for(unsigned long int i = 0; i < nFaces; i++)
	{
		FacePoints[i] = new unsigned long int[PointsPerFace];
	}

	//FacePoints = new unsigned long int [nFaces][4];

	facesFile >> dummyChar;

	facesFile >> FacePoints[0][0];

	facesFile >> FacePoints[0][1];

	facesFile >> FacePoints[0][2];

	facesFile >> FacePoints[0][3];

	facesFile >> dummyChar;

	cout <<"\nFace 0 has points " << FacePoints[0][0] << " " << FacePoints[0][1] << " " << FacePoints[0][2] << " " << FacePoints[0][3] << "\n" << endl;

	int dummyInt=0;

	for(unsigned long int facesCounter=2; facesCounter<=nFaces; facesCounter++)
	{
		index = facesCounter - 1;
		facesFile >> dummyInt;
		facesFile >> dummyChar;

		facesFile >> FacePoints[index][0];

		facesFile >> FacePoints[index][1];

		facesFile >> FacePoints[index][2];

		facesFile >> FacePoints[index][3];

		facesFile >> dummyChar;
	}

	facesFile.close();

	cout << "\nDisplaying Read Faces Data as Face# Point1 Point2 Point3 Point4:\n" << endl;

	for(unsigned long int dCounter=1; dCounter<=nFaces; dCounter++)
	{
		index = dCounter-1;

		cout << "Face" << index << "  " << FacePoints[index][0] << " " << FacePoints[index][1] << " " << FacePoints[index][2] << " " << FacePoints[index][3] << endl;
	}

	//pointsFile >> test;

	//cout << "\nThe first read decimal is " << test << " .\n" << endl;


	//Read Owners

	cout << "\n" << "Reading polyMesh owner file...\n" << endl;

	ifstream ownerFile ("./constant/polyMesh/owner");

	//skip file header to line 22 containing first Owner entry

	for(int counter3=1; counter3<=21; counter3++)
	{
		getline(ownerFile,dummyStr);
	}

	unsigned long int * FaceOwnerCell;

	FaceOwnerCell = new unsigned long int [nFaces];

	for(unsigned long int counter4=1; counter4<=nFaces; counter4++)
	{
		index = counter4-1;
		ownerFile >> FaceOwnerCell[index];
	}

	ownerFile.close();

	//Show Read Owner Data

	for(unsigned long int counter5=1; counter5<=nFaces; counter5++)
	{
		index = counter5-1;

		cout << "Face" << index << " is owned by Cell" << FaceOwnerCell[index] << endl;
	}

	//Read Neighbours

	cout << "\n" << "Reading polyMesh neighbour file...\n" << endl;

	ifstream neighbourFile ("./constant/polyMesh/neighbour");

	//skip file header to line 20 containing number of Internal Faces

	for(int counter6=1; counter6<=19; counter6++)
	{
		getline(neighbourFile,dummyStr);
	}

	unsigned long int nInternalFaces=0;

	neighbourFile >> nInternalFaces;

	cout << "\nThere are " << nInternalFaces << " Internal Faces.\n" << endl;

	getline(neighbourFile,dummyStr);
	getline(neighbourFile,dummyStr);

	unsigned long int * FaceNeighbourCell;

	FaceNeighbourCell = new unsigned long int [nInternalFaces];

	for(unsigned long int counter7=1; counter7<=nInternalFaces; counter7++)
	{
		index=counter7-1;
		neighbourFile >> FaceNeighbourCell[index];
	}

	neighbourFile.close();

	//Show Read Neighbour Data

	for(unsigned long int counter8=1; counter8<=nInternalFaces; counter8++)
	{
		index=counter8-1;

		cout << "Internal Face" << index << " has Neighbour Cell" << FaceNeighbourCell[index] << endl;
	}

	//Find IsBoundary Face

	bool * IsBoundaryFace;

	IsBoundaryFace = new bool [nFaces];

	for(unsigned long int counter9=1; counter9<=nFaces; counter9++)
	{
		index = counter9 - 1;

		IsBoundaryFace[index]=0;
	}

	for(unsigned long int c10=nInternalFaces; c10<=(nFaces-1); c10++)//Changed c10=(nInternalFaces-1) to c10=nInternalFaces to correct
	{
		IsBoundaryFace[c10]=1;
	}

	cout << endl;

	for(unsigned long int c11=0; c11<=(nFaces-1); c11++)
	{
		if(IsBoundaryFace[c11]==0)
		{
			cout << "Face " << c11 << " is an Internal Face.\n";
		}
		else
		{
			cout << "Face " << c11 << " is a Boundary Face.\n";
		}

		/*if(c11==10000)//Just to pause for testing, comment out in production
		{
			cin >> dummyChar;
		}*/
	}

	cout << endl;

	unsigned long int nCells=0;

	/*cout << "Enter number of Cells, nCells: ";

	cin >> nCells;*/

	//Calculate nCells, 6 Faces per Cell One-To-One Accounted For by nFaces and nInternalFaces

	nCells = (nFaces+nInternalFaces)/6;

	unsigned long int ** CellFaces = new unsigned long int *[nCells];

	for(unsigned long int c12 = 0; c12 < nCells; c12++)
	{
		CellFaces[c12] = new unsigned long int[6];
	}

	int * CellFaceIndex;

	CellFaceIndex = new int [nCells];

	for(unsigned long int c14=0; c14<nCells; c14++)
	{
		CellFaceIndex[c14]=0;
	}

	//Get Cell Faces from FaceOwnerCell List

	for(unsigned long int c15=0; c15<nFaces; c15++)
	{
		CellFaces[FaceOwnerCell[c15]][CellFaceIndex[FaceOwnerCell[c15]]]=c15;

		CellFaceIndex[FaceOwnerCell[c15]]+=1;
	}

	//Display Found Faces per Cell

	for(unsigned long int c16=0; c16<nCells; c16++)
	{
		cout << "Cell " << c16 << " has faces";

		for(int c17=0; c17<CellFaceIndex[c16]; c17++)
		{
			cout << " " << CellFaces[c16][c17];
		}

		cout << endl;
	}

	//Get Cell Faces from FaceNeighbourCell List

	for(unsigned long int c18=0; c18<nInternalFaces; c18++)
	{
		CellFaces[FaceNeighbourCell[c18]][CellFaceIndex[FaceNeighbourCell[c18]]]=c18;

		CellFaceIndex[FaceNeighbourCell[c18]]+=1;
	}

	//Display Total Found Faces per Cell

	for(unsigned long int c19=0; c19<nCells; c19++)
	{
		cout << "Cell " << c19 << " has faces";

		for(int c20=0; c20<CellFaceIndex[c19]; c20++)
		{
			cout << " " << CellFaces[c19][c20];
		}

		cout << endl;
	}

	//Check That All Cells Have 6 Faces

	cout << "\nVerifying that All Cells Have 6 Faces...\n" << endl;

	bool Test = 1;

	for(unsigned long int c21=0; c21<nCells; c21++)
	{
		if(CellFaceIndex[c21]!=6)
		{
			Test = 0;
		}
	}

	if(Test==1)
	{
		cout << "Test passed. All Cells have 6 Faces.\n" << endl;
	}
	else
	{
		cout << "Test failed.\n" << endl;
	}

	//Calculate Face Center Coordinates

	cout << "\nCalculating Face Center Coordinates...\n" << endl;

	double * FaceCenterX;

	double * FaceCenterY;

	double * FaceCenterZ;

	FaceCenterX = new double [nFaces];

	FaceCenterY = new double [nFaces];

	FaceCenterZ = new double [nFaces];

	for(unsigned long int c22=0; c22<nFaces; c22++)
	{
		FaceCenterX[c22] = ((pointsX[FacePoints[c22][0]])+(pointsX[FacePoints[c22][1]])+(pointsX[FacePoints[c22][2]])+(pointsX[FacePoints[c22][3]]))/4.0;

		FaceCenterY[c22] = ((pointsY[FacePoints[c22][0]])+(pointsY[FacePoints[c22][1]])+(pointsY[FacePoints[c22][2]])+(pointsY[FacePoints[c22][3]]))/4.0;

		FaceCenterZ[c22] = ((pointsZ[FacePoints[c22][0]])+(pointsZ[FacePoints[c22][1]])+(pointsZ[FacePoints[c22][2]])+(pointsZ[FacePoints[c22][3]]))/4.0;
	}

	cout << "Face Center Coordinates Calculated.\n" << endl;

	//Calculate Cell Center Coordinates

	cout << "\nCalculating Cell Center Coordinates...\n" << endl;

	double * CellCenterX;

	double * CellCenterY;

	double * CellCenterZ;

	CellCenterX = new double [nCells];

	CellCenterY = new double [nCells];

	CellCenterZ = new double [nCells];

	for(unsigned long int c23=0; c23<nCells; c23++)
	{
		CellCenterX[c23] = ((FaceCenterX[CellFaces[c23][0]])+(FaceCenterX[CellFaces[c23][1]])+(FaceCenterX[CellFaces[c23][2]])+(FaceCenterX[CellFaces[c23][3]])+(FaceCenterX[CellFaces[c23][4]])+(FaceCenterX[CellFaces[c23][5]]))/6.0;

		CellCenterY[c23] = ((FaceCenterY[CellFaces[c23][0]])+(FaceCenterY[CellFaces[c23][1]])+(FaceCenterY[CellFaces[c23][2]])+(FaceCenterY[CellFaces[c23][3]])+(FaceCenterY[CellFaces[c23][4]])+(FaceCenterY[CellFaces[c23][5]]))/6.0;

		CellCenterZ[c23] = ((FaceCenterZ[CellFaces[c23][0]])+(FaceCenterZ[CellFaces[c23][1]])+(FaceCenterZ[CellFaces[c23][2]])+(FaceCenterZ[CellFaces[c23][3]])+(FaceCenterZ[CellFaces[c23][4]])+(FaceCenterZ[CellFaces[c23][5]]))/6.0;
	}

	cout << "Cell Center Coordinates Calculated.\n" << endl;

	//Read Boundary file

	cout << "\nReading boundary file...\n" << endl;

	ifstream boundaryFile ("./constant/polyMesh/boundary");

	unsigned long int nInletFaces = 0;

	unsigned long int InletStartFace = 0;

	unsigned long int nOutletFaces = 0;

	unsigned long int OutletStartFace = 0;

	unsigned long int nWallsFaces = 0;

	unsigned long int WallsStartFace = 0;

	//Go to line 24

	for(int c24=1; c24<=23; c24++)
	{
		getline(boundaryFile,dummyStr);
	}

	//Skip 6 characters, >> 'Eats' the spaces.

	for(int c25=1; c25<=6; c25++)
	{
		boundaryFile >> dummyChar;
	}

	boundaryFile >> nInletFaces;

	cout << "There are " << nInletFaces << " Inlet Faces.\n" << endl;

	getline(boundaryFile,dummyStr);

	for(int c26=1; c26<=9; c26++)
	{
		boundaryFile >> dummyChar;
	}

	boundaryFile >> InletStartFace;

	cout << "The Inlet Start Face is Face " << InletStartFace << " .\n" << endl;

	//skip 6 lines

	for(int c27=1; c27<=6; c27++)
	{
		getline(boundaryFile,dummyStr);
	}

	//Skip 6 characters, >> 'Eats' the spaces.

	for(int c28=1; c28<=6; c28++)
	{
		boundaryFile >> dummyChar;
	}

	boundaryFile >> nOutletFaces;

	cout << "There are " << nOutletFaces << " Outlet Faces.\n" << endl;

	getline(boundaryFile,dummyStr);

	//Skip 9 characters

	for(int c29=1; c29<=9; c29++)
	{
		boundaryFile >> dummyChar;
	}

	boundaryFile >> OutletStartFace;

	cout << "The Outlet Start Face is Face " << OutletStartFace << " .\n" << endl;

	//skip 6 lines

	for(int c30=1; c30<=6; c30++)
	{
		getline(boundaryFile,dummyStr);
	}

	//Skip 6 characters, >> 'Eats' the spaces.

	for(int c31=1; c31<=6; c31++)
	{
		boundaryFile >> dummyChar;
	}

	boundaryFile >> nWallsFaces;

	cout << "There are " << nWallsFaces << " Walls Faces.\n" << endl;

	getline(boundaryFile,dummyStr);

	//Skip 9 characters

	for(int c32=1; c32<=9; c32++)
	{
		boundaryFile >> dummyChar;
	}

	boundaryFile >> WallsStartFace;

	cout << "The Walls Start Face is Face " << WallsStartFace << " .\n" << endl;

	boundaryFile.close();

	//Compute Face Types, 0 is Internal Face, 1 is Wall Face (Passive Boundary), 2 is Active Boundary Inlet, 3 is Active Boundary Outlet

	cout << "\nComputing Face Types...\n" << endl;

	int * FaceType;

	FaceType = new int [nFaces];

	//Set All FaceTypes to 0, Internal Face

	for(unsigned long int c33=0; c33<nFaces; c33++)
	{
		FaceType[c33]=0;
	}

	//Set Active Boundary Inlet FaceType, 2

	for(unsigned long int c34=InletStartFace; c34<(InletStartFace+nInletFaces); c34++)
	{
		FaceType[c34]=2;
	}

	//Set Active Boundary Outlet FaceType, 3

	for(unsigned long int c35=OutletStartFace; c35<(OutletStartFace+nOutletFaces); c35++)
	{
		FaceType[c35]=3;
	}

	//Set Walls Passive Boundary FaceType, 1

	for(unsigned long int c36=WallsStartFace; c36<(WallsStartFace+nWallsFaces); c36++)
	{
		FaceType[c36]=1;
	}

	cout << "Face Types Computed.\n" << endl;

	//Generate Cell Connectivity

	cout << "\nGenerating Cell Connectivity...\n" << endl;

	int * CellnConnectedCells;

	CellnConnectedCells = new int [nCells];

	//Set Cell nConnected Cells to 0

	for(unsigned long int c37=0; c37<nCells; c37++)
	{
		CellnConnectedCells[c37]=0;
	}

	//Scan number of Cell Connections for Each Cell

	for(unsigned long int c38=0; c38<nCells; c38++)
	{
		for(int c39=0; c39<6; c39++)
		{
			if(IsBoundaryFace[CellFaces[c38][c39]]==0)
			{
				CellnConnectedCells[c38]+=1;
			}
		}
	}

	cout << "Generated number of connected Cells for each Cell.\n" << endl;

	//Allocate Memory for Cell Connectivity List

	unsigned long int ** CellConnectedCells = new unsigned long int *[nCells];

	for(unsigned long int c40 = 0; c40 < nCells; c40++)
	{
		CellConnectedCells[c40] = new unsigned long int[CellnConnectedCells[c40]];
	}

	int * CellConnectionIndex;

	CellConnectionIndex = new int [nCells];

	for(unsigned long int c42=0; c42<nCells; c42++)
	{
		CellConnectionIndex[c42]=0;
	}

	//Build Cell Connectivity

	//unsigned long int FaceUnderConsideration=0;

	/*unsigned long int SearchCell=0;

	bool FaceFound=0;

	unsigned long int AnswerCell=0;*/

	for(unsigned long int c43=0; c43<nCells; c43++)
	{
		for(int c44=0; c44<6; c44++)
		{
			if(IsBoundaryFace[CellFaces[c43][c44]]==0)
			{
				for(unsigned long int c45=0; c45<nCells; c45++)
				{
					if(c45!=c43)
					{
						for(int c48=0; c48<6; c48++)
						{
							if(CellFaces[c45][c48]==CellFaces[c43][c44])
							{
								CellConnectedCells[c43][CellConnectionIndex[c43]]=c45;
								CellConnectionIndex[c43]+=1;
							}
						}
					}
				}
			}
		}
	}

	cout << "Cell Connectivity Built.\n" << endl;

	//Display Cell Connectivity

	for(unsigned long int c46=0; c46<nCells; c46++)
	{
		cout << "Cell " << c46 << " is connected to Cells";

		for(int c47=0; c47<CellnConnectedCells[c46]; c47++)
		{
			cout << " " << CellConnectedCells[c46][c47];
		}

		cout << endl;
	}

	//Diagnostics

	cout << "\n\nCellnConnected Cells for Cell 35999 is " << CellnConnectedCells[35999] << " .\n" << endl;

	cout <<"\nCell 35999 has faces " << CellFaces[35999][0] << " " << CellFaces[35999][1] << " " << CellFaces[35999][2] << " " << CellFaces[35999][3] << " " << CellFaces[35999][4] << " " << CellFaces[35999][5] << " .\n" << endl;

	cout << "IsBoundaryFace for these faces is " << IsBoundaryFace[CellFaces[35999][0]] << " " << IsBoundaryFace[CellFaces[35999][1]] << " " << IsBoundaryFace[CellFaces[35999][2]] << " " << IsBoundaryFace[CellFaces[35999][3]] << " " << IsBoundaryFace[CellFaces[35999][4]] << " " << IsBoundaryFace[CellFaces[35999][5]] << " .\n" << endl;

	//Okay

	//Simulation Initial and Boundary Conditions

	cout << "\nA practical example of the following simulation is a copper heat pipe inside a laptop or compact desktop, cooling a CPU or a GPU by conducting heat away from it and towards a heatsink near the exterior.\n" << endl;

	double T1initial=0;

	double T1inlet=0;

	double T1outlet=0;

	double T2initial=0;

	double T2inlet=0;

	double T2outlet=0;

	double Rho=0;

	double Cp=0;

	double k=0;

	double Alpha=0;

	cout << "\nEnter the initial uniform temperature of the heat pipe in degrees Celsius for Volume 1 (suggested ambient 20) : ";

	cin >> T1initial;

	cout << "\nEnter the initial uniform temperature of the heat pipe in degrees Celsius for Volume 2 (suggested ambient 30) : ";

	cin >> T2initial;

	cout << "\nEnter the Hot Inlet temperature of the heat pipe in degrees Celsius for Volume 1 (suggested hot CPU/GPU at 80): ";

	cin >> T1inlet;

	cout << "\nEnter the Hot Inlet temperature of the heat pipe in degrees Celsius for Volume 2 (suggested hot CPU/GPU at 100): ";

	cin >> T2inlet;

	cout << "\nEnter the Cold Outlet temperature of the heat pipe in degrees Celsius for Volume 1 (suggested heatsink at ambient 20): ";

	cin >> T1outlet;

	cout << "\nEnter the Cold Outlet temperature of the heat pipe in degrees Celsius for Volume 2 (suggested heatsink at ambient 30): ";

	cin >> T2outlet;

	cout << "\nEnter heat pipe density Rho in kilograms per m^3 (Suggested 8933 for Copper): ";

	cin >> Rho;

	cout << "\nEnter heat pipe heat capacity Cp in Joules per (Kg.K) (Suggested 385 for Copper): ";

	cin >> Cp;

	cout << "\nEnter heat pipe thermal conductivity k in Watts per (m.K) (Suggested 401 for Copper): ";

	cin >> k;

	cout << endl;

	Alpha=k/(Rho*Cp);

	//Simulation Parameters

	double TimeStep=0;

	double tEnd=0;

	unsigned long int SaveInterval=0;

	double DeltaX=0;

	cout << "\nEnter simulation TimeStep (suggested 0.001) in seconds: ";

	cin >> TimeStep;

	cout << "\nEnter duration of the simulation (suggested 120) in seconds: ";

	cin >> tEnd;

	cout << "\nEnter data save interval as a multiplication factor of TimeStep (suggested 1000): ";

	cin >> SaveInterval;

	cout << endl;

	//Begin Simulation

	cout <<"\n!!!!Warning!!!! Mesh Cells Are Assumed to be Cubic of Uniform Size.\n" << endl;

	cout << "Enter Uniform Cubic Mesh Cube Side in meters (suggested 0.001): ";

	cin >> DeltaX;

	cout << "\nStarting Simulation...\n" << endl;

	double * T1current; double * T1new; double * T2current; double * T2new;

	T1current = new double [nCells];

	T1new = new double [nCells];

	T2current = new double [nCells];

	T2new = new double [nCells];

	for(unsigned long int c49=0; c49<nCells; c49++)
	{
		T1current[c49]=T1initial;
	}

	for(unsigned long int c50=0; c50<nCells; c50++)
	{
		T1new[c50]=0;
	}

	for(unsigned long int c56=0; c56<nCells; c56++)
	{
		T2current[c56]=T2initial;
	}

	for(unsigned long int c57=0; c57<nCells; c57++)
	{
		T2new[c57]=0;
	}

	unsigned int TotalIterations=0;

	TotalIterations = tEnd / TimeStep;

	cout << "\nTotal number of iterations is " << TotalIterations << " iterations.\n" << endl;

	cout << "\nCreating 0 directory...\n" << endl;

	thread ZeroDir (Create0Dir);

	ZeroDir.join();

	cout << "\n0 directory created. Writing time 0 file for T1 and T2...\n" << endl;

	ofstream T10file;

	T10file.open ("./0/T1");

	T10file << "FoamFile\n";

	T10file << "{\n";

	T10file << "\tversion\t2.0;\n";

	T10file << "\tformat\tascii;\n";

	T10file << "\tclass\tvolScalarField;\n";

	T10file << "\tobject\tT1;\n";

	T10file << "}\n\n";

	T10file << "dimensions\t[0 0 0 1 0 0 0];\n\n";

	T10file << "internalField\tuniform " << T1initial << ";\n\n";

	T10file << "boundaryField\n";

	T10file << "{\n";

	T10file << "\tInlet\n";

	T10file << "\t{\n";

	T10file << "\t\ttype\t\tempty;\n";

	//T0file << "\t\tvalue\t\tuniform " << Tinlet << ";\n";

	T10file << "\t}\n\n";

	T10file << "\tOutlet\n";

	T10file << "\t{\n";

	T10file << "\t\ttype\t\tempty;\n";

	//T0file << "\t\tvalue\t\tuniform " << Toutlet << ";\n";

	T10file << "\t}\n\n";

	T10file << "\tWalls\n";

	T10file << "\t{\n";

	T10file << "\t\ttype\t\tempty;\n";

	T10file << "\t}\n";

	T10file << "}";

	T10file.close();

	ofstream T20file;

	T20file.open ("./0/T2");

	T20file << "FoamFile\n";

	T20file << "{\n";

	T20file << "\tversion\t2.0;\n";

	T20file << "\tformat\tascii;\n";

	T20file << "\tclass\tvolScalarField;\n";

	T20file << "\tobject\tT2;\n";

	T20file << "}\n\n";

	T20file << "dimensions\t[0 0 0 1 0 0 0];\n\n";

	T20file << "internalField\tuniform " << T2initial << ";\n\n";

	T20file << "boundaryField\n";

	T20file << "{\n";

	T20file << "\tInlet\n";

	T20file << "\t{\n";

	T20file << "\t\ttype\t\tempty;\n";

	//T0file << "\t\tvalue\t\tuniform " << Tinlet << ";\n";

	T20file << "\t}\n\n";

	T20file << "\tOutlet\n";

	T20file << "\t{\n";

	T20file << "\t\ttype\t\tempty;\n";

	//T0file << "\t\tvalue\t\tuniform " << Toutlet << ";\n";

	T20file << "\t}\n\n";

	T20file << "\tWalls\n";

	T20file << "\t{\n";

	T20file << "\t\ttype\t\tempty;\n";

	T20file << "\t}\n";

	T20file << "}";

	T20file.close();

	double Sum1=0.0; double Sum2=0.0;

	double DeltaXsquared=0.0;

	DeltaXsquared = DeltaX*DeltaX;

	unsigned long int SaveTime=0;

	string mkdirstr = "mkdir ";

	string command;

	const char* cmdstring;

	string PathStart = "./";

	string PathEnd1 = "/T1";
	string PathEnd2 = "/T2";

	string FilePathString1;
	string FilePathString2;

	const char* FilePath1;
	const char* FilePath2;

	//Main Simulation Loop

	for(unsigned int IterationCounter=1; IterationCounter<=TotalIterations; IterationCounter++)
	{
		cout << "\nStarting Iteration " << IterationCounter << " ...\n" << endl;

		for(unsigned long int c51=0; c51<nCells; c51++)
		{
			Sum1=0.0; Sum2=0.0;

			//Active Boundary Face Interactions

			for(int c52=0; c52<6; c52++)
			{	//Inlet
				if(FaceType[CellFaces[c51][c52]]==2)
				{
					Sum1+=(2*(T1inlet-T1current[c51]));
					Sum2+=(2*(T2inlet-T2current[c51]));
				}//Outlet
				else if(FaceType[CellFaces[c51][c52]]==3)
				{
					Sum1+=(2*(T1outlet-T1current[c51]));
					Sum2+=(2*(T2outlet-T2current[c51]));
				}
			}

			//Connected Cell Interactions

			for(int c53=0; c53<CellnConnectedCells[c51]; c53++)
			{
				Sum1+=(T1current[CellConnectedCells[c51][c53]]-T1current[c51]);
				Sum2+=(T2current[CellConnectedCells[c51][c53]]-T2current[c51]);
			}

			//Update Temperature

			T1new[c51] = T1current[c51] +((Alpha*TimeStep*Sum1)/(DeltaXsquared));
			T2new[c51] = T2current[c51] +((Alpha*TimeStep*Sum2)/(DeltaXsquared));

		}

		cout << "\nIteration " << IterationCounter << " Complete.\n" << endl;

		//Write Output At SaveInterval

		if(IterationCounter%SaveInterval==0)
		{
			SaveTime = (IterationCounter/SaveInterval);

			cout << "\nWriting Output for Time = " << SaveTime << " seconds.\n" << endl;

			command = mkdirstr + to_string(SaveTime); //std::to_string

			cmdstring = command.c_str();

			thread TimeDir (CreateDir,cmdstring);

			TimeDir.join();

			FilePathString1 = PathStart + to_string(SaveTime) + PathEnd1;
			FilePathString2 = PathStart + to_string(SaveTime) + PathEnd2;

			FilePath1 = FilePathString1.c_str();
			FilePath2 = FilePathString2.c_str();

			ofstream T1file;

			T1file.open(FilePath1);

			T1file << "FoamFile\n";

			T1file << "{\n";

			T1file << "\tversion\t2.0;\n";

			T1file << "\tformat\tascii;\n";

			T1file << "\tclass\tvolScalarField;\n";

			T1file << "\tlocation\t\"" << SaveTime << "\";\n";

			T1file << "\tobject\tT1;\n";

			T1file << "}\n\n";

			T1file << "dimensions\t[0 0 0 1 0 0 0];\n\n";

			T1file << "internalField\tnonuniform List<scalar>\n";

			T1file << nCells << "\n(\n";

			for(unsigned long int c55=0; c55<nCells; c55++)
			{
				T1file << T1new[c55] << endl;
			}

			T1file << ")\n";

			T1file << ";\n\n";

			T1file << "boundaryField\n";

			T1file << "{\n";

			T1file << "\tInlet\n";

			T1file << "\t{\n";

			T1file << "\t\ttype\t\tempty;\n";

			T1file << "\t}\n\n";

			T1file << "\tOutlet\n";

			T1file << "\t{\n";

			T1file << "\t\ttype\t\tempty;\n";

			T1file << "\t}\n\n";

			T1file << "\tWalls\n";

			T1file << "\t{\n";

			T1file << "\t\ttype\t\tempty;\n";

			T1file << "\t}" << endl;

			T1file << "}";

			T1file.close();

			ofstream T2file;

			T2file.open(FilePath2);

			T2file << "FoamFile\n";

			T2file << "{\n";

			T2file << "\tversion\t2.0;\n";

			T2file << "\tformat\tascii;\n";

			T2file << "\tclass\tvolScalarField;\n";

			T2file << "\tlocation\t\"" << SaveTime << "\";\n";

			T2file << "\tobject\tT2;\n";

			T2file << "}\n\n";

			T2file << "dimensions\t[0 0 0 1 0 0 0];\n\n";

			T2file << "internalField\tnonuniform List<scalar>\n";

			T2file << nCells << "\n(\n";

			for(unsigned long int c58=0; c58<nCells; c58++)
			{
				T2file << T2new[c58] << endl;
			}

			T2file << ")\n";

			T2file << ";\n\n";

			T2file << "boundaryField\n";

			T2file << "{\n";

			T2file << "\tInlet\n";

			T2file << "\t{\n";

			T2file << "\t\ttype\t\tempty;\n";

			T2file << "\t}\n\n";

			T2file << "\tOutlet\n";

			T2file << "\t{\n";

			T2file << "\t\ttype\t\tempty;\n";

			T2file << "\t}\n\n";

			T2file << "\tWalls\n";

			T2file << "\t{\n";

			T2file << "\t\ttype\t\tempty;\n";

			T2file << "\t}" << endl;

			T2file << "}";

			T2file.close();
		}

		//Prepare Variables for Next Iteration

		for(unsigned long int c54=0; c54<nCells; c54++)
		{
			T1current[c54]=T1new[c54];
			T1new[c54]=0.0;

			T2current[c54]=T2new[c54];
			T2new[c54]=0.0;
		}
	}

	//Free Dynamically Allocated Memory

	delete [] pointsX; delete [] pointsY; delete [] pointsZ;

	for(unsigned long int j =0; j < nFaces; j++)
	{
		delete [] FacePoints[j];
	}

	delete [] FacePoints;

	delete [] FaceOwnerCell;

	delete [] FaceNeighbourCell;

	delete [] IsBoundaryFace;

	for(unsigned long int c13 =0; c13 < nCells; c13++)
	{
		delete [] CellFaces[c13];
	}

	delete [] CellFaces;

	delete [] CellFaceIndex;

	delete [] FaceCenterX; delete [] FaceCenterY; delete [] FaceCenterZ;

	delete [] CellCenterX; delete [] CellCenterY; delete [] CellCenterZ;

	delete [] FaceType;

	delete [] CellnConnectedCells;

	for(unsigned long int c41=0; c41<nCells; c41++)
	{
		delete [] CellConnectedCells[c41];
	}

	delete [] CellConnectedCells;

	delete [] CellConnectionIndex;

	delete [] T1current; delete [] T1new; delete [] T2current; delete [] T2new;

	return 0;
}
