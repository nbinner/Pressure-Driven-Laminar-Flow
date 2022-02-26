/************************************************************************************************************
Mech 7225 Fluid Mechanics 2, Winter 2022
PRESSURE DRIVEN LAMINAR FLOW IN A RECTANGULAR CHANNEL

Purpose: The purpose of this project is to solve for the source term given a specific rectangular channel
         aspect ratio using a derived 2D Poisson's equation and using 2D using Finite-Difference (F-D) methods.
			Through finding the correct source term, the mean velocity may also be calculated and plotted on a 
			contour plot in Matlab. This velocity can then be used to find the vorticity vector field which is also
			plotted in Matlab.

Author(s):     Nathan Binner and Rehatbir Dhingra, previous project written by Nathan B., Muneer A. and Alex S.
Student ID(s): A01159743, A01087624

Declaration:
We, Rehatbir and Nathan declare that the following code was written by us.

Date Submitted: 2022-02-25
************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <float.h>

//------- GLOBAL CONSTANTS ----------------------------------------------------------------------------------

const char* ENDMESSAGE = "Thank you for running this program!";
const char* DASHES = "--------";
const char* TITLE_HEADER = "Pressure Driven Laminar Flow in a Rectangular Channel";

const int TABLE_COLUMN_WIDTH = 50;
const int TABLE_MARGIN_SIZE = 3;

const double PI = 3.141592653589793;
const int MAX_BUFF_SIZE = 1024;              // for reading lines from a file
const int MAX_ITER = 1000000;                // maximum iterations for F-D
const double MAX_RESIDUAL = 1.0e-8;          // solution accuracy
const double MAX_TRANSIENT_RESIDUAL = 0.01;  // transient solution accuracy

//--- Table Border Characters
const unsigned char HL = 196;  // horizontal border line
const unsigned char VL = 179;  // vertical border line
const unsigned char TL = 218;  // top left border symbol
const unsigned char TC = 194;  // top center border symbol
const unsigned char TR = 191;  // top right border symbol
const unsigned char CL = 195;  // left center border symbol
const unsigned char CC = 197;  // center center border symbol (cross)
const unsigned char CR = 180;  // right center border symbol
const unsigned char BL = 192;  // bottom left border symbol
const unsigned char BC = 193;  // bottom center border symbol
const unsigned char BR = 217;  // bottom right border symbol


//------- STRUCTURE DEFINITIONS -----------------------------------------------------------------------------
typedef struct PLATEPOINT
{
	double x, y;   // grid node physical position on the plate
	double W;      // grid node velocity
	double Vortx;  // vorticity in the x
	double Vorty;  // vorticity in the y
	double res;    // grid node residual for finite-difference solution
}
PLATEPOINT;


typedef struct SIMULATION_DATA     // holds data for each simulation
{
	double w, h, dx, dy, AR;        // plate width, plate height; x, y cellSizes and aspect ratio
	double S;                       // source term
	size_t I, J;                    // number of nodes in x and y directions
}
SIMULATION_DATA;


//------------------------- FUNCTION PROTOTYPES -------------------------------------------------------------
void waitForEnterKey();                                               // robust version of getchar (cleans input buffer)
bool flushInputBuffer2();                                             // return true if has garbage left in input buffer
void endProgram(const char*);                                         // terminates program with optional message
void printRepeatedChar(unsigned char, int);                           // prints a character repeatedly
void TitleBlock();                                                    // prints out title block to screen
SIMULATION_DATA GetSimulationData();                                  // reads a input file to obtain simulation data
int  printHorizontalBorder(char, char);                               // Prints the top or bottom border of the array display box
void GetNumericalSolution(PLATEPOINT**, SIMULATION_DATA);             // numerically calculates the solution of each case
void printSolution(PLATEPOINT**, const SIMULATION_DATA*);             // 2nd xmas present!  Prints contour plot data.
PLATEPOINT** initialize(SIMULATION_DATA, PLATEPOINT**);               // initializes the dynamic array
double DIntegral(PLATEPOINT**, SIMULATION_DATA);                      // sums the values for the double integral term
PLATEPOINT** IntergalacticZoomTheory(PLATEPOINT**, SIMULATION_DATA);  // calculates vorticity and stores it into the respective structure
void FreeMemory(PLATEPOINT**, SIMULATION_DATA*, size_t);              // frees the memory of the dynamically allocated arrays


//-----------------------------------------------------------------------------------------------------------
int main()
{
	PLATEPOINT** P = NULL;        // For 2D dynamically allocated array for a simulation
	SIMULATION_DATA SD;           // the array to hold simulation data for all cases in simulations.in

	SD = GetSimulationData();
	P = initialize(SD, P);
	GetNumericalSolution(P, SD);
	printSolution(P, &SD);
	FreeMemory(P, &SD, SD.I);
	waitForEnterKey();
	
	endProgram(NULL);
}

//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  This function reads from the input file to store the data uesd for the simulation
// ARGUMENTS:    *SD: the SIMULATION_DATA array, 
// RETURN VALUE: SD dynamic array, int NS
//               
//               
SIMULATION_DATA GetSimulationData()
{
	int ret;
	SIMULATION_DATA SD;

	TitleBlock();
	printf("\nEnter Aspect Ratio: ");
	ret = scanf_s("%lf", &SD.AR);
	if (ret != 1) EXIT_FAILURE;

	printf("\nEnter number of nodes in x direction: ");
	ret = scanf_s("%d", &SD.I);
	if (ret != 1) EXIT_FAILURE;

	// Calculate y nodes based on AR
	SD.J = (int)(SD.I / SD.AR);

	double AR = SD.AR;
	// Calculating w using AR
	double  w = SD.w;
	w = (AR + 1) / 2;
	SD.w = w;
	// Calculating h using AR
	double  h = SD.h;
	h = (1 + AR) / (2 * AR);
	SD.h = h;
	double dx = SD.w / (double)(SD.I - 1);
	SD.dx = dx;
	double dy = SD.h / (double)(SD.J - 1);
	SD.dy = dy;

	return SD;
}

//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  Dynamically allocates the memory for each array and intializes to 0 for each element 
//               and it intializes each x and y position of each node
// ARGUMENTS:    P:  the 2D PLATEPOINT array
//               SD: the simulation data for case B
// RETURN VALUE: PLATEPOINT P
PLATEPOINT** initialize(SIMULATION_DATA SD, PLATEPOINT** P)
{
	size_t n, i, j;
	size_t I = SD.I;
	size_t J = SD.J;
	P = (PLATEPOINT**)calloc(I, sizeof(PLATEPOINT*));
	if (P == NULL) exit(0);
	for (i = 0; i < I; i++)
	{
		P[i] = (PLATEPOINT*)calloc(J, sizeof(PLATEPOINT));
		if (P == NULL) exit(0);
	}
	if (I == 0 || J == 0) exit(0);
	for (i = 0; i < I; i++)
	{
		for (j = 0; j < J; j++)
		{	
			P[i][j].x = (double)i * SD.dx;
			P[i][j].y = (double)j * SD.dy;
		}
	}

	return P;
}

//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  Calculates the double integral term used to find velocity
//               
// ARGUMENTS:    P:  the 2D PLATEPOINT array
//               SD: the simulation data for case B
// RETURN VALUE: Sum Wsum of the double integral
double DIntegral(PLATEPOINT** P, SIMULATION_DATA SD)
{
	double Wsum = 0;
	size_t i, j;
	for (i = 1; i < SD.I - 1; i++)
	{
		for (j = 1; j < SD.J - 1; j++)
		{
			Wsum += P[i][j].W * SD.dx * SD.dy;
		}
	}

	return Wsum;
}

//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  Calculates the vorticity at each node of the channel
//               
// ARGUMENTS:    P:  the 2D PLATEPOINT array
//               SD: the simulation data for case B
// RETURN VALUE: PlATEPOINT P
PLATEPOINT** IntergalacticZoomTheory(PLATEPOINT** P, SIMULATION_DATA SD)
{
	size_t i = 0, j = 0;
	int I = SD.I;
	int J = SD.J;

	// Vorticity x component
	for (i = 1; i < I - 1; i++)
	{
		// Vorticity y component
		for (j = 1; j < J - 1; j++)
		{
			P[i][j].Vorty = -(P[i + 1][j].W - P[i - 1][j].W) / (2 * SD.dx);
			P[i][j].Vortx = ((P[i][j + 1].W - P[i][j - 1].W) / (2 * SD.dy));
		}
	}

	return P;
}

//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  Uses Finite-difference method to numerically solve for the temperature of each node 
//               Cycles through each node and finds the temperature based on the average of neighbouring nodes
// ARGUMENTS:    P:  the 2D PLATEPOINT array
//               SD: the simulation data for the selected case
// RETURN VALUE: none
void GetNumericalSolution(PLATEPOINT** P, SIMULATION_DATA SD)
{
	FILE* fConverge = NULL;
	errno_t err;
	char strConvergenceFile[MAX_BUFF_SIZE]; // convergence file string name

	double rmax = 0;   // defines and initializes rmax to zero
	int i = 0, j = 0;  // counters 

	double RMS = 0.0;  // variable holder for RMS value
	double dy = SD.dy; // defines dy from structure 
	double dx = SD.dx; // defines dx from structure 
	size_t I = SD.I;
	size_t J = SD.J;
	int iter = 0;      // iteration counter
	double Sg = -20.0; // hard coded source term guess
	double TOL = 1E-6;
	double Wsum = 0;

	double Wr, Wl, Wt, Wb;
	double lamda = (-2 / (dx * dx)) - (2 / (dy * dy));

	//-------------------Velocity calcs------------------------

	double Powdx = (1 / (dx * dx));
	double Powdy = (1 / (dy * dy));

	do
	{
		for (j = 1; j < J - 1; j++)
		{
			for (i = 1; i < I - 1; i++)
			{
				Wr = Powdx * P[i + 1][j].W;
				Wl = Powdy * P[i][j + 1].W;
				Wt = Powdx * P[i - 1][j].W;
				Wb = Powdy * P[i][j - 1].W;

				P[i][j].W = (SD.S - Wr - Wl - Wt - Wb) / lamda;
			}
		}

		RMS = 0.0;   // resets RMS to zero
		rmax = 0.0;  // resets rmax to zero

		for (j = 1; j < J - 1; j++)
		{
			for (i = 1; i < I - 1; i++)

			{
				Wr = Powdx * P[i + 1][j].W;
				Wl = Powdy * P[i][j + 1].W;
				Wt = Powdx * P[i - 1][j].W;
				Wb = Powdy * P[i][j - 1].W;
				P[i][j].res = fabs(P[i][j].W - ((SD.S - Wr - Wl - Wt - Wb) / lamda));
				if (P[i][j].res > rmax) rmax = P[i][j].res;
				RMS += (P[i][j].res * P[i][j].res);
			}
		}
		RMS = sqrt(RMS / (((double)I - 2) * ((double)J - 2))); // calculates RMS
		iter++; // iter increments 
		// do the loop while iter is less than or equal to MAX_ITER AND rmax is 
		//greater or eqal to MAX_RESIDUAL AND RMS greater or equal to MAX_RESIDUAL  

	} while (iter <= MAX_ITER && (rmax >= MAX_RESIDUAL && RMS >= MAX_RESIDUAL));

	//find the correct source term
	Wsum = DIntegral(P, SD);
	SD.S = (SD.S * SD.w * SD.h) / Wsum;
	printf("Correct source term = %lf", SD.S);
	Wsum = 0;

	do
	{

		for (j = 1; j < J - 1; j++)
		{
			for (i = 1; i < I - 1; i++)
			{

				Wr = Powdx * P[i + 1][j].W;
				Wt = Powdy * P[i][j + 1].W;
				Wl = Powdx * P[i - 1][j].W;
				Wb = Powdy * P[i][j - 1].W;

				P[i][j].W = (SD.S - Wr - Wl - Wt - Wb) / lamda;
			}
		}

		RMS = 0.0;  // resets RMS to zero
		rmax = 0.0; // resets rmax to zero

		for (j = 1; j < J - 1; j++)
		{
			for (i = 1; i < I - 1; i++)
			{
				Wr = Powdx * P[i + 1][j].W;
				Wl = Powdy * P[i][j + 1].W;
				Wt = Powdx * P[i - 1][j].W;
				Wb = Powdy * P[i][j - 1].W;
				P[i][j].res = fabs(P[i][j].W - ((SD.S - Wr - Wl - Wt - Wb) / lamda));
				if (P[i][j].res > rmax) rmax = P[i][j].res;
				RMS += (P[i][j].res * P[i][j].res);

			}
		}
		RMS = sqrt(RMS / (((double)I - 2) * ((double)J - 2))); // calculates RMS
		iter++; // iter increments 
		// do the loop while iter is less than or equal to MAX_ITER AND rmax is 
		//greater or eqal to MAX_RESIDUAL AND RMS greater or equal to MAX_RESIDUAL  

	} while (iter <= MAX_ITER && (rmax >= MAX_RESIDUAL && RMS >= MAX_RESIDUAL));
	//find the correct source term
	Wsum = DIntegral(P, SD);

	printf("\nIntegral of the velocity field using the correct source term = %lf", Wsum);
	printf("\nIntegral of velocity field = w*h = %lf", SD.w * SD.h);
	// prints to screen - the values of iter, rmax and RMS
	printf("\nNumber of iterations: %d", iter);
	printf("\nRmax = %.5le", rmax);
	printf("\nRMS = %.5le\n\n", RMS);

	// using the correct velocity, calculates vorticity
	P = IntergalacticZoomTheory(P, SD);
}

//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  Prints the solution to Matlab in order to display and graph the temperature distribution 
//               onto the steel plate
// ARGUMENTS:    P:  the 2D PLATEPOINT array
//               SD: the simulation data for case B
// RETURN VALUE: none

void printSolution(PLATEPOINT** P, const SIMULATION_DATA* pSD)
{
	size_t i, j;                                   // loop/temp variables
	FILE* fa = NULL, * ffd = NULL, * fres = NULL;  // for printing solutions to file for matlab
	char strFileNameAnalytical[MAX_BUFF_SIZE];     // buffer to hold analytical output file name
	char strFileNameFD[MAX_BUFF_SIZE];             // buffer to hold F-D output file name
	char strFileNameResidual[MAX_BUFF_SIZE];       // buffer to hold Residual output file name
	char strFileNameVorticity[MAX_BUFF_SIZE];      // buffer to hold Vorticity output file name
	errno_t err;                                   // check fopen

	// open finite different output file. Create name dynamically with sprintf_s
	sprintf_s(strFileNameFD, MAX_BUFF_SIZE, "%s Finite Difference.dat","Velocity - ");
	err = fopen_s(&ffd, strFileNameFD, "w");
	if (err != 0 || ffd == NULL)
	{
		printf("Cannot open \"%s\" for writing. Skipping printout...\n", strFileNameFD);
		return;
	}
	// open residual output file. Create name dynamically with sprintf_s
	sprintf_s(strFileNameResidual, MAX_BUFF_SIZE, "%s Residual.dat", "Velocity -");
	err = fopen_s(&fres, strFileNameResidual, "w");
	if (err != 0 || fres == NULL)
	{
		printf("Cannot open \"%s\" for writing. Skipping printout...\n", strFileNameResidual);
		return;
	}

	// open vorticity output file. Create name dynamically with sprintf_s
	sprintf_s(strFileNameVorticity, MAX_BUFF_SIZE, "%s.dat", "Vorticity");
	err = fopen_s(&fa, strFileNameVorticity, "w");
	if (err != 0 || fa == NULL)
	{
		printf("Cannot open \"%s\" for writing. Skipping printout...\n", strFileNameVorticity);
		return;
	}

	//--------------  print outputs to dat files for matlab -------------------------
	for (i = 0; i < pSD->I; i++)
	{
		for (j = 0; j < pSD->J; j++)
		{
			fprintf(fres, "%+12.5le,%+12.5le,%+12.5le\n", P[i][j].x, P[i][j].y, P[i][j].res);
			fprintf(ffd, "%+12.5le,%+12.5le,%+12.5le\n", P[i][j].x, P[i][j].y, P[i][j].W);
			fprintf(fa, "%+12.5le,%+12.5le,%+12.5le,%+12.5le\n", P[i][j].x, P[i][j].y, P[i][j].Vortx, P[i][j].Vorty);
		}
	}
	// close the files
	fclose(ffd);
	fclose(fres);
	fclose(fa);

	printf("Printed data to \"%s\"\n", strFileNameFD);
	printf("Printed data to \"%s\"\n", strFileNameResidual);
	printf("Printed data to \"%s\"\n", strFileNameVorticity);
}

//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  Frees the memory that was allocated to the dynamic arrays for the platepoint and SD struc 
// ARGUMENTS:    P:  the 2D PLATEPOINT array
//               SD: the simulation data for case B
//                I: total nodes in the i direction
// RETURN VALUE: none
void FreeMemory(PLATEPOINT** P, SIMULATION_DATA* SD, size_t I)
{
	size_t i;
	
	// Freeing Platepoint Array
	for (i = 0; i < SD->I; i++)
		free(P[i]);
	free(P);
}

//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  Prints a title block in the command window
// ARGUMENTS:    None
//              
// RETURN VALUE: None
void TitleBlock()
{
	printHorizontalBorder(TL, TR);
	printf("%c   %s   %c\n", VL, TITLE_HEADER, VL);
	printHorizontalBorder(BL, BR);
}

//--------------------------------------------------------------------------------------------
// DESCRIPTION:  Prints a horizontal border
// ARGUMENTS:    cLeft: the left border character, cRight: rght border character
// RETURN VALUE: int numChars
int printHorizontalBorder(char cLeft, char cRight)  // prints a character repeatedly
{
	int n;  // loop counter
	int numChars = 0; // number of characters printed by printfs

	numChars += printf("%c", cLeft);  // left corner
	for (n = 0; n < strlen(TITLE_HEADER) + 2 * TABLE_MARGIN_SIZE; n++) numChars += printf("%c", HL);  // mid section line
	numChars += printf("%c\n", cRight) - 1;  // right corner (-1 for newline)

	return numChars;
}

//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  This function flushes the input buffer to avoid scanf issues
//               ***** CALL THIS FUNCTION AFTER EVERY CALL TO SCANF!!! *****
// ARGUMENTS:    none
// RETURN VALUE: false if nothing or only '\n' is in the input buffer
//               true if extra keystrokes precede the '\n'.  Good for detecting left 
//               over garbage from scanf_s in the input buffer
bool flushInputBuffer2()
{
	unsigned char ch; // temp character variable
	bool bHasGarbage = false;

	// exit loop when all characters are flushed
	while ((ch = (unsigned char)getchar()) != '\n' && ch != EOF)
	{
		if (!bHasGarbage && !isspace(ch)) bHasGarbage = true;
	}
	return bHasGarbage;
}
//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  Waits for user to press enter.  flushes stdin if keystroke precede enter
// ARGUMENTS:    none
// RETURN VALUE: none
void waitForEnterKey()
{
	unsigned char ch;
	if ((ch = (unsigned char)getchar()) != EOF && ch != '\n') flushInputBuffer2();
}

//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  Display a message and terminate the program
// ARGUMENTS:    the message string
// RETURN VALUE: none
void endProgram(const char* message)
{
	printf("\n");
	if (message != NULL) printf(message);
	printf("\nPress ENTER to end this program...");
	waitForEnterKey();
	exit(0);
}

//-----------------------------------------------------------------------------------------------------------
// DESCRIPTION:  print a character repeatedly
// ARGUMENTS:    ch: the character, N: the number of repititions
// RETURN VALUE: none
void printRepeatedChar(unsigned char ch, int N)
{
	int n;
	for (n = 0; n < N; n++) printf("%c", ch);
}

