/*
 * main.c
 *
 * Code generation for function 'main'
 *
 */

 /*************************************************************************/
 /* This automatically generated example C main file shows how to call    */
 /* entry-point functions that MATLAB Coder generated. You must customize */
 /* this file for your application. Do not modify this file directly.     */
 /* Instead, make a copy of this file, modify it, and integrate it into   */
 /* your development environment.                                         */
 /*                                                                       */
 /* This file initializes entry-point function arguments to a default     */
 /* size and value before calling the entry-point functions. It does      */
 /* not store or use any values returned from the entry-point functions.  */
 /* If necessary, it does pre-allocate memory for returned values.        */
 /* You can use this file as a starting point for a main function that    */
 /* you can deploy in your application.                                   */
 /*                                                                       */
 /* After you copy the file, and before you deploy it, you must make the  */
 /* following changes:                                                    */
 /* * For variable-size function arguments, change the example sizes to   */
 /* the sizes that your application requires.                             */
 /* * Change the example values of function arguments to the values that  */
 /* your application requires.                                            */
 /* * If the entry-point functions return values, store these values or   */
 /* otherwise use them as required by your application.                   */
 /*                                                                       */
 /*************************************************************************/

 /* Include files */
#include "main.h"
#include "gyrotron.h"
#include "gyrotron_emxAPI.h"
#include "gyrotron_terminate.h"
#include "rt_nonfinite.h"

static void main_gyrotron(double Ne, int Nzz, double Lz, double Tend, double Delta,
	double I0, double R0, double Rb, double g, double ukv, double dz, double dt, double tol, double INTT, double INTZ, double SPLINE, double DK, double a0)
{
	emxArray_real_T* OUTFre;
	emxArray_real_T* OUTFim;
	emxArray_real_T* OUTJre;
	emxArray_real_T* OUTJim;
	emxArray_real_T* OUTZAxis;
	emxArray_real_T* OUTTAxis;
	emxArray_real_T* Eff;
	emxArray_real_T* Omega;
	double jout;
	emxInitArray_real_T(&OUTFre, 2);
	emxInitArray_real_T(&OUTFim, 2);
	emxInitArray_real_T(&OUTJre, 2);
	emxInitArray_real_T(&OUTJim, 2);
	emxInitArray_real_T(&OUTZAxis, 1);
	emxInitArray_real_T(&OUTTAxis, 1);
	emxInitArray_real_T(&Eff, 1);
	emxInitArray_real_T(&Omega, 1);

	FILE* fre, *fim, * ire, *iim, * feff, * fomega, * fileID;
	errno_t err;
	bool convert;
	double Ic, gamma, betta, betta2, betta_z, betta_z2, betta_perp2, gamma0,
		c, e, m, nu, w_op, ZetaEx, TauEnd;

	int Nz, Nt, OUTNz, OUTNt, idx0, idx1;
	static int iv2[2];

	convert = true;
	Ic = 0;
	if (Nzz == 0) {
		Nz = (int)(Lz / dz) + 1;
		Nt = (int)(Tend / dt) + 1;
		Ic = I0;
		convert = false;
	}

	gamma = 1.0 + ukv / 511.0;
	betta = sqrt(1.0 - 1.0 / (gamma * gamma));
	betta2 = 1.0 - 1.0 / (gamma * gamma);
	betta_z = betta / sqrt(g * g + 1.0);
	betta_z2 = betta2 / (g * g + 1.0);
	betta_perp2 = betta2 - betta_z2;
	gamma0 = sqrt(1 - betta_perp2 - betta_z);

	c = 29979245800; // [cm / s]
	e = 4.803e-10; // [וה.ֳׁׁ]
	m = 9.1093837015e-28; // [g]

	nu = 73.952055635763557;
	w_op = c * nu / R0;

	if (convert == true) {
		Nz = Nzz;
		ZetaEx = betta_perp2 / 2.0 / betta_z * w_op * Lz / c;
		TauEnd = betta_perp2 * betta_perp2 * w_op * Tend / 8 / betta_z2;
		dz = ZetaEx / (Nz-1);
		dt = DK * dz * dz;
	}
	else {
		ZetaEx = Lz;
		TauEnd = Tend;
	}

	Nt = (int)(Tend / dt) + 1;

	if (INTT < 1) {
		printf("Too small Tend\n");
		exit(-1);
	}
	if (INTZ < 1) {
		printf("Too small Zend\n");
		exit(-1);
	}

	if (INTT > 1 && INTZ > 1) {
		OUTNt = (int)(Nt / INTT) + 1;
		OUTNz = (int)(Nz / INTZ) + 1;
	}
	else if (INTT == 1 && INTZ > 1) {
		OUTNt = Nt;
		OUTNz = (int)(Nz / INTZ) + 1;
	}
	else if (INTT > 1 && INTZ == 1) {
		OUTNt = (int)(Nt / INTT) + 1;
		OUTNz = Nz;
	}
	else {
		OUTNt = Nt;
		OUTNz = Nz;
	}

	iv2[0] = OUTNz;
	iv2[1] = OUTNt;

	OUTFre = emxCreateND_real_T(2, iv2);
	OUTFim = emxCreateND_real_T(2, iv2);
	OUTJre = emxCreateND_real_T(2, iv2);
	OUTJim = emxCreateND_real_T(2, iv2);
	OUTZAxis = emxCreateND_real_T(1, &OUTNz);
	OUTTAxis = emxCreateND_real_T(1, &OUTNt);
	Eff = emxCreateND_real_T(1, &OUTNt);
	Omega = emxCreateND_real_T(1, &OUTNt);

	idx0 = 1;
	idx1 = 1;
	//printf("Number %f", OUTBvsT->data[idx0 + OUTBvsT->size[0] * idx1]);

	err = fopen_s(&fileID, "input_c.txt", "w");
	fprintf(fileID, "Nz = % i\n", Nz);
	fprintf(fileID, "Nt = % i\n", Nt);
	fprintf(fileID, "Ne = % i\n", (int)Ne);
	fprintf(fileID, "ZetaEx = % f\n", ZetaEx);
	fprintf(fileID, "TauEnd = % f\n", TauEnd);
	fprintf(fileID, "Delta = % f\n", Delta);
	fprintf(fileID, "I0 = % f\n", I0);
	fprintf(fileID, "R0 = % f\n", R0);
	fprintf(fileID, "Rb = % f\n", Rb);
	fprintf(fileID, "g = % f\n", g);
	fprintf(fileID, "ukv = % f\n", ukv);
	fprintf(fileID, "dz = % f\n", dz);
	fprintf(fileID, "dt = % f\n", dt);
	fprintf(fileID, "tol = % g\n", tol);
	fclose(fileID);

	/* Call the entry-point 'orotron'. */
	gyrotron(Ne, Nzz, Lz, Tend, Delta, I0, R0, Rb, g, ukv, dz, dt, tol, INTT, INTZ, SPLINE, DK, a0,
		OUTFre, OUTFim, OUTJre, OUTJim, OUTZAxis, OUTTAxis, Eff, Omega, &jout);

	err = fopen_s(&fre, "fre.dat", "w");
	if (err == 0)
	{
		printf("The file 'fre.dat' was opened\n");
	}
	else
	{
		printf("The file 'fre.dat' was not opened\n");
	}

	err = fopen_s(&fim, "fim.dat", "w");
	if (err == 0)
	{
		printf("The file 'fim.dat' was opened\n");
	}
	else
	{
		printf("The file 'fim.dat' was not opened\n");
	}

	err = fopen_s(&ire, "ire.dat", "w");
	if (err == 0)
	{
		printf("The file 'ire.dat' was opened\n");
	}
	else
	{
		printf("The file 'ire.dat' was not opened\n");
	}

	err = fopen_s(&iim, "iim.dat", "w");
	if (err == 0)
	{
		printf("The file 'iim.dat' was opened\n");
	}
	else
	{
		printf("The file 'iim.dat' was not opened\n");
	}

	err = fopen_s(&feff, "e.dat", "w");
	if (err == 0)
	{
		printf("The file 'e.dat' was opened\n");
	}
	else
	{
		printf("The file 'e.dat' was not opened\n");
	}

	err = fopen_s(&fomega, "w.dat", "w");
	if (err == 0)
	{
		printf("The file 'w.dat' was opened\n");
	}
	else
	{
		printf("The file 'w.dat' was not opened\n");
	}

	for (idx0 = 0; idx0 < OUTNz; idx0++) {
		fprintf(fre, "%e\t", OUTZAxis->data[idx0]);
		for (int idx1 = 0; idx1 < OUTNt; idx1++) {
			fprintf(fre, "%e\t", OUTFre->data[idx0 + OUTFre->size[0] * idx1]);
		}
		fprintf(fre, "\n");
	}

	for (idx0 = 0; idx0 < OUTNz; idx0++) {
		fprintf(fim, "%e\t", OUTZAxis->data[idx0]);
		for (int idx1 = 0; idx1 < OUTNt; idx1++) {
			fprintf(fim, "%e\t", OUTFim->data[idx0 + OUTFim->size[0] * idx1]);
		}
		fprintf(fim, "\n");
	}

	for (idx0 = 0; idx0 < OUTNz; idx0++) {
		fprintf(ire, "%e\t", OUTZAxis->data[idx0]);
		for (int idx1 = 0; idx1 < OUTNt; idx1++) {
			fprintf(ire, "%e\t", OUTJre->data[idx0 + OUTJre->size[0] * idx1]);
		}
		fprintf(ire, "\n");
	}

	for (idx0 = 0; idx0 < OUTNz; idx0++) {
		fprintf(iim, "%e\t", OUTZAxis->data[idx0]);
		for (int idx1 = 0; idx1 < OUTNt; idx1++) {
			fprintf(iim, "%e\t", OUTJim->data[idx0 + OUTJim->size[0] * idx1]);
		}
		fprintf(iim, "\n");
	}

	for (int idx1 = 0; idx1 < Nt; idx1++) {
		fprintf(feff, "%e\t%e\n", idx1*dt, Eff->data[idx1]);
	}

	for (int idx1 = 0; idx1 < Nt; idx1++) {
		fprintf(fomega, "%e\t%e\n", idx1*dt, Omega->data[idx1]);
	}


	if (ire)
	{
		err = fclose(ire);
		if (err == 0)
		{
			printf("The file 'ire.dat' was closed\n");
		}
		else
		{
			printf("The file 'ire.dat' was not closed\n");
		}
	}

	if (iim)
	{
		err = fclose(iim);
		if (err == 0)
		{
			printf("The file 'iim.dat' was closed\n");
		}
		else
		{
			printf("The file 'iim.dat' was not closed\n");
		}
	}

	if (fre)
	{
		err = fclose(fre);
		if (err == 0)
		{
			printf("The file 'fre.dat' was closed\n");
		}
		else
		{
			printf("The file 'fre.dat' was not closed\n");
		}
	}

	if (fim)
	{
		err = fclose(fim);
		if (err == 0)
		{
			printf("The file 'fim.dat' was closed\n");
		}
		else
		{
			printf("The file 'fim.dat' was not closed\n");
		}
	}

	if (feff)
	{
		err = fclose(feff);
		if (err == 0)
		{
			printf("The file 'e.dat' was closed\n");
		}
		else
		{
			printf("The file 'e.dat' was not closed\n");
		}
	}

	if (fomega)
	{
		err = fclose(fomega);
		if (err == 0)
		{
			printf("The file 'w.dat' was closed\n");
		}
		else
		{
			printf("The file 'w.dat' was not closed\n");
		}
	}

	emxDestroyArray_real_T(Omega);
	emxDestroyArray_real_T(Eff);
	emxDestroyArray_real_T(OUTTAxis);
	emxDestroyArray_real_T(OUTZAxis);
	emxDestroyArray_real_T(OUTJim);
	emxDestroyArray_real_T(OUTJre);
	emxDestroyArray_real_T(OUTFim);
	emxDestroyArray_real_T(OUTFre);
}

int main(int argc, const char* const argv[])
{
	int Nzz, INTT, INTZ;
	double Ne, Lz, Tend, Delta, I0, R0, Rb, g, ukv, dz, dt, tol, SPLINE, DK, a0;

	if (argc != 19) {
		printf("Expected 18 arguments: Ne, Nz, Lz, Tend, Delta, I0, R0, Rb, g, ukv, dz, dt, tol, INTT, INTZ, SPLINE, DK, a0\n");
		exit(-1);
	}

	Ne = atof(argv[1]);
	Nzz = atoi(argv[2]);
	Lz = atof(argv[3]);
	Tend = atof(argv[4]);
	Delta = atof(argv[5]);
	I0 = atof(argv[6]);
	R0 = atof(argv[7]);
	Rb = atof(argv[8]);
	g = atof(argv[9]);
	ukv = atof(argv[10]);
	dz = atof(argv[11]);
	dt = atof(argv[12]);
	tol = atof(argv[13]);
	INTT = atoi(argv[14]);
	INTZ = atoi(argv[15]);
	SPLINE = atof(argv[16]);
	DK = atof(argv[17]);
	a0 = atof(argv[18]);

	/* The initialize function is being called automatically from your entry-point function. So, a call to initialize is not included here. */
	/* Invoke the entry-point functions.
	   You can call entry-point functions multiple times. */
	main_gyrotron(Ne, Nzz, Lz, Tend, Delta, I0, R0, Rb, g, ukv, dz, dt, tol, INTT, INTZ, SPLINE, DK, a0);

	/* Terminate the application.
	   You do not need to do this more than one time. */
	gyrotron_terminate();
	return 0;
}

/* End of code generation (main.c) */
