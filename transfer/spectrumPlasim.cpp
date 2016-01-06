#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <vector>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_int.h>
#include <gsl/gsl_matrix.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "arlnsmat.h"
#include "arlsnsym.h"
#include "atio.hpp"


typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SpMatCSC;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMatCSR;


// Declarations

void writeSpectrum(FILE *, FILE *, double *, double *, double *, int, int);


int main(int argc, char * argv[])
{
  if (argc <= 1) {
    std::cerr << "argc <= 1!" << std::endl << "Usage: " << argv[0] << " nev" << std::endl;
    return -1;
  }
  double S = 1360;
  char restartState[] = "warm";
  int lastYear = 9999;
  int firstYear = 101;
  char processing[] = "_yearly";
  const char *indexChoice[1];
  indexChoice[0] = "nhemisic";
  size_t dim = 1;
  // const char *indexChoice[2];
  // indexChoice[0] = "nhemisic";
  // indexChoice[1] = "eqmst";
  // size_t dim = 2;

  // Grid definition
  int nx = 50;
  int nSTD = 5;

  // Lags
  int nLags = 1;
  gsl_vector *tauDimRng = gsl_vector_alloc(nLags);
  gsl_vector_set(tauDimRng, 0, 1.);

  // Eigen problem parameters
  int nev = atoi(argv[1]);
  bool solveAdjoint = true;
  const std::string& which = "LM"; // Look for eigenvalues of Largest Magnitude
  int ncv = 0;
  double tol =0.;
  int maxit = 0;
  double *resid = NULL;
  bool AutoShift = true;
  

  // Definitions

  // Lags
  double tauDim;

  // Names and Files
  char resDir[128], dstDir[128], obsName[128], caseName[128], mkdirCmd[128];
  char postfix[128], gridPostfix[128], cpyBuffer[128];
  char transitionMatrixFileName[128], EigValFileName[128], EigVecFileName[128];
  FILE *transitionMatrixFile, *EigValFile, *EigVecFile;
  int N = 1;

  sprintf(resDir, "%s_%d", restartState, (int) round(S*10));
  strcpy(dstDir, resDir);
  sprintf(caseName, "%s_%05d_%05d", resDir, firstYear, lastYear);
  sprintf(mkdirCmd, "mkdir %s %s/spectrum %s/spectrum/eigval %s/spectrum/eigvec 2> /dev/null",
	  dstDir, dstDir, dstDir, dstDir);
  system(mkdirCmd);

  strcpy(obsName, caseName);
  sprintf(gridPostfix, "N");
  for (size_t d = 0; d < dim; d++) {
    N *= nx;
    sprintf(obsName, "%s_%s", obsName, indexChoice[d]);
    strcpy(cpyBuffer, gridPostfix);
    if (d > 0)
      sprintf(gridPostfix, "%sx%d", cpyBuffer, nx);
    else
      sprintf(gridPostfix, "%s%d", cpyBuffer, nx);
  }
  strcpy(cpyBuffer, gridPostfix);
  sprintf(gridPostfix, "%s_%s_%s_%dstd", processing, obsName, cpyBuffer, nSTD);

  
  // Transition Matrix
  SpMatCSR *PCSR;
  SpMatCSC *PCSC;
  ARluNonSymMatrix<double, double> *P, *PT;

  // Eigen problem
  double *EigValReal = new double [nev+1];
  double *EigValImag = new double [nev+1];
  double *EigVec = new double [(nev+2)*N];
  ARluNonSymStdEig<double> EigProb;
  int nconv;

  
  // Get transition matrices for different lags
  for (int lag = 0; lag < nLags; lag++){
    tauDim = gsl_vector_get(tauDimRng, lag);

    // Open source and destination files
    sprintf(postfix, "%s_tau%02d", gridPostfix, (int) tauDim);
    sprintf(transitionMatrixFileName, "%s/transitionMatrix/transitionMatrix%s.csr",
	    dstDir, postfix);
    if ((transitionMatrixFile = fopen(transitionMatrixFileName, "r")) == NULL){
      fprintf(stderr, "Can't open %s for reading!\n", transitionMatrixFileName);
      return -1;
    }


    // Read transition matrix written in CSR as CSC to take the transpose
    std::cout << endl << "Reading transpose transition matrix for lag " << tauDim
	      << " in " << transitionMatrixFileName << std::endl;
    PT = Compressed2AR(transitionMatrixFile);

    // Declare eigen problem: real non-symetric (see examples/areig.h)
    std::cout << "Solving eigen problem for the first " << nev
	      << " eigenvalues" << std::endl;
    EigProb = ARluNonSymStdEig<double>(nev, *PT, which, ncv, tol, maxit,
				       resid, AutoShift);

    // Open destination files and write spectrum
    sprintf(EigValFileName, "%s/spectrum/eigval/eigval_nev%d%s.txt", dstDir, nev, postfix);
    sprintf(EigVecFileName, "%s/spectrum/eigvec/eigvec_nev%d%s.txt", dstDir, nev, postfix);
    if ((EigValFile = fopen(EigValFileName, "w")) == NULL){
      fprintf(stderr, "Can't open %s for writing!\n", EigValFileName);
      return -1;
    }
    if ((EigVecFile = fopen(EigVecFileName, "w")) == NULL){
      fprintf(stderr, "Can't open %s for writing!\n", EigVecFileName);
      return -1;
    }
  
    // Find eigenvalues and left eigenvectors
    EigProb.EigenValVectors(EigVec, EigValReal, EigValImag);
    nconv = EigProb.ConvergedEigenvalues();
    std::cout << "Found " << nconv << "/" << nev << " eigenvalues." << std::endl;

    // Write results
    std::cout << "Write eigenvalues to " << EigValFileName << std::endl;
    std::cout << "and eigenvectors to " << EigVecFileName << std::endl;
    writeSpectrum(EigValFile, EigVecFile, EigValReal, EigValImag, EigVec, nev, N);
    fclose(EigValFile);
    fclose(EigVecFile);

    delete PT;

    // Sove the adjoint problem on the transpose (of the transpose)
    if (solveAdjoint){
      // Read transition matrix written in CSR and convert to CSC
      std::cout << std::endl << "Reading transition matrix for lag " << tauDim
		<< " in " << transitionMatrixFileName << std::endl;
      if (fseek(transitionMatrixFile, 0, SEEK_SET) != 0){
	fprintf(stderr, "Can't rewind %s calling fseek!\n", transitionMatrixFileName);
	return -1;
      }

      PCSR = Compressed2Eigen(transitionMatrixFile);
      fclose(transitionMatrixFile);
      std::cout << "Convert from CSR to CSC" << std::endl;
      PCSC = CSR2CSC(PCSR);
      std::cout << "Convert from Eigen to Arpack" << std::endl;
      P = Eigen2AR(PCSC);

      // Declare eigen problem: real non-symetric (see examples/areig.h)
      std::cout << "Solving eigen problem for the first " << nev
		<< " eigenvalues" << std::endl;
      EigProb = ARluNonSymStdEig<double>(nev, *P, which, ncv, tol, maxit, resid, AutoShift);

      // Open destination files and write spectrum
      sprintf(EigValFileName, "%s/spectrum/eigval/eigvalAdjoint_nev%d%s.txt",
	      dstDir, nev, postfix);
      sprintf(EigVecFileName, "%s/spectrum/eigvec/eigvecAdjoint_nev%d%s.txt",
	      dstDir, nev, postfix);
      if ((EigValFile = fopen(EigValFileName, "w")) == NULL){
	fprintf(stderr, "Can't open %s for writing!\n", EigValFileName);
	return -1;
      }
      if ((EigVecFile = fopen(EigVecFileName, "w")) == NULL){
	fprintf(stderr, "Can't open %s for writing!\n", EigVecFileName);
	return -1;
      }
  
      // Find eigenvalues and left eigenvectors
      EigProb.EigenValVectors(EigVec, EigValReal, EigValImag);
      nconv = EigProb.ConvergedEigenvalues();
      std::cout << "Found " << nconv << "/" << nev << " adjoint eigenvalues." << std::endl;

      // Write results
      std::cout << "Write adjoint eigenvalues to " << EigValFileName << std::endl;
      std::cout << "and adjoint eigenvectors to " << EigVecFileName << std::endl;
      writeSpectrum(EigValFile, EigVecFile, EigValReal, EigValImag, EigVec, nev, N);
      fclose(EigValFile);
      fclose(EigVecFile);
      
      delete PCSR;
      delete PCSC;
      delete P;
    }
  }
  
  // Clean-up
  delete[] EigValReal;
  delete[] EigValImag;
  delete[] EigVec;
  
  return 0;
}


// Definitions

// Write complex eigenvalues and eigenvectors
void writeSpectrum(FILE *fEigVal, FILE *fEigVec, double *EigValReal, double *EigValImag,
		   double *EigVec, int nev, int N)
{
  int vecCount = 0;
  int ev =0;
  // Write real and imaginary parts of each eigenvalue on each line
  // Write on each pair of line the real part of an eigenvector then its imaginary part
  while (ev < nev) {
    // Always write the eigenvalue
    fprintf(fEigVal, "%lf %lf\n", EigValReal[ev], EigValImag[ev]);
    // Always write the real part of the eigenvector ev
    for (int i = 0; i < N; i++){
      fprintf(fEigVec, "%lf ", EigVec[vecCount*N+i]);
    }
    fprintf(fEigVec, "\n");
    vecCount++;
    
    // Write its imaginary part or the zero vector
    if (EigValImag[ev] != 0.){
      for (int i = 0; i < N; i++)
	fprintf(fEigVec, "%lf ", EigVec[vecCount*N+i]);
      vecCount++;
      // Skip the conjugate
      ev += 2;
    }
    else{
      for (int i = 0; i < N; i++)
	fprintf(fEigVec, "%lf ", 0.);
      ev += 1;
    }
    fprintf(fEigVec, "\n");
  }

  return;
}
	
