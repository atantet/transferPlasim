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
#include "transferOperator.hpp"
#include "atio.hpp"

typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SpMatCSC;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMatCSR;


int main(int argc, char * argv[])
{
  char srcDir[] = "../runPlasim/postprocessor/indices/";
  double S = 1360;
  char restartState[] = "warm";
  int lastYear = 9999;
  int firstYear = 101;
  size_t daysPerYear = 360;
  char processing[] = "_yearly";
  const char *indexChoice[1];
  indexChoice[0] = "nhemisic";
  size_t dim = 1;
  //const char *indexChoice[2];
  //indexChoice[0] = "nhemisic";
  //indexChoice[1] = "eqmst";
  //size_t dim = 2;

  size_t spinupYears = 200;
  size_t sampling = daysPerYear;
  size_t spinup = spinupYears * daysPerYear;
  double dt = 1.;

  // Grid definition
  int nx = 50;
  int nSTD = 5;

  // Lags
  int nLags = 1;
  gsl_vector *tauDimRng = gsl_vector_alloc(nLags);
  gsl_vector_set(tauDimRng, 0, 1.);

  // Names and Files
  char resDir[128], dstDir[128], indicesPath[128], obsName[128], caseName[128], mkdirCmd[128];
  char postfix[128], gridPostfix[128], cpyBuffer[128];
  char indexFileName[128], statesFileName[128], gridFileName[128],
    gridMemFileName[128], transitionMatrixFileName[128];
  FILE *indexFile, *statesFile, *gridFile, *gridMemFile, *transitionMatrixFile;
  
  // Vectors and matrices
  size_t tauStep, tauDim;
  double sample;
  gsl_vector_int *gridMem;
  SpMatCSR *P;
  size_t nt = (size_t) ((lastYear - firstYear + 1) * daysPerYear);
  size_t ntSample = (nt - spinup) / sampling;
  gsl_vector *index = gsl_vector_alloc(nt);
  gsl_matrix *states = gsl_matrix_alloc(ntSample, dim);
  gsl_vector *statesMean = gsl_vector_calloc(dim);
  gsl_vector *statesSTD = gsl_vector_calloc(dim);

  // Grid related
  int N = 1;
  double delta;
  gsl_matrix *gridLims;
  gsl_vector_int *nBox;
  std::vector<gsl_vector *> *gridBounds;
  gsl_vector *xmin = gsl_vector_alloc(dim);
  gsl_vector *xmax = gsl_vector_alloc(dim);

  
  // Define names and make directories
  sprintf(resDir, "%s_%d", restartState, (int) round(S*10));
  strcpy(dstDir, resDir);
  sprintf(indicesPath, "%s/%s/", srcDir, resDir);
  sprintf(caseName, "%s_%05d_%05d", resDir, firstYear, lastYear);
  sprintf(mkdirCmd, "mkdir %s %s/transitionMatrix %s/grid %s/obs 2> /dev/null",
	  dstDir, dstDir, dstDir, dstDir);
  system(mkdirCmd);
	 

  // Read indices
  strcpy(obsName, caseName);
  sprintf(gridPostfix, "N");
  for (size_t d = 0; d < dim; d++) {
    sprintf(indexFileName, "%s/%s_%s.txt", indicesPath, indexChoice[d], caseName);
    sprintf(obsName, "%s_%s", obsName, indexChoice[d]);
    std::cout << "Reading index file " << indexFileName << std::endl;
    if ((indexFile = fopen(indexFileName, "r")) == NULL){
      fprintf(stderr, "Can't open %s for reading!\n", indexFileName);
      return -1;
    }
    
    // Read trajectory
    gsl_vector_fscanf(indexFile, index);
    fclose(indexFile);

    // Take the average for each year k starting after the spinup
    // and calculate mean and std used to define the grid
    for (size_t k = 0; k < ntSample; k++){
      // Average over each day of the year
      sample = 0;
      for (size_t s = 0; s < sampling; s++)
	sample += gsl_vector_get(index, (k + spinupYears)*sampling + s);
      sample /= sampling;

      // Set state 
      gsl_matrix_set(states, k, d, sample);

      // Calculate mean and STD
      gsl_vector_set(statesMean, d,
		     gsl_vector_get(statesMean, d) + gsl_matrix_get(states, k, d));
      gsl_vector_set(statesSTD, d,
		     gsl_vector_get(statesSTD, d) + pow(gsl_matrix_get(states, k, d), 2));
    }
    gsl_vector_set(statesMean, d, gsl_vector_get(statesMean, d) / ntSample);
    gsl_vector_set(statesSTD, d, sqrt(gsl_vector_get(statesSTD, d) / ntSample
				      - pow(gsl_vector_get(statesMean, d), 2)));
    gsl_vector_set(xmin, d, gsl_vector_get(statesMean, d) - nSTD * gsl_vector_get(statesSTD, d));
    gsl_vector_set(xmax, d, gsl_vector_get(statesMean, d) + nSTD * gsl_vector_get(statesSTD, d));
    strcpy(cpyBuffer, gridPostfix);
    if (d > 0)
      sprintf(gridPostfix, "%sx%d", cpyBuffer, nx);
    else
      sprintf(gridPostfix, "%s%d", cpyBuffer, nx);
  }
  strcpy(cpyBuffer, gridPostfix);
  sprintf(gridPostfix, "%s_%s_%s_%dstd", processing, obsName, cpyBuffer, nSTD);


  // Save states
  sprintf(statesFileName, "%s/obs/obs%s.txt", dstDir, gridPostfix);
  std::cout << "Writing states in " << statesFileName << std::endl;
  if ((statesFile = fopen(statesFileName, "w")) == NULL){
    fprintf(stderr, "Can't open %s for writing!\n", statesFileName);
    return -1;
  }
  gsl_matrix_fprintf(statesFile, states, "%lf");
  fclose(statesFile);

  // Open grid file
  sprintf(gridFileName, "%s/grid/grid%s.txt", dstDir, gridPostfix);
  if ((gridFile = fopen(gridFileName, "w")) == NULL){
    fprintf(stderr, "Can't open %s for writing!\n", gridFileName);
    return -1;
  }
  // Define grid
  gridLims = gsl_matrix_alloc(dim, 2);
  nBox = gsl_vector_int_alloc(dim);
  gridBounds = new std::vector<gsl_vector *>(dim);
  std::cout << "Domain grid (min, max, n):" << std::endl;
  for (size_t d = 0; d < dim; d++) {
    gsl_vector_int_set(nBox, d, nx);
    N *= gsl_vector_int_get(nBox, d);
    gsl_matrix_set(gridLims, d, 0, gsl_vector_get(xmin, d));
    gsl_matrix_set(gridLims, d, 1, gsl_vector_get(xmax, d));
    std::cout << "dim " << d+1 << ": (" << gsl_matrix_get(gridLims, d, 0) << ", "
	      << gsl_matrix_get(gridLims, d, 1) << ", " <<gsl_vector_int_get(nBox, d) << ")"
	      << std::endl;
    // Alloc one dimensional box boundaries vector
    (*gridBounds)[d] = gsl_vector_alloc(gsl_vector_int_get(nBox, d) + 1);
    // Get spatial step
    delta = (gsl_matrix_get(gridLims, d, 1) - gsl_matrix_get(gridLims, d, 0))
      / gsl_vector_int_get(nBox, d);
    gsl_vector_set((*gridBounds)[d], 0, gsl_matrix_get(gridLims, d, 0));
    fprintf(gridFile, "%lf ", gsl_vector_get((*gridBounds)[d], 0));
    for (int i = 1; i < gsl_vector_int_get(nBox, d) + 1; i++){
      gsl_vector_set((*gridBounds)[d], i, gsl_vector_get((*gridBounds)[d], i-1) + delta);
      fprintf(gridFile, "%lf ", gsl_vector_get((*gridBounds)[d], i));
    }
    fprintf(gridFile, "\n");
  }
  fclose(gridFile);

  
  // Open grid membership file, get and write
  sprintf(gridMemFileName, "%s/transitionMatrix/gridMem%s.txt", dstDir, gridPostfix);
  if ((gridMemFile = fopen(gridMemFileName, "w")) == NULL){
    fprintf(stderr, "Can't open %s for writing!\n", gridMemFileName);
    return -1;
  }

  std::cout << "Getting grid membership vector..." << std::endl;
  gridMem = getGridMembership(states, gridBounds);
  gsl_vector_int_fprintf(gridMemFile, gridMem, "%d");
  fclose(gridMemFile);

  
  // Get transition matrices for ifferent lags
  for (int lag = 0; lag < nLags; lag++){
    tauDim = gsl_vector_get(tauDimRng, lag);
    tauStep = tauDim * daysPerYear / (dt * sampling);

    // Open transition matrix
    sprintf(postfix, "%s_tau%02d", gridPostfix, (int) tauDim);
    sprintf(transitionMatrixFileName, "%s/transitionMatrix/transitionMatrix%s.csr",
	    dstDir, postfix);
    if ((transitionMatrixFile = fopen(transitionMatrixFileName, "w")) == NULL){
      fprintf(stderr, "Can't open %s for writing!\n", transitionMatrixFileName);
      return -1;
    }

    // Get transition matrix as CSR
    std::cout << "Getting transition matrix" << std::endl;
    P = getTransitionMatrix(gridMem, N, tauStep);
    
    // Write transition matrix as CSR
    std::cout << "Writing transition matrix to " << transitionMatrixFileName << std::endl;
    Eigen2Compressed(transitionMatrixFile, P);
    fclose(transitionMatrixFile);
  }

  
  for (size_t d = 0; d < dim; d++)
    gsl_vector_free((*gridBounds)[d]);
  delete gridBounds;
  gsl_vector_int_free(nBox);
  gsl_matrix_free(gridLims);
  gsl_vector_free(index);
  gsl_matrix_free(states);
  gsl_vector_int_free(gridMem);
  gsl_vector_free(tauDimRng);
  gsl_vector_free(statesMean);
  gsl_vector_free(statesSTD);
  gsl_vector_free(xmin);
  gsl_vector_free(xmax);
  delete P;
		
  return 0;
}
