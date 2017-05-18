/*########################################################################
Copyright (c) 2014-2016, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.

Created by Geoffrey Oxberry (oxberry1@llnl.gov, goxberry@gmail.com),
Lluis-Miquel Munguia Conejero (lluis.munguia@gatech.edu), and Deepak
Rajan (rajan3@llnl.gov). LLNL-CODE-699387. All rights reserved.

This file is part of PIPS-SBB. For details, see
https://github.com/llnl/PIPS-SBB.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (as
published by the Free Software Foundation) version 2.1, February 1999.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
conditions of the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
########################################################################*/
#include "BBSMPSHeuristicScenDecom.hpp"

using namespace std;

// Concatenate first-stage vector of bounds or objective function coefficients
// with the analogous second-stage quantities for a given scenario
double* concatOneScenVec(int scen, stochasticInput& input,
			 std::vector<double> (stochasticInput::*first)(),
			 std::vector<double> (stochasticInput::*second)(int))
{
  // allocate size first
  const std::vector<double>& vec1 = (input.*first)();
  const std::vector<double>& vec2 = (input.*second)(scen);
  double* out = new double[vec1.size() + vec2.size()];

  unsigned i = 0;
  for (i=0; i < vec1.size(); i++) {
    out[i] = vec1.at(i);
  }
  for (unsigned k = 0; k < vec2.size(); k++) {
    out[i++] = vec2.at(k);
  }

  return out;
}

// Concatenate the matrices of constraint coefficients from the first-stage
// constraints plus the second-stage constraints from a given scenario
CoinPackedMatrix* concatOneScenMat(int scen, stochasticInput& input)
{
  const CoinPackedMatrix& Amat = input.getFirstStageConstraints();
  const CoinPackedMatrix& Tmat = input.getLinkingConstraints(scen);
  const CoinPackedMatrix& Wmat = input.getSecondStageConstraints(scen);

  int nvar1 = input.nFirstStageVars();
  int nvar2 = input.nSecondStageVars(scen);
  int ncons1 = input.nFirstStageCons();
  int ncons2 = input.nSecondStageCons(scen);

  int totalVar = nvar1 + nvar2;
  int totalCons = ncons1 + ncons2;

  CoinBigIndex totalNnz = Amat.getNumElements()
    + Tmat.getNumElements() + Wmat.getNumElements();

  // CoinPackedMatrix owns pointers passed to it
  CoinBigIndex *starts = new CoinBigIndex[totalVar+1];
  double *elts = new double[totalNnz];
  int *rowIdx = new int[totalNnz];

  CoinBigIndex nnz = 0, start, end;
  // put first-stage variables first, as is customary
  int const *Aidx = Amat.getIndices();
  double const *Aelts = Amat.getElements();
  int const *Tidx = Tmat.getIndices();
  double const *Telts = Tmat.getElements();

  // Build columns of the constraint matrix corresponding to first
  // stage variables
  for (int c = 0; c < nvar1; c++) {
    // Build A matrix part (first-stage constraints)
    starts[c] = nnz;
    start = Amat.getVectorFirst(c);
    end = Amat.getVectorLast(c);
    for (CoinBigIndex j = start; j < end; j++) {
      elts[nnz] = Aelts[j];
      rowIdx[nnz++] = Aidx[j];
    }

    // Build T matrix part (linking constraints for second stage)
    start = Tmat.getVectorFirst(c);
    end = Tmat.getVectorLast(c);
    for (CoinBigIndex j = start; j < end; j++) {
      elts[nnz] = Telts[j];
      rowIdx[nnz++] = Tidx[j]+ncons1;
    }

  }

  // Copy W blocks
  int rowOffset = ncons1;
  int colOffset = nvar1;
  int const *Widx = Wmat.getIndices();
  double const *Welts = Wmat.getElements();

  for (int c = 0; c < nvar2; c++) {
    starts[colOffset++] = nnz;
    start = Wmat.getVectorFirst(c);
    end = Wmat.getVectorLast(c);
    for (CoinBigIndex j = start; j < end; j++) {
      elts[nnz] = Welts[j];
      rowIdx[nnz++] = Widx[j]+rowOffset;
    }
  }
  starts[totalVar] = nnz;
  assert(nnz == totalNnz);

  CoinPackedMatrix *constr = new CoinPackedMatrix();
  int *lens = NULL;
  constr->assignMatrix(true,totalCons,totalVar,totalNnz,elts,rowIdx,starts,lens);
  constr->verifyMtx(); // debug matrix

  return constr;
}

void checkFSBinary(stochasticInput& input, bool &fsInteger, bool &fsBinary)
{

  int nvar1 = input.nFirstStageVars();
  const vector<double>& lb = input.getFirstStageColLB();
  const vector<double>& ub = input.getFirstStageColUB();

  // initialize
  fsInteger = true;
  fsBinary = true;

  for (unsigned col=0; col < nvar1; col++) {
    // not integer or binary
    if (!input.isFirstStageColInteger(col)) {
      fsInteger = false;
      fsBinary = false;
      return;
    }
    // not binary
    if (!isBinary(lb[col],ub[col],1e-6)) {
      fsBinary = false;
    }
  }

  return;

}

void makeOneScenOsiSolverInterface(int scen, stochasticInput& input,
				   OsiSolverInterface *si) {

  // OsiClpSolverInterface takes ownership of these
  CoinPackedMatrix *constraints = concatOneScenMat(scen, input);

  double *collb = concatOneScenVec(scen, input, &stochasticInput::getFirstStageColLB,
                                   &stochasticInput::getSecondStageColLB);
  double *colub = concatOneScenVec(scen, input, &stochasticInput::getFirstStageColUB,
                                   &stochasticInput::getSecondStageColUB);
  double *obj = concatOneScenVec(scen, input, &stochasticInput::getFirstStageObj,
                                 &stochasticInput::getSecondStageObj);
  double *rowlb = concatOneScenVec(scen, input, &stochasticInput::getFirstStageRowLB,
                                   &stochasticInput::getSecondStageRowLB);
  double *rowub = concatOneScenVec(scen, input, &stochasticInput::getFirstStageRowUB,
                                   &stochasticInput::getSecondStageRowUB);

  int nvar1 = input.nFirstStageVars();
  int nvar2 = input.nSecondStageVars(scen);

  si->assignProblem(constraints, collb, colub, obj, rowlb, rowub);

  // Get integer information as well
  std::vector<int> integerIdx;
  for (int i = 0; i < nvar1; i++) {
    if(input.isFirstStageColInteger(i)) { integerIdx.push_back(i); }
  }
  for (int i = 0; i < nvar2; i++) {
    if(input.isSecondStageColInteger(scen, i)) { integerIdx.push_back(i + nvar1); }
  }
  // Could also use integerIdx.data() in C++11 and later
  si->setInteger(&integerIdx[0], integerIdx.size());
}

BBSMPSHeuristicScenDecom::~BBSMPSHeuristicScenDecom(){

  BBSMPS_ALG_LOG_SEV(debug) << "Cleaning up scen_wrap vector" << endl;
  // One less element in scen_wrap than size.
  for (unsigned i=0; i<scen_wrap.size(); i++) {
    delete scen_wrap[i];
  }

}

bool BBSMPSHeuristicScenDecom::runHeuristic(BBSMPSNode* node,
					    denseBAVector &nodeSolution){

  bool   success = false;
  double *solVec;

  // get current time
  double startTimeStamp = MPI_Wtime();

  // get context and input
  BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
  SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
  // get original dimensions
  const BADimensionsSlacks &dimsOriginalSlacks
    = BBSMPSSolver::instance()->getOriginalBADimensionsSlacks();

  // problem size
  int nvar1 = input.nFirstStageVars();
  int ncon1 = input.nFirstStageCons();
  int nscen = input.nScenarios();

  // MPI information
  int nproc = ctx.nprocs();
  int mype=BBSMPSSolver::instance()->getMype();

  BBSMPS_ALG_LOG_SEV(debug) << "Problem size: " << nvar1 << " " << ncon1 << " " << nscen;

  BBSMPS_ALG_LOG_SEV(debug) << "MPI info: " << nproc << " " << mype;
  if (0 == mype) {
    BBSMPS_ALG_LOG_SEV(debug) << "Performing the ScenDecom Heuristic.";
  }
  timesCalled++;

  // get the local scenarios
  vector<int> const &localScen = ctx.localScenarios();
  int nlocalscen = localScen.size();
  BBSMPS_ALG_LOG_SEV(debug) << "Number of local scenarios: " << nlocalscen;

  // get solution pool and current UB
  int originalSolutionPoolSize=BBSMPSSolver::instance()->getSolPoolSize();
  double objUB=COIN_DBL_MAX;
  if (BBSMPSSolver::instance()->getSolPoolSize()>0)
    objUB=BBSMPSSolver::instance()->getSoln(0).getObjValue();

  // Steps for the heuristic

  // get bounds for current node, so that in later nodes,
  // we use current branching information
  denseBAVector lb(BBSMPSSolver::instance()->getOriginalLB());
  denseBAVector ub(BBSMPSSolver::instance()->getOriginalUB());
  node->getAllBranchingInformation(lb,ub);

  // Build single scenario object for objects.

  // To reuse scen_wrap, we need to
  // 1. Not clear scen_wrap. Done.
  // 2. Initialize scen_wrap only once, using additional indicator for first time. Done.
  // 3. Do not update node variable bounds in initialization. Done.
  // 4. Refactor bound update code. Not to be refactored, since only used in one place.
  // 5. Update scen_wrap using new bound update code every time. Done.

  // TODO: Can I use existing solution pool to provide incumbents to CBC?

  // Useful to access CbcModel directly
  CbcModel *cbcModel;

  // First time entering the heuristic
  if (firstTime) {

    // Resize vector to the size it will have

    BBSMPS_ALG_LOG_SEV(debug) << "Size of scen_wrap: " << nlocalscen-1;

    scen_wrap.resize(nlocalscen-1);

    BBSMPS_ALG_LOG_SEV(debug) << "Resized scen_wrap vector";

    checkFSBinary(input, fsInteger, fsBinary);

    BBSMPS_ALG_LOG_SEV(debug) << "Problem is first stage integer/binary: "
				<< fsInteger << " " << fsBinary;

    // resizes localFsSolns
    localFsSolns.resize((nlocalscen-1)*nvar1);
    // assigns aveFsSoln and sets to 0
    aveFsSoln.assign(nvar1,0);

    for(unsigned localscen = 1; localscen < nlocalscen; localscen++) {

      int scen = localScen[localscen];

      // uses the default constructor
      OsiCbcSolverInterface *wrap = new OsiCbcSolverInterface(NULL,NULL);

      // should I be using calls from rootSolver instead. Yes, if there is presolve.
      // To be revisited.
      // Similarly OSi problem wrapper could be given the presolved problem.
      // many dimensions could break if there is presolve.
      // will have to stop using input.

      makeOneScenOsiSolverInterface(scen, input, wrap);

      // Add node and time limits
      wrap->setMaximumNodes(nodeLim);
      wrap->setMaximumSeconds(timeLim);

      // Add to queue.
      scen_wrap[localscen-1]=wrap;

    }

  }

  BBSMPS_ALG_LOG_SEV(debug) << "scen_wrap is now initialized";

  // can set firstTime = false. scen_wrap is now initialized.
  firstTime = false;

  for(unsigned localscen = 1; localscen < nlocalscen; localscen++) {

    int scen = localScen[localscen];

    int nvar2 = input.nSecondStageVars(scen);
    int ncon2 = input.nSecondStageCons(scen);

    // reset bounds to original bounds
    /*
    if (localscen==1) {

      BBSMPS_ALG_LOG_SEV(debug) << "Problem size: " << nvar2 << " " << ncon2 << " " << localscen;

      // Print original lb
      const double *testlb = scen_wrap[localscen-1]->getColLower();
      cout << "Printing previous lower bound" << endl;
      copy(testlb, testlb+nvar1+nvar2, ostream_iterator<double>(cout, " "));
      cout << endl;
      //TODO: Something is messed up here, but we reset bounds, so for now it is okay.

      cout << "Printing node lb, including slacks" << endl;
      lb.getFirstStageVec().print();
      lb.getSecondStageVec(scen).print();

    }
    */

    // reset to original bounds
    // first stage
    for (unsigned col = 0; col < nvar1; col++) {
      scen_wrap[localscen-1]->setColLower(col,input.getFirstStageColLB()[col]);
      scen_wrap[localscen-1]->setColUpper(col,input.getFirstStageColUB()[col]);
    }
    // second stage
    for (unsigned col = 0; col < nvar2; col++) {
      scen_wrap[localscen-1]->setColLower(nvar1+col,input.getSecondStageColLB(scen)[col]);
      scen_wrap[localscen-1]->setColUpper(nvar1+col,input.getSecondStageColUB(scen)[col]);
    }

    /*
    if (localscen==1) {

      // Print original lb
      const double *testlb = scen_wrap[localscen-1]->getColLower();
      cout << "Printing previous lower bound" << endl;
      copy(testlb, testlb+nvar1+nvar2, ostream_iterator<double>(cout, " "));
      cout << endl;
    }
    */

    // modify bounds based on node bounds

    if (updateBnds) {

      // first stage
      for (unsigned col = 0; col < nvar1; col++) {
	scen_wrap[localscen-1]->setColLower(col,lb.getFirstStageVec()[col]);
	scen_wrap[localscen-1]->setColUpper(col,ub.getFirstStageVec()[col]);
      }

      // second stage
      for (unsigned col = 0; col < nvar2; col++) {
	scen_wrap[localscen-1]->setColLower(nvar1+col,lb.getSecondStageVec(scen)[col]);
	scen_wrap[localscen-1]->setColUpper(nvar1+col,ub.getSecondStageVec(scen)[col]);
      }

      // Print modified lb
      /*
      if (localscen>0) {
	const double *modlb = scen_wrap[localscen-1]->getColLower();
	cout << "lb after modifying" << endl;
	copy(modlb, modlb+nvar1+nvar2, ostream_iterator<double>(cout, " "));
	cout << endl;
      }
      */

    }

    // Important note: Need to offset scen_wrap by -1 since localscen starts at 1

  }

  // if first stage is binary, I can add no-good cut
  if (fsBinary) {

    // uniqueness check already done
    // read from fsSet object

    set<vector<bool> >::iterator it;
    for (it = fsSet.begin(); it != fsSet.end(); ++it) {

      // add to Cbc objects

      // build row object
      vector<double> elts(nvar1);
      double rowlb = 1;
      double rowub = nvar1;

      for (unsigned col = 0; col < nvar1; col++) {

	// if the col element is true
	if ((*it)[col]) {
	  rowlb -= 1;
	  rowub -= 1;
	  elts[col] = -1.0;
	}
	else {
	  elts[col] = 1.0;
	}

      }

      CoinPackedVector vec(nvar1, &elts[0]);

      // for each local scenario
      for(unsigned localscen = 1; localscen < nlocalscen; localscen++) {
	scen_wrap[localscen-1]->addRow(vec, rowlb, rowub);
      }

      /*
      cout << "First stage solution cut added" << endl;
      copy((*it).begin(), (*it).end(), ostream_iterator<double>(cout, " "));
      cout << endl;
      */

    }
  }

  // clearing the fs cut set
  fsSet.clear();

  BBSMPS_ALG_LOG_SEV(debug) << "Done building single scenario objects";

  // size vector of first stage solutions for allgather
  // one per rank, each of length nvar1
  // first stage solutions for each rank
  vector<double> fsSolutions(nproc*nvar1);
  // indicator for whether solution was found
  vector<int> foundSol(nproc);

  // To modify objective
  // Replace c_i with c_i (1 + \rho (\hat x - \bar x)). Done.
  // We need to store \hat x for each local scenario: using localFsSolns
  // TODO: What if it was infeasible?
  // \bar x = old \bar x + (new \hat x - old \hat x)/nscen. Done.
  // We keep track of \bar x in each rank, use AveFsSoln
  // For each scenario solved del = (new \hat x - old \hat x)/nscen is calculated
  // All ranks reduce del and get sum across ranks. Done.
  // All ranks add reduced del to \bar x. Done.

  // For one scenario in rank, run CBC with settings

  BBSMPS_ALG_LOG_SEV(debug) << "Pre increment Local/actual scenario: "
					   << local_scen_num << " "
					   << localScen[local_scen_num];

  // increment if no sameScen
  if (!sameScen) {
    local_scen_num++; // increment
    if (local_scen_num == nlocalscen) {

      local_scen_num = 1; // handle overflow

      allSolvedOnce = true; // all have been solved once

    }
  }

  // indicator for finding feasible solution
  int found_sol = 1;

  // Write mps; just for debugging
  //  scen_wrap[local_scen_num-1]->writeMps("temp","mps");

  BBSMPS_ALG_LOG_SEV(debug) << "Post increment Local/actual scenario: "
					   << local_scen_num << " "
					   << localScen[local_scen_num];

  // get pointer to underlying model, easier to query
  cbcModel = scen_wrap[local_scen_num-1]->getModelPtr();
  cbcModel->setLogLevel(0);

  // Code for modifying objective funtion.
  // replace c_i with c_i (1 + \rho (\hat x - \bar x))
  // test code path with \bar x = \hat x
  //  memcpy(&aveFsSoln[0], &localFsSolns[(local_scen_num-1)*nvar1], nvar1*sizeof(double));
  double rho = 1.0;
  vector<double> diffobj(nvar1);
  // get objective coefficients
  const double *oldobjVec = scen_wrap[local_scen_num-1]->getObjCoefficients();
  // Print objective vector
  /*
  cout << "Old Objective vector" << endl;
  copy(oldobjVec, oldobjVec+nvar1, ostream_iterator<double>(cout, " "));
  cout << endl;
  */
  // TODO: Ideally, should wait till all ranks are solved once. But
  // it is a heuristic, so for now changing more often
  // Need to check I am being correct, since they update at different times

  // updating delta in objective function if all have been solved once
  for (unsigned col=0; col<nvar1; col++) {
    diffobj[col] = (allSolvedOnce) ? oldobjVec[col] * rho *
      (localFsSolns[(local_scen_num-1)*nvar1+col] - aveFsSoln[col]) : 0;
  }

  // Just to be safe, resetting model solve
  cbcModel->resetModel();

  // Modify objective with penalty
  for (unsigned col=0; col<nvar1; col++) {
    scen_wrap[local_scen_num-1]->setObjCoeff(col, oldobjVec[col]+diffobj[col]);
  }

  // Print objective vector
  /*
  const double *printobjVec = scen_wrap[local_scen_num-1]->getObjCoefficients();
  cout << "Modified Objective vector" << endl;
  copy(printobjVec, printobjVec+nvar1, ostream_iterator<double>(cout, " "));
  cout << endl;
  */
  /*
  cout << "All local First stage solution before solving" << endl;
  copy(localFsSolns.begin(), localFsSolns.end(), ostream_iterator<double>(cout, " "));
  cout << endl;
  */

  // Solving directly with Cbc
  cbcModel->branchAndBound();

  // Modify objective back
  const double *newobjVec = scen_wrap[local_scen_num-1]->getObjCoefficients();
  for (unsigned col=0; col<nvar1; col++) {
    scen_wrap[local_scen_num-1]->setObjCoeff(col, newobjVec[col]-diffobj[col]);
  }

  // Print objective vector

  const double *resetobjVec = scen_wrap[local_scen_num-1]->getObjCoefficients();
  cout << "Reset back to original Objective vector" << endl;
  copy(resetobjVec, resetobjVec+nvar1, ostream_iterator<double>(cout, " "));
  cout << endl;


  // clear up diffobj
  diffobj.clear();

  // initialize to NULL
  solVec = NULL;
  // get best solution
  solVec = cbcModel->bestSolution();
  vector<double> fsSolVec;
  // abort if all no solution, mark as no solution
  if (solVec == NULL) {
    // assign dummy entries
    fsSolVec.assign(nvar1,0);
    BBSMPS_ALG_LOG_SEV(debug) << "Found no feasible solution";
    found_sol = 0;
  }
  else {
    // assign to local FS soln vector
    fsSolVec.assign(solVec, solVec+nvar1);
  }

  // Need to calculate new \hat x - old \hat x, else 0
  vector<double> diffFsSol(nvar1);
  for (unsigned col=0; col<nvar1; col++) {
    diffFsSol[col] = (allSolvedOnce) ? (fsSolVec[col] - localFsSolns[(local_scen_num-1)*nvar1+col])/ (nscen + 0.0) : fsSolVec[col] / (nscen + 0.0);
  }

  /*
  cout << "Diff First stage solution" << endl;
  copy(diffFsSol.begin(), diffFsSol.end(), ostream_iterator<double>(cout, " "));
  cout << endl;
  */

  // assigns to localFsSolns
  memcpy(&localFsSolns[(local_scen_num-1)*nvar1], &fsSolVec[0], nvar1*sizeof(double));

  // Print first stage solution vector
  /*
  cout << "All local First stage solution" << endl;
  copy(localFsSolns.begin(), localFsSolns.end(), ostream_iterator<double>(cout, " "));
  cout << endl;
  */

  // MPI call to share which ranks found first-stage solutions
  MPI_Allgather(&found_sol, 1, MPI_INT,
		&foundSol[0], 1, MPI_INT,
		ctx.comm());

  // MPI allreduce call to check if anybody failed.
  MPI_Allreduce(MPI_IN_PLACE, &found_sol, 1,
		MPI_INT, MPI_SUM,
		ctx.comm());

  // TODO: Maybe only solve some of them by guessing which one is best

  // Abort if all ranks found no solution
  BBSMPS_ALG_LOG_SEV(debug) << "Number of ranks that found a solution " << found_sol;
  /*
  cout << "Vector of feasible solution indicators" << endl;
  copy(foundSol.begin(), foundSol.end(), ostream_iterator<double>(cout, " "));
  cout << endl;
  */

  // Continue as long as someone has a solution
  // Be able to continue with fewer/more solutions than ranks.
  // Done with foundSol indicators
  // only quit if all processes found no solution
  if (found_sol == 0) {

    BBSMPS_ALG_LOG_SEV(debug) << "All ranks found no feasible solution. Aborting";
    cumulativeTime+=(MPI_Wtime()-startTimeStamp);
    return success;

  }
  // Joint first stage solution vector pre all gather
  /*
  cout << "First stage solution before allgather" << endl;
  copy(fsSolutions.begin(), fsSolutions.end(), ostream_iterator<double>(cout, " "));
  cout << endl;
  */

  // MPI call to share all first stage solutions.
  MPI_Allgather(&fsSolVec[0], nvar1, MPI_DOUBLE,
		&fsSolutions[0], nvar1, MPI_DOUBLE,
		ctx.comm());
  // Some of these may be all zeros.

  // MPI call to share all delta first stage solutions.
  // Now each rank has the sum of all the deltas
  MPI_Allreduce(MPI_IN_PLACE, &diffFsSol[0], nvar1, MPI_DOUBLE, MPI_SUM,
		ctx.comm());

  /*
  cout << "Diff First stage solution after reduction" << endl;
  copy(diffFsSol.begin(), diffFsSol.end(), ostream_iterator<double>(cout, " "));
  cout << endl;
  */

  // TODO: Change to <algorithm> calls
  // Update average Fs soln
  for (unsigned col=0; col<nvar1; col++) {
    aveFsSoln[col] += diffFsSol[col];
    // TODO: Other places where there might be a numerical issue.
    if (isZero(aveFsSoln[col],1e-6)) aveFsSoln[col] = 0.0;
  }

  /*
  cout << "Average First stage solution after update" << endl;
  copy(aveFsSoln.begin(), aveFsSoln.end(), ostream_iterator<double>(cout, " "));
  cout << endl;
  */

  diffFsSol.clear();

  // Solutions may not be unique. Mark uniqueness first
  // check uniqueness by putting in set for FS integer?
  // FS binary is done.
  if (fsBinary) {


    for(unsigned soln = 0; soln < nproc; soln++) {

      // ignore if no solution found.
      if (!foundSol[soln]) continue;

      // construct bool vector
      vector<bool> fsbool(nvar1);
      for (unsigned col = 0; col < nvar1; col++) {
	fsbool[col] = (isOne(fsSolutions[soln*nvar1+col],1e-6)) ? true : false;
      }

      // check if solution soln is unique
      if (fsSet.find(fsbool) == fsSet.end()) {

	// insert
	fsSet.insert(fsbool);

	/*
	cout << "First stage solution " << soln << " added to be cut later" << endl;
	copy(fsbool.begin(), fsbool.end(), ostream_iterator<double>(cout, " "));
	cout << endl;
	*/

      }

      else {

	// mark as duplicate
	BBSMPS_ALG_LOG_SEV(debug) << "Solution is duplicate: " << soln;

	foundSol[soln] = 0;
	// decrement
	found_sol --;

      }
    }

  }

  // Joint first stage solution vector post all gather
  /*
  cout << "First stage solution after allgather" << endl;
  copy(fsSolutions.begin(), fsSolutions.end(), ostream_iterator<double>(cout, " "));
  cout << endl;
  */

  // Now, all ranks have all the solutions.
  // And all ranks know which solution should be used.

  // size vector for second stage solutions, one vector for each scenario
  ssSolutions.resize(nlocalscen);

  // size vector for obj values for each solution
  // second stage will be reduced across ranks
  // also sets to 0.
  fsObjvalues.assign(nproc,0);
  ssObjvalues.assign(nproc,0);
  // Need to collect across all scenarios to evaluate a solution.

  // TODO: different nvar2 for each scenario?

  BBSMPS_ALG_LOG_SEV(debug) << "Ready to solve second stage problems";

  // Doesn't matter which is the outer loop as long as we wait till all solves are
  // done for a solution before computing its objective
  // Useful to do scenario in outer loop beccause of Osi object

  // for now, we wait till all solves are done
  // Could be more intelligent doing one solution at a time - may not need to
  // sync certain solutions if not of good quality, but not worth the extra wait
  // time for checking. Better off solving everything before checking.

  // each ssSolutions holds everything only for its rank
  // ssSolutions [ one for each local scenario ]
  // for each local scenario, the second stage variables for each solution
  // Important note: need to offset by -1 since localscen starts at 1.

  // Don't need to sync ssSolution, but I need to store it, just in case
  // this solution is the best.

  // resets value
  found_sol = 1;

  // for each local scenario
  for(unsigned localscen = 1; localscen < nlocalscen; localscen++) {

    // abort if some rank found no feasible solution
    if (!found_sol) break;

    int scen = localScen[localscen];
    int nvar2 = input.nSecondStageVars(scen);
    int ncon2 = input.nSecondStageCons(scen);

    BBSMPS_ALG_LOG_SEV(debug) << "Local/actual scenario: " << localscen << " "
					     << scen;

    // nproc solutions in fsSolutions.
    // Need to solve second stage problem for each of them.
    // for each scenario, and each solution vector is nvar2 long
    ssSolutions[localscen-1].resize(nproc*nvar2);

    // define vector for objective coefficients
    const double *objVec;

    // Important Note: Probabilities are not scaled into first stage problem

    // For each solution, solve second stage problem for current scenario,
    for(unsigned soln = 0; soln < nproc; soln++) {

      // ignore solving second stage if not feasible
      if (!foundSol[soln]) continue;

      BBSMPS_ALG_LOG_SEV(debug) << "Considering solution: " << soln;

      /*
      cout << "First stage solution being fixed" << endl;
      copy(fsSolutions.begin(), fsSolutions.end(), ostream_iterator<double>(cout, " "));
      cout << endl;
      */

      // get pointer to underlying model, easier to query
      cbcModel = scen_wrap[localscen-1]->getModelPtr();
      cbcModel->setLogLevel(0);

      // Just to be safe, resetting model solve
      cbcModel->resetModel();

      // fix first stage variables
      // maybe fix only integer variables?
      for (unsigned col = 0; col < nvar1; col++) {
	scen_wrap[localscen-1]->setColLower(col,fsSolutions[soln*nvar1+col]);
	scen_wrap[localscen-1]->setColUpper(col,fsSolutions[soln*nvar1+col]);
      }

      // sanity check to print lb
      /*
      const double *lb = scen_wrap[localscen-1]->getColLower();
      cout << "Printing set lower bound" << endl;
      copy(lb, lb+nvar1, ostream_iterator<double>(cout, " "));
      cout << endl;
      */

      // Solving directly with Cbc
      cbcModel->branchAndBound();

      // initialize to NULL
      solVec = NULL;
      // get best solution
      solVec = cbcModel->bestSolution();
      // abort if no solution
      // Can't about without an AllReduce Call, exit loop
      if (solVec == NULL) {
	BBSMPS_ALG_LOG_SEV(debug) << "Found no feasible solution";
	found_sol = 0;
	break;
      }

      // TODO: Use the solutions we still have

      // Print first and second stage solution vector
      /*
      cout << "First and second Stage solution" << endl;
      copy(solVec, solVec+nvar1+nvar2, ostream_iterator<double>(cout, " "));
      cout << endl;
      */

      // put second stage part of solVec in soln position in ssSolutions vector.
      memcpy(&ssSolutions[localscen-1][soln*nvar2],
	     &solVec[nvar1], nvar2 * sizeof(double));

      // Second stage solution vector
      /*
      cout << "Second stage solution for scenario " << localscen << endl;
      copy(ssSolutions[localscen-1].begin(), ssSolutions[localscen-1].end(),
	   ostream_iterator<double>(cout, " "));
      cout << endl;
      */

      // Need to get it fresh after the solve, else it seems corrupted.
      // get objective coefficients
      objVec = scen_wrap[localscen-1]->getObjCoefficients();

      // Print objective vector
      /*
      cout << "Objective vector" << endl;
      copy(objVec, objVec+nvar1+nvar2, ostream_iterator<double>(cout, " "));
      cout << endl;
      */

      // Can get FS objective as long as we do it only once per solution
      // evaluate, dot product of solVec and objVec, to get objective contribution
      // of first (fixed) stage variables (only need to do first stage once)
      // for each solution
      // add to current ssObjvalues at end.
      if (localscen == 1) {
	fsObjvalues[soln] = inner_product(&objVec[0], &objVec[nvar1],
					  &fsSolutions[soln*nvar1],
					  fsObjvalues[soln]);
      }

      // First stage objective value vector
      /*
      cout << "First stage objective values after scenario/soln " << localscen
	   << " " << soln << endl;
      copy(fsObjvalues.begin(), fsObjvalues.end(),
	   ostream_iterator<double>(cout, " "));
      cout << endl;
      */

      // evaluate, dot product of solVec and objVec, to get objective contribution
      // of second stage variables
      // for each solution
      // add to current ssObjvalues

      ssObjvalues[soln] = inner_product(&objVec[nvar1], &objVec[nvar1+nvar2],
					&solVec[nvar1], ssObjvalues[soln]);

      // Second stage objective value vector
      /*
      cout << "Second stage objective values after scenario/soln " << localscen
	   << " " << soln << endl;
      copy(ssObjvalues.begin(), ssObjvalues.end(),
	   ostream_iterator<double>(cout, " "));
      cout << endl;
      */

    }

  }

  // MPI allreduce call to check if anybody failed.
  MPI_Allreduce(MPI_IN_PLACE, &found_sol, 1,
		MPI_INT, MPI_MIN,
		ctx.comm());

  // TODO: Don't need to abort, can just discard solution

  // Abort if some rank found no solution
  if (!found_sol) {

    BBSMPS_ALG_LOG_SEV(debug) << "Some FS solution had no recourse. Aborting";
    cumulativeTime+=(MPI_Wtime()-startTimeStamp);
    return success;

  }

  BBSMPS_ALG_LOG_SEV(debug) << "Now all solutions are fully populated";

  /*
  // Joint second stage solution vector pre all gather for all solutions
  cout << "Second stage solution before allgather for scenario" << endl;
  copy(ssSolutions[localscen-1].begin(), ssSolutions[localscen-1].end(),
       ostream_iterator<double>(cout, " "));
  cout << endl;
  */

  // TODO: Can I use getObjValue() from OsiCbc?

  BBSMPS_ALG_LOG_SEV(debug) << "Combining objective value from all scenarios";

  // Joint second stage objective vector pre all reduce for all solutions
  /*
  cout << "Second stage objective before allreduce for all ranks" << endl;
  copy(ssObjvalues.begin(), ssObjvalues.end(),
       ostream_iterator<double>(cout, " "));
  cout << endl;
  */

  // For each solution, collect second stage objective contribution
  // MPI_Allreduce, where you add each elemint

  // send values to all, which sum them all up.
  MPI_Allreduce(MPI_IN_PLACE, &ssObjvalues[0], nproc,
		MPI_DOUBLE, MPI_SUM,
		ctx.comm());

  // Joint second stage objective vector post all reduce for all solutions
  /*
  cout << "Second stage objective after allreduce for all ranks" << endl;
  copy(ssObjvalues.begin(), ssObjvalues.end(),
       ostream_iterator<double>(cout, " "));
  cout << endl;
  */

  BBSMPS_ALG_LOG_SEV(debug) << "Now adding first stage contributions for each solution";

  for(unsigned soln = 0; soln < nproc; soln++) {
    // sets to large value if infeasible
    ssObjvalues[soln] = (foundSol[soln]==1) ? ssObjvalues[soln] + fsObjvalues[soln]
      : DBL_MAX;
  }

  // Objective vector for all solutions

  /*
  cout << "Objective values for all solutions" << endl;
  copy(ssObjvalues.begin(), ssObjvalues.end(),
       ostream_iterator<double>(cout, " "));
  cout << endl;
  */

  // Pick best overall solution add to solution pool
  // Pick the best.
  vector<double>::iterator smallest
    = min_element(ssObjvalues.begin(), ssObjvalues.end());

  // index
  unsigned bestsol = smallest - ssObjvalues.begin();
  BBSMPS_ALG_LOG_SEV(debug) << "Best solution is solution: "
					   << bestsol << " with value: "
					   << ssObjvalues[bestsol];


  // Add
  // need to construct densebavector from fssolutions[bestsol]
  // and sssolutions from all ranks.
  // do not need mpi call since populating a denseBAvector

  BBSMPS_ALG_LOG_SEV(debug) << "Populating first stage solution";

  // allocating memory for primal solution
  ubPrimalSolution.allocate(dimsOriginalSlacks, ctx, PrimalVector);
  // dummy slack vector
  vector<double> firstslacks;
  firstslacks.assign(ncon1,0);

  // first stage variables
  // pick appropriate column from fsSolutions
  ubPrimalSolution.getFirstStageVec().copyBeginning(&fsSolutions[bestsol*nvar1],
						    nvar1);
  // assign slacks of 0
  ubPrimalSolution.getFirstStageVec().copyToPosition(&firstslacks[0],nvar1,ncon1);

  // print
  // ubPrimalSolution.getFirstStageVec().print();

  BBSMPS_ALG_LOG_SEV(debug) << "Populating second stage solution";

  // second stage
  for(unsigned localscen = 1; localscen < nlocalscen; localscen++) {

    int scen = localScen[localscen];
    int nvar2 = input.nSecondStageVars(scen);
    int ncon2 = input.nSecondStageCons(scen);

    BBSMPS_ALG_LOG_SEV(debug) << "Local/actual scenario: " << localscen << " "
					     << scen;

    // dummy slack vector
    vector<double> secondslacks;
    secondslacks.assign(ncon2,0);

    // pick appropriate column from ssSolutions
    ubPrimalSolution.getSecondStageVec(scen).copyBeginning(&ssSolutions[localscen-1][bestsol*nvar2], nvar2);

    // second stage slacks. What should slacks be? 0? setting to secondslacks
    ubPrimalSolution.getSecondStageVec(scen).copyToPosition(&secondslacks[0],
							    nvar2,ncon2);

    // ubPrimalSolution.getSecondStageVec(scen).print();

  }

  // For sanity check, one could add these values to rootSolver and do an LP solve

  // TO DO: Check for feasibility - what if there is no recourse

  // improving solution
  if(ssObjvalues[bestsol]<objUB){
    // check for feasibility
    if (isLPIntFeas(ubPrimalSolution)){
      BBSMPSSolution sol(ubPrimalSolution,ssObjvalues[bestsol]);
      sol.setTimeOfDiscovery(BBSMPSSolver::instance()->getWallTime());
      BBSMPSSolver::instance()->addSolutionToPool(sol);
    }
  }

  //return success
  success= (ssObjvalues[bestsol]<objUB);
  timesSuccessful+=success;

  if (0 == mype && success)
    BBSMPS_ALG_LOG_SEV(warning)
      << "The scenario decomposition heuristic was successful.";

  cumulativeTime+=(MPI_Wtime()-startTimeStamp);

  return success;

}

bool BBSMPSHeuristicScenDecom::shouldItRun(BBSMPSNode* node, denseBAVector &nodeSolution){

  BBSMPS_ALG_LOG_SEV(debug)
    << "Should the scenario decomposition heuristic run?";

  // Run if periodicity says okay.
  return true;

}

