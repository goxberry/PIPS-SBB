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

  BBSMPS_ALG_LOG_SEV(warning) << "Problem size: " << nvar1 << " " << ncon1 << " " << nscen;

  BBSMPS_ALG_LOG_SEV(debug) << "MPI info: " << nproc << " " << mype;
  if (0 == mype) {
    BBSMPS_ALG_LOG_SEV(warning) << "Performing the ScenDecom Heuristic.";
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

  // TODO: Resize vector to the size it will have

  // TODO: Can I use existing solution pool to provide incumbents to CBC?

  // Useful to access CbcModel directly
  CbcModel *cbcModel;

  // First time entering the heuristic
  if (firstTime) {

    for(unsigned localscen = 1; localscen < nlocalscen; localscen++) {

      int scen = localScen[localscen];

      // uses the default constructor
      OsiCbcSolverInterface wrap(NULL,NULL);

      // should I be using calls from rootSolver instead. Yes, if there is presolve.
      // To be revisited.
      // Similarly OSi problem wrapper could be given the presolved problem.
      // many dimensions could break if there is presolve.
      // will have to stop using input.

      makeOneScenOsiSolverInterface(scen, input, &wrap);

      // Add node and time limits
      wrap.setMaximumNodes(nodeLim);
      wrap.setMaximumSeconds(timeLim);

      // Add to queue.
      scen_wrap.push_back(wrap);

    }

  }

  // can set firstTime = false. scen_wrap is now initialized.
  firstTime = false;

  for(unsigned localscen = 1; localscen < nlocalscen; localscen++) {

    int scen = localScen[localscen];

    int nvar2 = input.nSecondStageVars(scen);
    int ncon2 = input.nSecondStageCons(scen);

    if (localscen==1) {
      BBSMPS_ALG_LOG_SEV(warning) << "Problem size: " << nvar2 << " " << ncon2 << " " << scen;
    }

    // modify bounds based on node bounds

    // Print original lb
    if (localscen==1) {

      const double *testlb = scen_wrap[localscen-1].getColLower();
      cout << "Printing previous lower bound" << endl;
      copy(testlb, testlb+nvar1+nvar2, ostream_iterator<double>(cout, " "));
      cout << endl;

      cout << "Printing node lb, including slacks" << endl;
      lb.getFirstStageVec().print();
      lb.getSecondStageVec(scen).print();

    }

    // first stage
    for (unsigned col = 0; col < nvar1; col++) {
      scen_wrap[localscen-1].setColLower(col,lb.getFirstStageVec()[col]);
      scen_wrap[localscen-1].setColUpper(col,ub.getFirstStageVec()[col]);
    }

    // second stage
    for (unsigned col = 0; col < nvar2; col++) {
      scen_wrap[localscen-1].setColLower(nvar1+col,lb.getSecondStageVec(scen)[col]);
      scen_wrap[localscen-1].setColUpper(nvar1+col,ub.getSecondStageVec(scen)[col]);
    }

    // Print modified lb
    if (localscen==1) {
      const double *modlb = scen_wrap[localscen-1].getColLower();
      cout << "lb after modifying" << endl;
      copy(modlb, modlb+nvar1+nvar2, ostream_iterator<double>(cout, " "));
      cout << endl;
    }

    // Important note: Need to offset scen_wrap by -1 since localscen starts at 1

  }

  BBSMPS_ALG_LOG_SEV(debug) << "Done building single scenario objects";

  // size vector of first stage solutions for allgather
  // one per rank, each of length nvar1
  fsSolutions.resize(nproc*nvar1);

  // For one scenario in rank, run CBC with settings

  // increment if no sameScen
  if (!sameScen) {
    local_scen_num++; // increment
    if (local_scen_num == nlocalscen) local_scen_num = 1; // handle overflow
  }

  // indicator for finding feasible solution
  int found_sol = 1;

  // Write mps; just for debugging
  //  scen_wrap[local_scen_num-1].writeMps("temp","mps");

  BBSMPS_ALG_LOG_SEV(debug) << "Local/actual scenario: "
					   << local_scen_num << " "
					   << localScen[local_scen_num];

  // get pointer to underlying model, easier to query
  cbcModel = scen_wrap[local_scen_num-1].getModelPtr();
  cbcModel->setLogLevel(0);

  // Just to be safe, resetting model solve
  cbcModel->resetModel();

  // Solving directly with Cbc
  cbcModel->branchAndBound();

  // initialize to NULL
  solVec = NULL;
  // get best solution
  solVec = cbcModel->bestSolution();
  // abort if no solution, mark as no solution
  if (solVec == NULL) {
    cout << "Found no feasible solution" << endl;
    found_sol = 0;
  }

  // Print first stage solution vector
  /*
  cout << "First stage solution" << endl;
  copy(solVec, solVec+nvar1, ostream_iterator<double>(cout, " "));
  cout << endl;
  */

  // MPI allreduce call to check if anybody failed.
  MPI_Allreduce(MPI_IN_PLACE, &found_sol, 1,
		MPI_INT, MPI_MIN,
		ctx.comm());

  // TODO: Continue as long as someone has a solution
  // TODO: Be able to continue with fewer/more solutions than ranks
  // TODO: Maybe only solve some of them by guessing which one is best

  // Abort if some rank found no solution
  if (!found_sol) return success;

  // TODO: Solutions may not be unique

  // Joint first stage solution vector pre all gather
  /*
  cout << "First stage solution before allgather" << endl;
  copy(fsSolutions.begin(), fsSolutions.end(), ostream_iterator<double>(cout, " "));
  cout << endl;
  */

  // MPI call to share all first stage solutions.
  MPI_Allgather(&solVec[0], nvar1, MPI_DOUBLE,
		&fsSolutions[0], nvar1, MPI_DOUBLE,
		ctx.comm());

  // Joint first stage solution vector post all gather
  /*
  cout << "First stage solution after allgather" << endl;
  copy(fsSolutions.begin(), fsSolutions.end(), ostream_iterator<double>(cout, " "));
  cout << endl;
  */

  // Now, all ranks have all the solutions.

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

      BBSMPS_ALG_LOG_SEV(debug) << "Considering solution: " << soln;
      /*
      cout << "First stage solution being fixed" << endl;
      copy(fsSolutions.begin(), fsSolutions.end(), ostream_iterator<double>(cout, " "));
      cout << endl;
      */
      // get pointer to underlying model, easier to query
      cbcModel = scen_wrap[localscen-1].getModelPtr();
      cbcModel->setLogLevel(0);

      // Just to be safe, resetting model solve
      cbcModel->resetModel();

      // fix first stage variables
      // maybe fix only integer variables?
      for (unsigned col = 0; col < nvar1; col++) {
	scen_wrap[localscen-1].setColLower(col,fsSolutions[soln*nvar1+col]);
	scen_wrap[localscen-1].setColUpper(col,fsSolutions[soln*nvar1+col]);
      }

      // sanity check to print lb
      /*
      const double *lb = scen_wrap[localscen-1].getColLower();
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
	cout << "Found no feasible solution" << endl;
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
      objVec = scen_wrap[localscen-1].getObjCoefficients();

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
  if (!found_sol) return success;

  BBSMPS_ALG_LOG_SEV(debug) << "Now all solutions are fully populated";

  /*
  // Joint second stage solution vector pre all gather for all solutions
  cout << "Second stage solution before allgather for scenario" << endl;
  copy(ssSolutions[localscen-1].begin(), ssSolutions[localscen-1].end(),
       ostream_iterator<double>(cout, " "));
  cout << endl;
  */

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
    ssObjvalues[soln] += fsObjvalues[soln];
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

  // Turns out this is not very good solution. Maybe const scenarios?

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

