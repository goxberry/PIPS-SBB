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
// ----------------------------------------------------------------------------
/**
   File: BBSMPSSolver.hpp

   Description: Singleton class that manages the solver and the problem model. It
                contains:
                - The PIPSSInterface
                - The SMPSInput related to the problem model
                - The MPI context
                - The original lower and upper bounds for the problem variables
                - Problem dimensions and slacks

   Limitations: When accessing the solver data, please use the references with caution.
                To instantiate and initialize:

                BBSMPSSolver::initialize(smps); (Where smps is a SMPSInput object)

                Accessing examples:

                PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
                BAContext &ctx=BBSMPSSolver::instance()->getBAContext();

*/
// ----------------------------------------------------------------------------


#ifndef BBSMPSSOLVER_H
#define BBSMPSSOLVER_H

#include <vector>

#include "SMPSInput.hpp"
#include "BAData.hpp"
#include "PIPSSInterface.hpp"
#include "Presolve.hpp"
#include "BBSMPSLogging.hpp"
#include "BBSMPSSolution.hpp"
#include <set>

struct solutionComparison {
  bool operator() (const BBSMPSSolution& lhs, const BBSMPSSolution& rhs) const
  {return (lhs.getObjValue()<rhs.getObjValue());}
};

class BBSMPSSolver {

public:

  ~BBSMPSSolver();

  BAContext& getBAContext();
  int& getMype();
  SMPSInput& getSMPSInput();
  PIPSSInterface& getPIPSInterface();
  BADimensions& getBADimensions();
  const BADimensionsSlacks& getOriginalBADimensionsSlacks();
  const BADimensionsSlacks& getBADimensionsSlacks();
  const denseBAVector& getOriginalLB();
  const denseBAVector& getOriginalUB();
  const BAFlagVector<variableState>& getOriginalWarmStart();
  void setOriginalWarmStart(BAFlagVector<variableState>&warmStart);
  denseBAVector& getLPRelaxation();
  void setLPRelaxation(denseBAVector &_LPRelaxation);
  void setLPRelaxationObjectiveValue(double lpRelObjVal);
  double getLPRelaxationObjectiveValue();
  static BBSMPSSolver *instance();
  static BBSMPSSolver *initialize(const SMPSInput &_input);
  static bool isInitialized();
  static void deInitialize();
  void printPresolveStatistics();
  void addSolutionToPool(BBSMPSSolution &sol);
  void printSolutionStatistics(double objLB);
  const BBSMPSSolution &getSoln(int index);
  const BBSMPSSolution &getSolnBySolNumber(int number);
  void commitNewColsAndRows();
  int getSolPoolSize();
  double getWallTime();
void resetSolver();

protected:

private:
  BBSMPSSolver(const SMPSInput &_input);
  BAContext ctx; // MPI communication context for PIPS-S
  int mype; // MPI rank of process storing tree (relative to comm in ctx)
  SMPSInput input; // SMPS input file for reading in block angular MILP
   BAData problemData;
  PIPSSInterface *rootSolver; // PIPS-S instance for root LP relaxation
   // Presolve pre;
  BADimensions dims; // Dimension object for instantiating warm start information
  BADimensionsSlacks originalDims;
  BADimensionsSlacks dimsSlacks; // Dimension object for warm start info
  denseBAVector lb;
  denseBAVector ub;
  denseBAVector lbModifiableWithCuts;
  denseBAVector ubModifiableWithCuts;
  denseBAVector LPRelaxation;
  BAFlagVector<variableState> originalWarmStart;
  BAFlagVector<variableState> originalWarmStartModifiableWithCuts;
  std::set<BBSMPSSolution,solutionComparison> solutionPool;

  double LPRelaxationObjectiveValue;
  double startTimeStamp;
  static BBSMPSSolver *solverInstance;
};

#endif
