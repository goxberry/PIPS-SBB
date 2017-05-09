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
   File: BBSMPSHeuristicScenDecom.hpp

   Description: ScenDecom PRIMAL HEURISTIC: Single scenario problems are solved
   to obtain feasible solutions. Single scenario MIPS solved using CBC

   Initially, each MPI rank will only solve for one scenario in its rank.
   May solve all scenarios in its rank in the future.

   The set of first stage feasible solutions will be distributed among all MPI ranks.
   Each rank will solve each first stage solution for all its second stage problems.

   Each problem may be optimized up to a certain time limit and/or node limit.
   In the future, only some first stage solutions may be evaluated.

   Initially, will only run at root node, and only once.
   Future versions may run multiple times, and with lambda updates between runs.

*/
// ----------------------------------------------------------------------------


#ifndef BBSMPSHEURISTICSCENDECOM_H
#define BBSMPSHEURISTICSCENDECOM_H

// Needed for Cbc interfaces
#include "OsiCbcSolverInterface.hpp"

#include "BBSMPSNode.hpp"
#include "BBSMPSHeuristic.hpp"
#include "BBSMPSUtils.hpp"

#include <numeric>

class BBSMPSHeuristicScenDecom: public BBSMPSHeuristic {

public:
  BBSMPSHeuristicScenDecom(int offset, int depth, const char *_name,
			   int _nodeLim, double _timeLim, bool _sameScen):
    nodeLim(_nodeLim), timeLim(_timeLim), sameScen(_sameScen),
    firstTime(true),
    local_scen_num(1),
    BBSMPSHeuristic(offset,depth,_name){};
  bool runHeuristic(BBSMPSNode* node, denseBAVector &LPRelaxationSolution);
  bool shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution);

private:

  int    nodeLim;  // nodes to run B&B on single scen problems
  double timeLim; // timelimit on single scen problem
  bool   sameScen;  // solve all scenarios in rank? Currently not used.

  bool   firstTime; // indicator for first call

  int    local_scen_num;  // which local scenario to consider

  std::vector<OsiCbcSolverInterface> scen_wrap;
  // first stage solutions for each rank
  std::vector<double> fsSolutions;
  // second stage solution for each local scenario (outer)
  // for each fs solution (inner)
  std::vector<std::vector<double> > ssSolutions;
  // vector of objective value for second stage solutions for each rank
  std::vector<double> ssObjvalues;
  // vector of objective value for first stage solutions for each rank
  std::vector<double> fsObjvalues;
  // primal solution vector to hold best solution
  denseBAVector ubPrimalSolution;

};


#endif

