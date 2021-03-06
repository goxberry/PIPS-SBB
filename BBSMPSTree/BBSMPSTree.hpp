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
   File: BBSMPSTree.hpp

   Description: Header file for the tree class. The tree class in the vertebral piece
                of the algorithm. It handles and coordinates the execution of the branch
                and bound. It incudes:
                - The node tree, materialized in a heap of BBSMPSNode pointers.
                - The branching rule manager
                - The heuristics manager

   Limitations: Due to design restrictions, the tree shouldn't access the internals of any
                object, but merely act as the coordinator/manager.

*/
// ----------------------------------------------------------------------------



#ifndef BBSMPSTREE_H
#define BBSMPSTREE_H

#include "SMPSInput.hpp"
#include "BAData.hpp"
#include "PIPSSInterface.hpp"
#include "BBSMPSNode.hpp"
#include <boost/scoped_ptr.hpp>
#include <cstdlib>

// For branch-and-bound code:
#include <queue> // priority queue
#include <cassert> // C-style assertions
#include <cmath> // for floor, ceil, abs functions
#include <algorithm> // for min
#include <set>
#include <utility>



#include "BBSMPSMaxFracBranchingRule.hpp"
#include "BBSMPSBranchingRuleManager.hpp"
#include "BBSMPSHeuristicsManager.hpp"
#include "BBSMPSHeuristicRounding.hpp"
#include "BBSMPSHeuristicRENS.hpp"
#include "BBSMPSUtils.hpp"
#include "BBSMPSSolver.hpp"
#include "BBSMPSSolverState.hpp"
#include "BBSMPSLogging.hpp"
#include "BBSMPSSolution.hpp"
#include "BBSMPSHeuristicRINS.hpp"
#include "BBSMPSHeuristicFixAndDive.hpp"
#include "BBSMPSHeuristicFixAndDiveLocks.hpp"
#include "BBSMPSHeuristicCrossover.hpp"
#include "BBSMPSHeuristicLockRounding.hpp"
#include "BBSMPSHeuristicMagic.hpp"
#include "BBSMPSHeuristicSolutionRINS.hpp"
#include "BBSMPSHeuristicBestRINSJump.hpp"
#include "BBSMPSHeuristicSolutionPolishing.hpp"
#include "BBSMPSPseudoCostBranchingRule.hpp"
#include "BBSMPSCuttingPlane.hpp"
#include "BBSMPSCuttingPlaneGenerator01KP.hpp"
#include "BBSMPSCuttingPlanesManager.hpp"
// Outputs solver status:
void outputLPStatus(solverState lpStatus);

// Overload the "less than" operator so that priority_queue can use it
// for heapifying comparisons. Make BBSMPSNode a templated
// class when adding additional branching heuristics?
bool operator< (const BBSMPSNode& left,
  const BBSMPSNode& right);

bool operator> (const BBSMPSNode& left,
  const BBSMPSNode& right);


class nodePtrComparison
{
public:
  bool operator() (const BBSMPSNode* lhs, const BBSMPSNode* rhs) const
  {
  if(lhs->getParentObjective() == rhs->getParentObjective()){
    return (lhs->getNodeNumber() < rhs->getNodeNumber());
  }
  return (lhs->getParentObjective() > rhs->getParentObjective());

 }
};


enum nodeSelectionRule { BestBound, DepthFirst };

// TODO: Check if PIPS-S always minimizes; otherwise, must change logic.

class BBSMPSTree {


public:
  // Upon constructing the B&B tree:
  // - instantiate MPI communicator context for block angular objects
  // - Get rank of process relative to ctx.comm()
  // - instantiate input object
  // - instantiate PIPS-S for root LP relaxation
  // - instantiate dimensions object for holding problem
  //      dimension information for allocating vectors...
  // - instantiate objective function upper bound to +infty (or a value close to that)
  // - set integrality tolerance
  // - set optimality gap tolerances
  // - set LP solver tolerances
  // - set comparison tolerance to primal tolerance for now
  // - set solver status to "LoadedFromFile" because this interface forces MILP to
  //   be loaded from an SMPS file

  BBSMPSTree(const SMPSInput& smps);

  BBSMPSTree(BBSMPSNode *node, double lb=COIN_DBL_MIN, double ub=COIN_DBL_MAX);

  // Default destructor
  ~BBSMPSTree();

  void setTimeLimit(int _tiLim);

  void setNodeLimit(int _nodeLim);

  void setSolLimit(int _solLim);

  void setGAPTolLimit( double _GAPTolLim);

  bool retrieveBestSolution(BBSMPSSolution &solution);

  void branchAndBound();

  void loadSimpleHeuristics();

  void loadMIPHeuristics();

  void loadCuttingPlanes();

  void setVerbosity(bool verbose);

  void removeCuts();


void loadLPHeuristic(BBSMPSHeuristic *heur);
  void loadMIPHeuristic(BBSMPSHeuristic *heur);
  void generateIncrementalWarmState(BBSMPSNode* node, const BAFlagVector<variableState> & originalState, const BAFlagVector<variableState> &currentState);


  BBSMPSNode* topOfHeap();

  static BBSMPSNode* getRootNode(){
    return rootNode;
  }


  void setLB(double lb){ objLB=lb;};
  void setUB(double ub){ objUB=ub;};


private:


    // Node selection rule, determines which node is chosen next. By default, it is bestbound,
  // by using the min-heap. To handle other rules, we will need to refactor priority_queue to vector,
  // with make_heap, etc. (Except for depth-first, which doesn't need to make_heap)
  nodeSelectionRule nodesel;

  double objUB; // best upper bound on objective function value
  denseBAVector ubPrimalSolution; // primal solution for best UB on obj

  double objLB; // best lower bound on objective function value

  //double intTol; // tolerance on integrality checks
  double optGapTol; // tolerance on optimality gap between UB and LB.
  double lpPrimalTol; // tolerance on LP primal problems
  double lpDualTol; // tolerance on LP dual problems
  double compTol; // tolerance on LP objective function comparisons

  int bbIterationCounter;

  int solsDiscoveredInit;
  int solsDiscoveredLimit;
  int tiLim;

  int nodeLim;
  double LPRelaxationTime;
  double PreProcessingTime;
  double LPRelaxationValue;

  int nodesFathomed;
  int nodesBecameInteger;
  bool verbosityActivated;

  vector<int> currentlyAppliedPlanes;
  BBSMPSBranchingRuleManager branchingRuleManager;
  BBSMPSHeuristicsManager heuristicsManager;
  BBSMPSCuttingPlanesManager cuttingPlanesManager;
  // max-heap data structure with nodes
  // TODO: Refactor to vector<BranchAndBound> & replace w/ make_heap, push_heap, pop_heap
  //std::priority_queue<BBSMPSNode, std::vector<BBSMPSNode>, std::less<BBSMPSNode> > heap; // max-heap
  std::priority_queue<BBSMPSNode*, std::vector<BBSMPSNode*>, nodePtrComparison > heap; // min-heap

  // Solver status; can only be in the set {LoadedFromFile, Initialized,
  // PrimalFeasible, Optimal, ProvenUnbounded, ProvenInfeasible, Stopped}
  // because there is no duality theory, and currently, the only interface
  // to the solver has to load problem data from a file as one of the first
  // steps.
  BBSMPSSolverState status;

  static BBSMPSNode* rootNode;

    // Auxiliary functions for branching
  int getFirstStageMinIntInfeasCol(const denseBAVector& primalSoln);

  void setStatusToPrimalFeasible() ;

  void setStatusToOptimal() ;

  void setStatusToProvenInfeasible();

  void setStatusToStopped();




  // Make default constructor impossible to call.
  BBSMPSTree();

} ;
#endif
