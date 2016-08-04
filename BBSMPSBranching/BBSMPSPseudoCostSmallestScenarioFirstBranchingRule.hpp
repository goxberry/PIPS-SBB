// ----------------------------------------------------------------------------
/**
   File: BBSMPSPseudoCostSmallestScenarioFirstBranchingRule.hpp

   Description: Branching rule that branches on the maximum fractional variable first. 


*/ 
// ----------------------------------------------------------------------------

#ifndef BBSMPSPSEUDOCOSTSMALLESTSCENARIOFIRSTBRANCHINGRULE_H
#define BBSMPSPSEUDOCOSTSMALLESTSCENARIOFIRSTBRANCHINGRULE_H

#include "BBSMPSBranchingRule.hpp"
#include "BBSMPSUtils.hpp"
#include "BBSMPSSolver.hpp"
#include "BBSMPSLogging.hpp"

class BBSMPSPseudoCostSmallestScenarioFirstBranchingRule: public BBSMPSBranchingRule {

public:

   	virtual bool branch(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes, const denseBAVector& primalSoln);
   	
   	BBSMPSPseudoCostSmallestScenarioFirstBranchingRule(int priority);

    void getCostHistory(denseBAVector &_downPseudoCost,denseBAVector &_upPseudoCost,denseBAVector &_downBranchingHistory,denseBAVector &_upBranchingHistory){
      _downPseudoCost=downPseudoCost;
      _upPseudoCost=upPseudoCost;
      _downBranchingHistory=downBranchingHistory;
      _upBranchingHistory=upBranchingHistory;
    }

    void setCostHistory(denseBAVector &_downPseudoCost,denseBAVector &_upPseudoCost,denseBAVector &_downBranchingHistory,denseBAVector &_upBranchingHistory){
      downPseudoCost=_downPseudoCost;
      upPseudoCost=_upPseudoCost;
      downBranchingHistory=_downBranchingHistory;
      upBranchingHistory=_upBranchingHistory;
    }

private:

	// Auxiliary functions for branching
   	int getFirstStageMinIntInfeasCol( const denseBAVector& primalSoln,  SMPSInput& input);

   	int getFirstStageMaxFint( const denseBAVector& primalSoln,  SMPSInput& input) ;

   	int getFirstStageMaxFracPartCol( const denseBAVector& primalSoln,  SMPSInput& input);

	// NOTE: MPI standard requires passing ints, not bools
   	int isFirstStageIntFeas( const denseBAVector& primalSoln,  SMPSInput& input) ;

   	void branchOnFirstStage(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes,  const denseBAVector& primalSoln,  SMPSInput& input) ;

   	int getSecondStageMinIntInfeasCol( const denseBAVector& primalSoln, int scen,  SMPSInput& input);

   	int isSecondStageIntFeas( const denseBAVector& primalSoln, int scen, SMPSInput & input) ;

   	void branchOnSecondStage(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes, const denseBAVector& primalSoln,  SMPSInput& input,BAContext &ctx,int mype, int scen) ;

    void initializeVariable(const denseBAVector &nodeRelaxation,BAFlagVector<variableState> &warmstart, double lpRelaxationObjValue, denseBAVector &lb, denseBAVector &ub, int scen, int col);

    bool performRoundOfFirstStageInitializations(const denseBAVector &nodeRelaxation,BAFlagVector<variableState> &warmstart, denseBAVector &lb, denseBAVector &ub, double lpRelaxationObjValue);

    bool performRoundOfSecondStageInitializations(const denseBAVector &nodeRelaxation,BAFlagVector<variableState> &warmstart, denseBAVector &lb, denseBAVector &ub, double lpRelaxationObjValue, int scen);

    denseBAVector downPseudoCost;
    denseBAVector upPseudoCost;
    denseBAVector downBranchingHistory;
    denseBAVector upBranchingHistory;
    bool everythingFirstStageInitialized;
    vector<bool> everythingSecondStageInitialized;
    int reliabilityFactor;

};

#endif

