// ----------------------------------------------------------------------------
/**
   File: BBSMPSMinFracSmallestScenarioFirstBranchingRule.hpp

   Description: Branching rule that branches on the minimum fractional variable first, of the smallest scenario. 


*/ 
// ----------------------------------------------------------------------------

#ifndef BBSMPSMINFRACSMALLESTSCENARIOFIRSTBRANCHINGRULE_H
#define BBSMPSMINFRACSMALLESTSCENARIOFIRSTBRANCHINGRULE_H

#include "BBSMPSBranchingRule.hpp"
#include "BBSMPSUtils.hpp"
#include "BBSMPSSolver.hpp"
#include "BBSMPSLogging.hpp"

class BBSMPSMinFracSmallestScenarioFirstBranchingRule: public BBSMPSBranchingRule {

public:

   	virtual bool branch(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes, const denseBAVector& primalSoln);
   	
   	BBSMPSMinFracSmallestScenarioFirstBranchingRule(int priority): BBSMPSBranchingRule(priority){name="Min Fractional Branching Rule, smallest scenario first";};


private:

	// Auxiliary functions for branching
   	int getFirstStageMinIntInfeasCol( const denseBAVector& primalSoln,  SMPSInput& input);

   	int getFirstStageMaxFint( const denseBAVector& primalSoln,  SMPSInput& input) ;

   	int getFirstStageMinFracPartCol( const denseBAVector& primalSoln,  SMPSInput& input);

      int getSecondStageMinFracPartCol( const denseBAVector& primalSoln,  SMPSInput& input, int scen);

      int getNSecondStageFracVars( const denseBAVector& primalSoln,  SMPSInput& input, int scen);

	// NOTE: MPI standard requires passing ints, not bools
   	int isFirstStageIntFeas( const denseBAVector& primalSoln,  SMPSInput& input) ;

   	void branchOnFirstStage(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes,  const denseBAVector& primalSoln,  SMPSInput& input) ;

   	int getSecondStageMinIntInfeasCol( const denseBAVector& primalSoln, int scen,  SMPSInput& input);

   	int isSecondStageIntFeas( const denseBAVector& primalSoln, int scen, SMPSInput & input) ;

   	void branchOnSecondStage(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes, const denseBAVector& primalSoln,  SMPSInput& input,BAContext &ctx,int mype) ;

};

#endif

