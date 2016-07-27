// ----------------------------------------------------------------------------
/**
   File: BBSMPSMaxFracSmallestScenarioFirstBranchingRule.hpp

   Description: Branching rule that branches on the maximum fractional variable first, of the smallest scenario. 


*/ 
// ----------------------------------------------------------------------------

#ifndef BBSMPSMAXFRACSMALLESTSCENARIOFIRSTBRANCHINGRULE_H
#define BBSMPSMAXFRACSMALLESTSCENARIOFIRSTBRANCHINGRULE_H

#include "BBSMPSBranchingRule.hpp"
#include "BBSMPSUtils.hpp"
#include "BBSMPSSolver.hpp"
#include "BBSMPSLogging.hpp"

class BBSMPSMaxFracSmallestScenarioFirstBranchingRule: public BBSMPSBranchingRule {

public:

   	virtual bool branch(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes, const denseBAVector& primalSoln);
   	
   	BBSMPSMaxFracSmallestScenarioFirstBranchingRule(int priority): BBSMPSBranchingRule(priority){name="Max Fractional Branching Rule, smallest scenario first";};


private:

	// Auxiliary functions for branching
   	int getFirstStageMinIntInfeasCol( const denseBAVector& primalSoln,  SMPSInput& input);

   	int getFirstStageMaxFint( const denseBAVector& primalSoln,  SMPSInput& input) ;

   	int getFirstStageMaxFracPartCol( const denseBAVector& primalSoln,  SMPSInput& input);

      int getSecondStageMaxFracPartCol( const denseBAVector& primalSoln,  SMPSInput& input, int scen);

      int getNSecondStageFracVars( const denseBAVector& primalSoln,  SMPSInput& input, int scen);

	// NOTE: MPI standard requires passing ints, not bools
   	int isFirstStageIntFeas( const denseBAVector& primalSoln,  SMPSInput& input) ;

   	void branchOnFirstStage(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes,  const denseBAVector& primalSoln,  SMPSInput& input) ;

   	int getSecondStageMinIntInfeasCol( const denseBAVector& primalSoln, int scen,  SMPSInput& input);

   	int isSecondStageIntFeas( const denseBAVector& primalSoln, int scen, SMPSInput & input) ;

   	void branchOnSecondStage(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes, const denseBAVector& primalSoln,  SMPSInput& input,BAContext &ctx,int mype) ;

};

#endif

