// ----------------------------------------------------------------------------
/**
   File: BBSMPSMinFracBranchingRule.hpp

   Description: Branching rule that branches on the minimum fractional variable first. 


*/ 
// ----------------------------------------------------------------------------

#ifndef BBSMPSMINFRACBRANCHINGRULE_H
#define BBSMPSMINFRACBRANCHINGRULE_H

#include "BBSMPSBranchingRule.hpp"
#include "BBSMPSUtils.hpp"
#include "BBSMPSSolver.hpp"
#include "BBSMPSLogging.hpp"



class BBSMPSMinFracBranchingRule: public BBSMPSBranchingRule {

public:

   	virtual bool branch(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes, const denseBAVector& primalSoln);
   	
   	BBSMPSMinFracBranchingRule(int priority): BBSMPSBranchingRule(priority){name="Min Fractional Branching Rule";};


private:

	// Auxiliary functions for branching
   	int getFirstStageMinIntInfeasCol( const denseBAVector& primalSoln,  SMPSInput& input);

   	int getFirstStageMinFint( const denseBAVector& primalSoln,  SMPSInput& input) ;

   	int getFirstStageMinFracPartCol( const denseBAVector& primalSoln,  SMPSInput& input);

      int getSecondStageMinFracPartCol( const denseBAVector& primalSoln,  SMPSInput& input, int scen);


	// NOTE: MPI standard requires passing ints, not bools
   	int isFirstStageIntFeas( const denseBAVector& primalSoln,  SMPSInput& input) ;

   	void branchOnFirstStage(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes,  const denseBAVector& primalSoln,  SMPSInput& input) ;

   	int getSecondStageMinIntInfeasCol( const denseBAVector& primalSoln, int scen,  SMPSInput& input);

   	int isSecondStageIntFeas( const denseBAVector& primalSoln, int scen, SMPSInput & input) ;

   	void branchOnSecondStage(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes, const denseBAVector& primalSoln,  SMPSInput& input,BAContext &ctx,int mype) ;

};

#endif

