// ----------------------------------------------------------------------------
/**
   File: BBSMPSBranchingRule.hpp

   Description: Base virtual class of a branching rule.

   Limitations: An object of this class should not be instantiated, it should always
   be subclassed. A rule must be initialized with a priority.

*/ 
// ----------------------------------------------------------------------------

#ifndef BBSMPSBRANCHINGRULE_H
#define BBSMPSBRANCHINGRULE_H

#include "BBSMPSNode.hpp"
#include "BBSMPSLogging.hpp"
#include "BAData.hpp"

typedef struct comp
{
  double result;
  int rank;
} comp;

typedef struct compIntInt
{
  int result;
  int rank;
} compIntInt;

class BBSMPSBranchingRule {

public:
	BBSMPSBranchingRule(int _priority);
	~BBSMPSBranchingRule();
	int getPriority() const;
	void setPriority(int _priority);
	virtual bool branch(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes,const  denseBAVector& primalSoln){std::cout<<"Well, this is an error...\n";};
	virtual void printStatistics();
	void setAccumulatedRunTime(double time){accumulatedRunTime=time;};
	double getAccumulatedRunTime(){return accumulatedRunTime;};
private:
	int priority;

protected:
	int timesCalled;
	int timesSuccessful;
	double accumulatedRunTime;
	std::string name;


};

#endif
