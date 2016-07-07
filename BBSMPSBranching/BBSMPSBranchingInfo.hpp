// ----------------------------------------------------------------------------
/**
   File: BBSMPSBranchingInfo.hpp

   Description: Container class for storing the information relative to a single branching
   				decision.

*/ 
// ----------------------------------------------------------------------------


#ifndef BBSMPSBRANCHINGINFO_H
#define BBSMPSBRANCHINGINFO_H

#include "BAVector.hpp"
class BBSMPSBranchingInfo {

public:
	BBSMPSBranchingInfo(int _varNumber, double _bound, char _direction, int _stage, int _scenario=-1);
	BBSMPSBranchingInfo(int *intVector, double *dblVector);
	BBSMPSBranchingInfo();
	~BBSMPSBranchingInfo();

	int getVarNumber();
	double getBound();
	char getDirection();
	int getStageNumber();
	int getScenarioNumber();
	void applyBranchingInfo(denseBAVector &lb,denseBAVector &ub);
	void serialize(int *intVector, double *dblVector);
	static void getSerializationSize(int &intVectorSize, int &dblVectorSize);

private:
	int varNumber;
	double bound;
	char direction;
	int stage;
	int scenario;
	
};

#endif