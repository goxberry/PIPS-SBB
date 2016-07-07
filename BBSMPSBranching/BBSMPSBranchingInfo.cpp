#include "BBSMPSBranchingInfo.hpp"

using namespace std;


BBSMPSBranchingInfo::BBSMPSBranchingInfo(int _varNumber, double _bound, char _direction, int _stage, int _scenario):
varNumber(_varNumber),
bound(_bound),
direction(_direction),
stage(_stage),
scenario(_scenario){}



BBSMPSBranchingInfo::BBSMPSBranchingInfo(int *intVector, double *dblVector)
{
	varNumber=intVector[0];
	stage=intVector[2];
	scenario=intVector[3];
	switch ( intVector[1] ) {
		case 1:
		  direction='L';
		  break;
		case 2:
		  direction='U';
		  break;
		default:
		  direction='E';
	}
}
	

BBSMPSBranchingInfo::~BBSMPSBranchingInfo(){

}

BBSMPSBranchingInfo::BBSMPSBranchingInfo(){

}

int BBSMPSBranchingInfo::getVarNumber(){
	return varNumber;
}
double BBSMPSBranchingInfo::getBound(){
	return bound;
}
char BBSMPSBranchingInfo::getDirection(){
	return direction;
}

int BBSMPSBranchingInfo::getStageNumber(){
	return stage;
}

int BBSMPSBranchingInfo::getScenarioNumber(){
	return scenario;
}

void BBSMPSBranchingInfo::applyBranchingInfo(denseBAVector &lb,denseBAVector &ub){
	if (getStageNumber()==1){
    	int varN=getVarNumber();
     	double bd=getBound();
	    if (getDirection()=='L'){
	    	lb.getFirstStageVec()[varN] = std::max(bd,lb.getFirstStageVec()[varN]);
	    }
	    else if (getDirection()=='U'){
	    	ub.getFirstStageVec()[varN] =  std::min(bd,ub.getFirstStageVec()[varN]);
	    }
	    else if(getDirection()=='E'){
	      lb.getFirstStageVec()[varN] = std::max(bd,lb.getFirstStageVec()[varN]);
	      ub.getFirstStageVec()[varN] =  std::min(bd,ub.getFirstStageVec()[varN]);
	    }
    }
   	else{
   		int scenario=getScenarioNumber();
	    int varN=getVarNumber();
	    double bd=getBound();
	    if (getDirection()=='L'){
	      lb.getSecondStageVec(scenario)[varN] = std::max(bd,lb.getSecondStageVec(scenario)[varN]);
	    }
	    else if (getDirection()=='U'){
	      ub.getSecondStageVec(scenario)[varN] =  std::min(bd,ub.getSecondStageVec(scenario)[varN]);
	    }
	    else if(getDirection()=='E'){
	      lb.getSecondStageVec(scenario)[varN] = std::max(bd,lb.getSecondStageVec(scenario)[varN]);
	      ub.getSecondStageVec(scenario)[varN] =  std::min(bd,ub.getSecondStageVec(scenario)[varN]);
	    }
	}
}

void BBSMPSBranchingInfo::serialize(int *intVector, double *dblVector){
	intVector[0]=varNumber;
	intVector[1]=(direction=='L')*1+(direction=='U')*2+(direction=='E')*3;
	intVector[2]=stage;
	intVector[3]=scenario;
	dblVector[0]=bound;
	
}

 void BBSMPSBranchingInfo::getSerializationSize(int &intVectorSize, int &dblVectorSize){
	intVectorSize=4;//var number, scenario, direction, stage
	dblVectorSize=1;//the bound
}
