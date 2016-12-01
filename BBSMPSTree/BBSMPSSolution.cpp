#include "BBSMPSSolution.hpp"

using namespace std;
	int BBSMPSSolution::solCounter=0;

	bool BBSMPSSolution::operator==(const BBSMPSSolution &other) const {
	    return (objValue==other.getObjValue());
	  }

	BBSMPSSolution::BBSMPSSolution(const denseBAVector &_solutionVector, double _objValue,double _timeOfDiscovery):
	objValue(_objValue),timeOfDiscovery(_timeOfDiscovery),solutionVector(_solutionVector){
		solNumber=(++solCounter);
	}



	BBSMPSSolution::~BBSMPSSolution(){}

	void BBSMPSSolution::setSolutionVector(const denseBAVector &_solutionVector){
		solutionVector.copyFrom(_solutionVector);
	}
	void BBSMPSSolution::getSolutionVector(denseBAVector &_solutionVector)const{
		
		_solutionVector = solutionVector;
	}
	void BBSMPSSolution::setObjValue(double _objValue){
		objValue=_objValue;
	}
	double BBSMPSSolution::getObjValue()const{
		return objValue;
	}

	double BBSMPSSolution::getTimeOfDiscovery() const{
		return timeOfDiscovery;
	}

	void BBSMPSSolution::setTimeOfDiscovery(const double _timeOfDiscovery){
	  timeOfDiscovery=_timeOfDiscovery;
	}

	int BBSMPSSolution::getSolNumber()const{
		return solNumber;
	}

	int BBSMPSSolution::getSerializationSize() const{
		int counter=1+1+1;
		const std::vector<int>& scens =solutionVector.localScenarios();
		for (int i=0; i< scens.size();i++)counter=counter+solutionVector.getVec(scens[i]).length();
		counter=counter+solutionVector.getVec(-1).length();
		return counter;
	}

	void BBSMPSSolution::serialize(std::vector<double> &vec)const{
		vec[0]=solNumber;
		vec[1]=objValue;
		vec[2]=timeOfDiscovery;
		int ptr=3;
		const std::vector<int>& scens =solutionVector.localScenarios();
		for (int i=0; i< solutionVector.getVec(-1).length(); i++){
			vec[ptr]=solutionVector.getVec(-1)[i];
			ptr++;
		}
		for (int j=0; j< scens.size();j++){
			for (int i=0; i< solutionVector.getVec(scens[j]).length(); i++){
				vec[ptr]=solutionVector.getVec(scens[j])[i];
				ptr++;
			}
		}
	}
	void BBSMPSSolution::initializeFromVector(std::vector<double> &vec,const BAContext &ctx,const BADimensions &dims){
		solNumber=vec[0];
		objValue=vec[1];
		timeOfDiscovery=vec[2];
		solutionVector.allocate(dims, ctx, PrimalVector); 
		int ptr=3;
		const std::vector<int>& scens =solutionVector.localScenarios();
		for (int i=0; i< solutionVector.getVec(-1).length(); i++){
			solutionVector.getVec(-1)[i]=vec[ptr];
			ptr++;
		}
		for (int j=0; j< scens.size();j++){
			for (int i=0; i< solutionVector.getVec(scens[j]).length(); i++){
				solutionVector.getVec(scens[j])[i]=vec[ptr];
				ptr++;
			}
		}
	}