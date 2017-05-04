#include "BBSMPSHeuristic.hpp"

using namespace std;

BBSMPSHeuristic::BBSMPSHeuristic(int _offset, int _depth, const char *_name){
	offset=_offset;
	depth=_depth;
	timesCalled=0;
	name=_name;
	timesSuccessful=0;
	cumulativeTime=0;
}
BBSMPSHeuristic::~BBSMPSHeuristic(){};

bool BBSMPSHeuristic::checkPeriodicity(BBSMPSNode* node){

	int nodeDepth=node->getNodeDepth();
	BBSMPS_ALG_LOG_SEV(debug) << "Node depth, heur offset, heur depth: " << nodeDepth << ", " << offset << ", " << depth;
	assert(nodeDepth>=0);
	if (offset>nodeDepth) return false;
	nodeDepth-=offset;
	return (nodeDepth%depth==0);
}


void BBSMPSHeuristic::printStatistics(){
	BBSMPS_ALG_LOG_SEV(warning)<<"Heuristic:"<<name<<":Times Called:"<<timesCalled<<":Times successful:"<<timesSuccessful<<":Execution Time:"<<cumulativeTime;
}
double BBSMPSHeuristic::objContribution(double valueToRound, double objCoefficient, int roundingDirection){

	double obj=fabs(objCoefficient);
	if (roundingDirection==0){//We are rounding down

		double differential= floor(valueToRound) - valueToRound;
		return obj*differential;
	}

		double differential= ceil(valueToRound) - valueToRound;
		return obj*differential;

}

void BBSMPSHeuristic::generateLocks(denseBAVector & upLocks, denseBAVector &downLocks){



	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();

	BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();

   	int firstStageVars=input.nFirstStageVars();
    int firstStageRows=input.nFirstStageCons();
   	const BAFlagVector<constraintType> varTypes = rootSolver.getVariableTypes();
   	for (int scen = 0; scen < input.nScenarios(); scen++)
	{
		if(ctx.assignedScenario(scen)) {
			for (int c = 0; c < input.nSecondStageCons(scen); c++)
			{

				const CoinShallowPackedVector row=rootSolver.retrieveTRow(c,scen);
				int nElems=row.getNumElements();
				const int*indices=row.getIndices();
				const double *elems=row.getElements();


		    	const CoinShallowPackedVector row2=rootSolver.retrieveWRow(c,scen);
		    	int nElems2=row2.getNumElements();
				const int* indices2=row2.getIndices();
				const double *elems2=row2.getElements();

				if (varTypes.getSecondStageVec(scen)[input.nSecondStageVars(scen)+c]==LB){
					for (int el=0; el<nElems; el++){

			    		if (elems[el]<0)upLocks.getFirstStageVec()[indices[el]]++;
			    		else downLocks.getFirstStageVec()[indices[el]]++;

			    	}

			    	for (int el=0; el<nElems2; el++){
			    		if (elems2[el]<0)upLocks.getSecondStageVec(scen)[indices2[el]]++;
			    		else downLocks.getSecondStageVec(scen)[indices2[el]]++;

			    	}
			    }
			    else if (varTypes.getSecondStageVec(scen)[input.nSecondStageVars(scen)+c]==UB){
			    	for (int el=0; el<nElems; el++){

			    		if (elems[el]>0)upLocks.getFirstStageVec()[indices[el]]++;
			    		else downLocks.getFirstStageVec()[indices[el]]++;

			    	}

			    	for (int el=0; el<nElems2; el++){
			    		if (elems2[el]>0)upLocks.getSecondStageVec(scen)[indices2[el]]++;
			    		else downLocks.getSecondStageVec(scen)[indices2[el]]++;

			    	}
			    }
			    else if(varTypes.getSecondStageVec(scen)[input.nSecondStageVars(scen)+c]==Fixed){
			    	for (int el=0; el<nElems; el++){

			    		upLocks.getFirstStageVec()[indices[el]]++;
			    	    downLocks.getFirstStageVec()[indices[el]]++;

			    	}

			    	for (int el=0; el<nElems2; el++){
			    		upLocks.getSecondStageVec(scen)[indices2[el]]++;
			    		downLocks.getSecondStageVec(scen)[indices2[el]]++;

			    	}
			    }


		    }

		}
    }

    //At this point each vector has its own count of first stage locks. Let's reduce
    double *upLock1stStagePtr = upLocks.getFirstStageVec().getPointer();
    double *downLock1stStagePtr = downLocks.getFirstStageVec().getPointer();
    MPI_Allreduce(MPI_IN_PLACE,upLock1stStagePtr,upLocks.getFirstStageVec().length(),MPI_DOUBLE,MPI_SUM,ctx.comm());
    MPI_Allreduce(MPI_IN_PLACE,downLock1stStagePtr,downLocks.getFirstStageVec().length(),MPI_DOUBLE,MPI_SUM,ctx.comm());

    for (int c=0; c< firstStageRows ; c++){

    	const CoinShallowPackedVector row=rootSolver.retrieveARow(c);
    	int nElems=row.getNumElements();
		const int*indices=row.getIndices();
		const double *elems=row.getElements();

		if (varTypes.getFirstStageVec()[input.nFirstStageVars()+c]==LB){
			for (int el=0; el<nElems; el++){
	    		if (elems[el]<0)upLocks.getFirstStageVec()[indices[el]]++;
	    		else downLocks.getFirstStageVec()[indices[el]]++;

	    	}
	    }
	    else if (varTypes.getFirstStageVec()[input.nFirstStageVars()+c]==UB){
	    	for (int el=0; el<nElems; el++){
    			if (elems[el]>0)upLocks.getFirstStageVec()[indices[el]]++;
    			else downLocks.getFirstStageVec()[indices[el]]++;
    		}
	    }
	    else if(varTypes.getFirstStageVec()[input.nFirstStageVars()+c]==Fixed){
	    	for (int el=0; el<nElems; el++){
				upLocks.getFirstStageVec()[indices[el]]++;
	    	    downLocks.getFirstStageVec()[indices[el]]++;

	    	}
    	}

    }
}
