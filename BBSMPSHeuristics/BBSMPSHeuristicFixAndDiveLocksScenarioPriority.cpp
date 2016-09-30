#include "BBSMPSHeuristicFixAndDiveLocksScenarioPriority.hpp"

using namespace std;






bool BBSMPSHeuristicFixAndDiveLocksScenarioPriority::runHeuristic(BBSMPSNode* node, denseBAVector &LPRelaxationSolution){
	
	double startTimeStamp = MPI_Wtime();
	int mype=BBSMPSSolver::instance()->getMype();
	if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Performing the Fix and Dive Locks heuristic with scenario priority.";
	timesCalled++;


	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	
	const BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getOriginalBADimensionsSlacks();
    BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
    const denseBAVector varObjectives = rootSolver.getVarObjective();
    

    const BAFlagVector<constraintType> varTypes = rootSolver.getVariableTypes();
    
    int firstStageVars=input.nFirstStageVars();
    int firstStageRows=input.nFirstStageCons();
   
   	int MAX_ITERS=firstStageVars;
   	for (int scen = 0; scen < input.nScenarios(); scen++)
	{
		MAX_ITERS+=input.nSecondStageVars(scen);
	}


	denseBAVector upLocks;
    denseBAVector downLocks;
    upLocks.allocate(dimsSlacks, ctx, PrimalVector);
	downLocks.allocate(dimsSlacks, ctx, PrimalVector);
   	upLocks.clear();
   	downLocks.clear();
    generateLocks(upLocks,downLocks);
	denseBAVector variableObjectives = rootSolver.getVarObjective();
    denseBAVector lb(BBSMPSSolver::instance()->getOriginalLB());
   	denseBAVector ub(BBSMPSSolver::instance()->getOriginalUB());
	node->getAllBranchingInformation(lb,ub);
	denseBAVector auxSolution(LPRelaxationSolution);
	//Get all variables
	//BAFlagVector<variableState> ps(BBSMPSSolver::instance()->getOriginalWarmStart());
	//node->reconstructWarmStartState(ps);
	//rootSolver.setStates(ps);
    double bestLockIndex=-1;
    int bestLock=-1;
    double bestVarObj=COIN_DBL_MAX;
	double objCont=COIN_DBL_MAX;//We are trying to minimize
	
	//Order by most to least fractional
    for (int i=0; i< input.nFirstStageVars(); i++){
       	
    	if (((upLocks.getFirstStageVec()[i]> bestLock)|| (upLocks.getFirstStageVec()[i]== bestLock && objCont> objContribution(auxSolution.getFirstStageVec()[i],varObjectives.getFirstStageVec()[i],0))) && input.isFirstStageColInteger(i) && !isIntFeas(auxSolution.getFirstStageVec()[i], intTol)){
    		bestLockIndex=i;
    		bestLock=upLocks.getFirstStageVec()[i];
    		bestVarObj=variableObjectives.getFirstStageVec()[i];
    		objCont= objContribution(auxSolution.getFirstStageVec()[i],varObjectives.getFirstStageVec()[i],0);
    		
    	}
    	if (((downLocks.getFirstStageVec()[i]> bestLock)|| (downLocks.getFirstStageVec()[i]== bestLock && objCont> objContribution(auxSolution.getFirstStageVec()[i],varObjectives.getFirstStageVec()[i],1))) && input.isFirstStageColInteger(i) && !isIntFeas(auxSolution.getFirstStageVec()[i], intTol)){
    		bestLockIndex=i;
    		bestLock=downLocks.getFirstStageVec()[i];
    		bestVarObj=variableObjectives.getFirstStageVec()[i];
    		objCont= objContribution(auxSolution.getFirstStageVec()[i],varObjectives.getFirstStageVec()[i],1);
    	
    	}
    }
    bool allInteger=(bestLockIndex==-1);
	int iter=0;
	while(!allInteger && iter<MAX_ITERS){
		
		if (bestLock==upLocks.getFirstStageVec()[bestLockIndex] && bestLock==downLocks.getFirstStageVec()[bestLockIndex]){

			lb.getFirstStageVec()[bestLockIndex]=roundToNearestInteger(auxSolution.getFirstStageVec()[bestLockIndex]);
			ub.getFirstStageVec()[bestLockIndex]=roundToNearestInteger(auxSolution.getFirstStageVec()[bestLockIndex]);
		
		}
		else if (bestLock==upLocks.getFirstStageVec()[bestLockIndex]){

			lb.getFirstStageVec()[bestLockIndex]=floor(auxSolution.getFirstStageVec()[bestLockIndex]);
			ub.getFirstStageVec()[bestLockIndex]=floor(auxSolution.getFirstStageVec()[bestLockIndex]);	
		}
		else{

			lb.getFirstStageVec()[bestLockIndex]=ceil(auxSolution.getFirstStageVec()[bestLockIndex]);
			ub.getFirstStageVec()[bestLockIndex]=ceil(auxSolution.getFirstStageVec()[bestLockIndex]);
		}
		
		rootSolver.setLB(lb);
		rootSolver.setUB(ub);

		rootSolver.commitStates();	

		rootSolver.go();

		solverState lpStatus = rootSolver.getStatus();
		bool otherThanOptimal = (Optimal != lpStatus); 
		/*if (otherThanOptimal){
			cumulativeTime+=(MPI_Wtime()-startTimeStamp);
			return false;
		} */
		auxSolution=rootSolver.getPrimalSolution();

		bestLockIndex=-1;
  		bestLock=-1;
		bestVarObj=COIN_DBL_MAX;
	    objCont=COIN_DBL_MAX;//We are trying to minimize
			
		//Order by most to least fractional
	    for (int i=0; i< input.nFirstStageVars(); i++){
	    	if (((upLocks.getFirstStageVec()[i]> bestLock)|| (upLocks.getFirstStageVec()[i]== bestLock && objCont> objContribution(auxSolution.getFirstStageVec()[i],varObjectives.getFirstStageVec()[i],0))) && input.isFirstStageColInteger(i) && !isIntFeas(auxSolution.getFirstStageVec()[i], intTol)){
	    		bestLockIndex=i;
	    		bestLock=upLocks.getFirstStageVec()[i];
	    		bestVarObj=variableObjectives.getFirstStageVec()[i];
	    		objCont= objContribution(auxSolution.getFirstStageVec()[i],varObjectives.getFirstStageVec()[i],0);
	    		
	    	}
	    	if (((downLocks.getFirstStageVec()[i]> bestLock)|| (downLocks.getFirstStageVec()[i]== bestLock && objCont> objContribution(auxSolution.getFirstStageVec()[i],varObjectives.getFirstStageVec()[i],1))) && input.isFirstStageColInteger(i) && !isIntFeas(auxSolution.getFirstStageVec()[i], intTol)){
	    		bestLockIndex=i;
	    		bestLock=downLocks.getFirstStageVec()[i];
	    		bestVarObj=variableObjectives.getFirstStageVec()[i];
	    		objCont= objContribution(auxSolution.getFirstStageVec()[i],varObjectives.getFirstStageVec()[i],1);
	    	
	    	}
	    }
	    iter++;
	    allInteger=(bestLockIndex==-1);

	}
 
 	for (int i=0; i< input.nFirstStageVars(); i++){
 		lb.getFirstStageVec()[i]=auxSolution.getFirstStageVec()[i];
		ub.getFirstStageVec()[i]=auxSolution.getFirstStageVec()[i];
	}


	int lowestFracValue=9999999;
	int lowestScenario=-1;
	for (int scen = 0; scen < input.nScenarios(); scen++)
	{
		if(ctx.assignedScenario(scen)) {
			int numElems=0;
			for (int i = 0; i < input.nSecondStageVars(scen); i++)
			{
				if (input.isSecondStageColInteger(scen,i) && !isIntFeas(auxSolution.getSecondStageVec(scen)[i], intTol))numElems++;
			}
			if (numElems>0&&(lowestScenario==-1 || numElems<lowestFracValue)){
				lowestScenario=scen;
				lowestFracValue=numElems;
			}
		}
	}

	intIntComp localComp;
    intIntComp globalComp;

    localComp.number=lowestFracValue;
    localComp.rank=lowestScenario;
    MPI_Allreduce (&localComp, &globalComp, 1, MPI_2INT, MPI_MINLOC, ctx.comm());
		    		


	int bestLockScen=-1;
	bestLockIndex=-1;
	bestLock=-1;
	bestVarObj=COIN_DBL_MAX;
	 objCont=COIN_DBL_MAX;//We are trying to minimize
	
	int scen =globalComp.rank;
	if(scen!=-1&&ctx.assignedScenario(scen)) {
		for (int i = 0; i < input.nSecondStageVars(scen); i++)
		{
			if (((upLocks.getSecondStageVec(scen)[i] >bestLock) || (upLocks.getSecondStageVec(scen)[i] ==bestLock && objCont> objContribution(auxSolution.getSecondStageVec(scen)[i],varObjectives.getSecondStageVec(scen)[i],0))) && input.isSecondStageColInteger(scen,i) && !isIntFeas(auxSolution.getSecondStageVec(scen)[i], intTol)){
	    		bestLockIndex=i;
	    		bestLockScen=scen;
	    		bestLock=upLocks.getSecondStageVec(scen)[i];
	    		objCont= objContribution(auxSolution.getSecondStageVec(scen)[i],varObjectives.getSecondStageVec(scen)[i],0);
	    	}
	    	if (((downLocks.getSecondStageVec(scen)[i] >bestLock)||(downLocks.getSecondStageVec(scen)[i]==bestLock && objCont> objContribution(auxSolution.getSecondStageVec(scen)[i],varObjectives.getSecondStageVec(scen)[i],1))) && input.isSecondStageColInteger(scen,i) && !isIntFeas(auxSolution.getSecondStageVec(scen)[i], intTol)){
	    		bestLockIndex=i;
	    		bestLockScen=scen;
	    		bestLock=downLocks.getSecondStageVec(scen)[i];
	    		objCont= objContribution(auxSolution.getSecondStageVec(scen)[i],varObjectives.getSecondStageVec(scen)[i],1);
	    	}
	    }
	}
	
	
	
	int iteration=0;
	iter =0;
	while (scen!=-1 && iteration<MAX_ITERS){
		iteration++;
		if (ctx.assignedScenario(scen)){
			if (bestLock==upLocks.getSecondStageVec(bestLockScen)[bestLockIndex] && bestLock==downLocks.getSecondStageVec(bestLockScen)[bestLockIndex]){
				lb.getSecondStageVec(bestLockScen)[bestLockIndex]=roundToNearestInteger(auxSolution.getSecondStageVec(bestLockScen)[bestLockIndex]);
				ub.getSecondStageVec(bestLockScen)[bestLockIndex]=roundToNearestInteger(auxSolution.getSecondStageVec(bestLockScen)[bestLockIndex]);
		
			}

	
			if (bestLock==upLocks.getSecondStageVec(bestLockScen)[bestLockIndex] ){
	
				lb.getSecondStageVec(bestLockScen)[bestLockIndex]=floor(auxSolution.getSecondStageVec(bestLockScen)[bestLockIndex]);
				ub.getSecondStageVec(bestLockScen)[bestLockIndex]=floor(auxSolution.getSecondStageVec(bestLockScen)[bestLockIndex]);
			
			}
			else{
	
				lb.getSecondStageVec(bestLockScen)[bestLockIndex]=ceil(auxSolution.getSecondStageVec(bestLockScen)[bestLockIndex]);
				ub.getSecondStageVec(bestLockScen)[bestLockIndex]=ceil(auxSolution.getSecondStageVec(bestLockScen)[bestLockIndex]);
			}
		
			
		}
		rootSolver.setLB(lb);
		rootSolver.setUB(ub);

		rootSolver.commitStates();	

		rootSolver.go();

		solverState lpStatus = rootSolver.getStatus();
		bool otherThanOptimal = (Optimal != lpStatus); 
		int anyOtherThanOptimal=0;
	   //  MPI_Allreduce(&otherThanOptimal, &anyOtherThanOptimal, 1, MPI_INT,  MPI_MAX, ctx.comm());

		/*if (anyOtherThanOptimal!=0){
			cumulativeTime+=(MPI_Wtime()-startTimeStamp);
			return false;
		} */
		auxSolution=rootSolver.getPrimalSolution();

		lowestFracValue=9999999;
		lowestScenario=-1;
		for (int scenn = 0; scenn < input.nScenarios(); scenn++)
		{
			if(ctx.assignedScenario(scenn)) {
				int numElems=0;
				for (int i = 0; i < input.nSecondStageVars(scenn); i++)
				{
					if (input.isSecondStageColInteger(scenn,i) && !isIntFeas(auxSolution.getSecondStageVec(scenn)[i], intTol))numElems++;
				}
				if (numElems>0&&(lowestScenario==-1 || numElems<lowestFracValue)){
					lowestScenario=scenn;
					lowestFracValue=numElems;
				}
			}
		}

		localComp;
	    globalComp;

	    localComp.number=lowestFracValue;
	    localComp.rank=lowestScenario;
	    MPI_Allreduce (&localComp, &globalComp, 1, MPI_2INT, MPI_MINLOC, ctx.comm());


		bestLockScen=-1;
		bestLockIndex=-1;
		bestLock=-1;
	
		bestVarObj=COIN_DBL_MAX;
		objCont=COIN_DBL_MAX;//We are trying to minimize
		
		scen =globalComp.rank;
		if(scen!=-1&&ctx.assignedScenario(scen)) {
			for (int i = 0; i < input.nSecondStageVars(scen); i++)
			{
				if (((upLocks.getSecondStageVec(scen)[i] >bestLock) || (upLocks.getSecondStageVec(scen)[i] ==bestLock && objCont> objContribution(auxSolution.getSecondStageVec(scen)[i],varObjectives.getSecondStageVec(scen)[i],0))) && input.isSecondStageColInteger(scen,i) && !isIntFeas(auxSolution.getSecondStageVec(scen)[i], intTol)){
		    		bestLockIndex=i;
		    		bestLockScen=scen;
		    		bestLock=upLocks.getSecondStageVec(scen)[i];
		    		objCont= objContribution(auxSolution.getSecondStageVec(scen)[i],varObjectives.getSecondStageVec(scen)[i],0);
		    	}
		    	if (((downLocks.getSecondStageVec(scen)[i] >bestLock)||(downLocks.getSecondStageVec(scen)[i]==bestLock && objCont> objContribution(auxSolution.getSecondStageVec(scen)[i],varObjectives.getSecondStageVec(scen)[i],1))) && input.isSecondStageColInteger(scen,i) && !isIntFeas(auxSolution.getSecondStageVec(scen)[i], intTol)){
		    		bestLockIndex=i;
		    		bestLockScen=scen;
		    		bestLock=downLocks.getSecondStageVec(scen)[i];
		    		objCont= objContribution(auxSolution.getSecondStageVec(scen)[i],varObjectives.getSecondStageVec(scen)[i],1);
		    	}
		    }
		}
	    iter++;
  
	}

   
	solverState lpStatus = rootSolver.getStatus();
	bool otherThanOptimal = (Optimal != lpStatus); 
	double objUB=COIN_DBL_MAX;
	if (BBSMPSSolver::instance()->getSolPoolSize()>0)objUB=BBSMPSSolver::instance()->getSoln(0).getObjValue();

	//return if success

	
	if(!otherThanOptimal && rootSolver.getObjective()<objUB){

		denseBAVector solVector=rootSolver.getPrimalSolution();

		if (isLPIntFeas(solVector)){
					BBSMPSSolution sol(solVector,rootSolver.getObjective());
					sol.setTimeOfDiscovery(BBSMPSSolver::instance()->getWallTime());
					BBSMPSSolver::instance()->addSolutionToPool(sol);
		
		}
		
		 
	}
	//return if success
	bool success= (!otherThanOptimal && rootSolver.getObjective()<objUB);
	timesSuccessful+=(success);
	
	cumulativeTime+=(MPI_Wtime()-startTimeStamp);
	return false;
	return success;

}

bool BBSMPSHeuristicFixAndDiveLocksScenarioPriority::shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution){
	
	if (node->getNodeDepth()<5)return true;
	//else depth=25;

	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	BAContext &ctx= BBSMPSSolver::instance()->getBAContext();
	
	int numberOfFractionalVariables=0;
	int nIntVars=0;
	for (int col = 0; col < input.nFirstStageVars(); col++)
	{	
		if(input.isFirstStageColInteger(col)){
			if (!isIntFeas(LPRelaxationSolution.getFirstStageVec()[col],intTol)) numberOfFractionalVariables++;
			nIntVars++;
		}
		
	}
	
	int numberOfFractionalVariables2=0;
	int nIntVars2=0;

	
	for (int scen = 0; scen < input.nScenarios(); scen++)
	{
		if(ctx.assignedScenario(scen)) {
			for (int col = 0; col < input.nSecondStageVars(scen); col++)
			{

				if(input.isSecondStageColInteger(scen,col)){
					if (!isIntFeas(LPRelaxationSolution.getSecondStageVec(scen)[col],intTol)) numberOfFractionalVariables2++;
					nIntVars2++;
					
				}
			}
		}
	}

	int totalCount2;
	int errorFlag = MPI_Allreduce(&numberOfFractionalVariables2,
		&totalCount2,
		1,
		MPI_INT, 
		MPI_SUM,
		ctx.comm());

	int totalIntVars2;
	errorFlag = MPI_Allreduce(&nIntVars2,
		&totalIntVars2,
		1,
		MPI_INT, 
		MPI_SUM,
		ctx.comm());


	nIntVars+=totalIntVars2;
	numberOfFractionalVariables=+totalCount2;
	if (nIntVars==0)return false;
	return ((numberOfFractionalVariables*100/nIntVars)<400 );

}