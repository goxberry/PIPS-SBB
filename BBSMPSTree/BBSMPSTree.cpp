#include "BBSMPSTree.hpp"

using namespace std;

double intTol=1e-6;
BBSMPSNode *BBSMPSTree::rootNode =NULL;
bool boundVectorComparator (pair <double, pair<int, int> > i, pair <double, pair<int, int> > j) {
	if (abs(i.first-j.first)<intTol){
		if (i.second.first==j.second.first)return i.second.second<j.second.second;
		return i.second.first<j.second.first;
	}
 return (i.first<j.first); 
 }

// Outputs solver status:
void outputLPStatus(solverState lpStatus) {
	
	string status;
	switch(lpStatus) {
		case Uninitialized:
		status = "Uninitialized";
		break;
		case LoadedFromFile:
		status = "LoadedFromFile";
		break;
		case Initialized:
		status = "Initialized";
		break;
		case PrimalFeasible:
		status = "PrimalFeasible";
		break;
		case DualFeasible:
		status = "DualFeasible";
		break;
		case Optimal:
		status = "Optimal";
		break;
		case ProvenUnbounded:
		status = "ProvenUnbounded";
		break;
		case ProvenInfeasible:
		status = "ProvenInfeasible";
		break;
		case Stopped:
		status = "Stopped";
		break;
	}
	BBSMPS_ALG_LOG_SEV(info) <<"PIPS-S has returned "<< status;
}



// Overload the "less than" operator so that priority_queue can use it
// for heapifying comparisons. Make BBSMPSNode a templated
// class when adding additional branching heuristics?
bool operator< (const BBSMPSNode& left,
	const BBSMPSNode&right)
{
	return left.getParentObjective() < right.getParentObjective();// left.parentObj < right.parentObj;
}

bool operator> (const BBSMPSNode& left,
	const BBSMPSNode& right)
{
	return  left.getParentObjective() > right.getParentObjective();//left.parentObj > right.parentObj;
}


BBSMPSTree::BBSMPSTree(const SMPSInput& smps,int nSolvers): 
objUB(COIN_DBL_MAX),
objLB(-COIN_DBL_MAX),
optGapTol(1e-6),
lpPrimalTol(1e-6),
lpDualTol(1e-6),
commTol(0.05),
compTol(lpPrimalTol),
status(LoadedFromFile),
nodesel(BestBound),
tiLim(COIN_INT_MAX),
nodeLim(COIN_INT_MAX),
solsDiscoveredLimit(COIN_INT_MAX),
solsDiscoveredInit(0),
verbosityActivated(true),
cuttingPlanesManager(5),
initializationTime(0)
{
	double timeStart=MPI_Wtime();

	BBSMPSSolver::initialize(smps,nSolvers); 
	double timeStampPreProc=MPI_Wtime();
	PreProcessingTime=timeStampPreProc-timeStart;
   // if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Calling B&B tree constructor.";

    /* Initialize branch-and-bound tree/heap */
    // Get {lower, upper} bounds on decision variables, lower bound on objective function
    // value from parent LP, initialize a node, and insert onto heap to start.
    assert (heap.empty()); // heap should be empty to start

	
    
    PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
    const BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
    BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    int mype=BBSMPSSolver::instance()->getMype();

    rootSolver.setPrimalTolerance(lpPrimalTol);
    rootSolver.setDualTolerance(lpDualTol);

    // TODO: Replace with real presolve.
    // For now, "cheat" by solving root LP before populating root node.
    // A warm start of the root node means that the root node solve
    // inside the B&B tree does not cost very much -- just the PIPS-S overhead.
    // However, this step is necessary in order to properly allocate primal
    // variables (due to slacks) and to determine which variables are basic.
    // This step is also currently required to instantiate lower & upper
    // bounds. In theory, the bounds could be obtained from the SMPS file.
    // In practice, getting the bounds from the SMPS file is cumbersome,
    // because first stage bounds and second stage bounds must be
    // queried separately and are returned as std::vector<double>s.
    // In addition, distributed data structures dictate some care in
    // how the assignments are performed: there must be checks to
    // ensure that the data to be assigned is owned by the "right" process.

    rootSolver.go();

   // BBSMPSCuttingPlane bcp;
		//for (int i=0; i<=bbIterationCounter; i++)	
	//bcp.applyCuttingPlane();

    LPRelaxationTime=MPI_Wtime()-timeStampPreProc;
    // Get lower & upper bounds on decision variables in LP.
    if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Getting bounds for root node from presolve.";
    denseBAVector lb(rootSolver.getLB()), ub(rootSolver.getUB());
       // Allocate current best primal solution; normally this primal solution
    // is for the upper bound, but here, we have only the solution to an
    // LP relaxation, which may not be primal feasible. We don't check
    // primal/integer feasibility here.
    //pwdif (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Allocating primal solution.";
    ubPrimalSolution.allocate(dimsSlacks, ctx, PrimalVector);
    //if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Getting primal solution";
    ubPrimalSolution.copyFrom(rootSolver.getPrimalSolution());
    //if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "MIP Primal solution updated"; 
    BBSMPSSolver::instance()->setLPRelaxation(ubPrimalSolution);
    // Update state of primal variables + slacks; slacks must be included
    // because these are used in a reformulation of the problem to standard
    // form: Ax + s = b, s >= 0.

    BAFlagVector<variableState> states(dimsSlacks, ctx, PrimalVector);
    
    rootSolver.getStates(states);
    BBSMPSSolver::instance()->setOriginalWarmStart(states);
    // Update global lower bound; really need to check feasibility, etc.
    // here, but I'm going to move this code to the B&B tree.
    double lpObj = rootSolver.getObjective();
    LPRelaxationValue=lpObj;
    if ((lpObj - compTol) >= objLB) objLB = lpObj;

    std::vector< std::pair < BAIndex, variableState > > emptyStates;

int BBSMPSMyPe=BBSMPSSolver::instance()->getSBBMype();
    BBSMPSNode *rootNode= new BBSMPSNode(lpObj,emptyStates,BBSMPSSolver::instance()->getSBBMype());

    //if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Pushing root node onto B&B tree.";
    heap.insert(rootNode);
    //if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Exiting B&B constructor.";

   //BBSMPSMaxFracBranchingRule *mfbr= new BBSMPSMaxFracBranchingRule(10);
   //branchingRuleManager.addBranchingRule(mfbr);


    BBSMPSPseudoCostBranchingRule *mfbr2= new BBSMPSPseudoCostBranchingRule(100);
  	  branchingRuleManager.addBranchingRule(mfbr2);

    bbIterationCounter=0;
    nodesFathomed=0;
    nodesBecameInteger=0;
    totalBbIterationsDone=0;
    BBSMPSTree::rootNode=rootNode;
    inTermination=false;
    request0=new MPI_Request();
    request1=new MPI_Request();
    request2=new MPI_Request();
    request3=new MPI_Request();
    request4=new MPI_Request();
    *request0 = MPI_REQUEST_NULL;
    *request1 = MPI_REQUEST_NULL;
    *request2 = MPI_REQUEST_NULL;
    *request3 = MPI_REQUEST_NULL;
    *request4 = MPI_REQUEST_NULL;
    communicationPending=false;

    communicationActivated=true;
    iterationsBetweenCommunication=MIN_COMM_ITERS;
}



BBSMPSTree::BBSMPSTree(BBSMPSNode *node, double lb, double ub): 
objUB(ub),
objLB(lb),
optGapTol(1e-6),
lpPrimalTol(1e-6),
lpDualTol(1e-6),
commTol(0.05),
compTol(lpPrimalTol),
status(LoadedFromFile),
nodesel(BestBound),
tiLim(COIN_INT_MAX),
nodeLim(COIN_INT_MAX),
solsDiscoveredLimit(COIN_INT_MAX),
solsDiscoveredInit(0),
verbosityActivated(true),
cuttingPlanesManager(2),
initializationTime(0)
{

	//This initialization assumes the solver class has already been intialized
	assert(BBSMPSSolver::isInitialized());

    //if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Calling B&B tree constructor.";

    /* Initialize branch-and-bound tree/heap */
    // Get {lower, upper} bounds on decision variables, lower bound on objective function
    // value from parent LP, initialize a node, and insert onto heap to start.
    assert (heap.empty()); // heap should be empty to start
    
    PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
   
    double lpObj = -INFINITY;
    if ((lpObj - compTol) >= objLB) objLB = lpObj;

    
    //if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Pushing root node onto B&B tree.";
    heap.insert(node);
    //if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Exiting B&B constructor.";

    BBSMPSMaxFracBranchingRule *mfbr= new BBSMPSMaxFracBranchingRule(10);
    branchingRuleManager.addBranchingRule(mfbr);


    bbIterationCounter=0;
     nodesFathomed=0;
   nodesBecameInteger=0;
   totalBbIterationsDone=0;
   removeCuts();
	//node->getAllCuttingUids(currentlyAppliedPlanes);
	inTermination=false;
	request0=new MPI_Request();
    request1=new MPI_Request();
    request2=new MPI_Request();
    request3=new MPI_Request();
    request4=new MPI_Request();
    *request0 = MPI_REQUEST_NULL;
    *request1 = MPI_REQUEST_NULL;
    *request2 = MPI_REQUEST_NULL;
    *request3 = MPI_REQUEST_NULL;
    *request4 = MPI_REQUEST_NULL;
    communicationPending=false;
    communicationActivated=false;
    iterationsBetweenCommunication=MIN_COMM_ITERS;
}

BBSMPSTree::~BBSMPSTree(){
	//Dismantling nodes
	while (!heap.empty()){
		BBSMPSNode *currentNode_ptr=*(heap.begin());
		currentNode_ptr->eliminate();
		
		heap.erase(heap.begin());
	}

	delete request0;
	delete request1;
	delete request2;
	delete request3;
	delete request4;

	//Dismantling Branching Rules
  	branchingRuleManager.freeResources();

	//Dismantling Heuristics
	heuristicsManager.freeResources();

	//Dismantling Cutting Planes
	cuttingPlanesManager.freeResources();
}



void BBSMPSTree::generateIncrementalWarmState(BBSMPSNode* node, const BAFlagVector<variableState> & originalState, const BAFlagVector<variableState> &currentState){

	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	BAContext &ctx=BBSMPSSolver::instance()->getBAContext();

	std::vector< std::pair < BAIndex, variableState > > changes;
	

	const denseFlagVector<variableState> &stageVec = currentState.getFirstStageVec();
	const denseFlagVector<variableState> &stageVec2 = originalState.getFirstStageVec();
	for (int j = 0; j < stageVec.length(); j++) {
		if (stageVec[j]!=stageVec2[j]){
			BAIndex aux;
			aux.scen=-1;
			aux.idx=j;
			changes.push_back(std::pair < BAIndex, variableState > (aux,stageVec[j]));

		}
	}


	for (int scen = 0; scen < input.nScenarios(); scen++) {
		if(ctx.assignedScenario(scen)) {
			const denseFlagVector<variableState> &stageVec = currentState.getSecondStageVec(scen);
			const denseFlagVector<variableState> &stageVec2 = originalState.getSecondStageVec(scen);
			for (int j = 0; j < stageVec.length(); j++) {
				if (stageVec[j]!=stageVec2[j]){
					BAIndex aux;
					aux.scen=scen;
					aux.idx=j;
					changes.push_back(std::pair < BAIndex, variableState > (aux,stageVec[j]));

				}
			}
		}
	}

	

	node->setIncrementalWarmStartState(changes);



}
void BBSMPSTree::checkSequentialTerminationConditions(){

	int mype=BBSMPSSolver::instance()->getMype();

	if ((BBSMPSSolver::instance()->getWallTime())>tiLim){
		if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Time Limit reached.";

		inTermination=true;
		
		return ;
	}
	if (BBSMPSSolver::instance()->getSolPoolSize()-solsDiscoveredInit>=solsDiscoveredLimit){
		if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Solution Limit reached.";
		inTermination=true;
		return ;
	}

	if (bbIterationCounter>nodeLim){
		if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Node Limit reached.";
		inTermination=true;
		
		return;
	}

	if (heap.size()==0){
		if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Heap is empty.";
		inTermination=true;
	
		return ;
	}
}



void BBSMPSTree::removeCuts(){


	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
	PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
	const BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
	const BADimensionsSlacks &originalDimsSlacks= BBSMPSSolver::instance()->getOriginalBADimensionsSlacks();
	
	denseBAVector lb(BBSMPSSolver::instance()->getOriginalLB());
	denseBAVector ub(BBSMPSSolver::instance()->getOriginalUB());
	rootSolver.setLB(lb);
	rootSolver.setUB(ub);

	bool modelChanged=false;

	for (int scen = 0; scen < input.nScenarios(); scen++) {
		if(ctx.assignedScenario(scen)) {
			if (originalDimsSlacks.inner.numSecondStageCons(scen) < dimsSlacks.inner.numSecondStageCons(scen)){
				rootSolver.deleteLastSecondStageConsecutiveRows(scen,dimsSlacks.inner.numSecondStageCons(scen)-originalDimsSlacks.inner.numSecondStageCons(scen));	
				modelChanged=true;
			}
			if (originalDimsSlacks.inner.numSecondStageVars(scen) < dimsSlacks.inner.numSecondStageVars(scen)){
				rootSolver.deleteLastSecondStageConsecutiveColumns(scen,dimsSlacks.inner.numSecondStageVars(scen)-originalDimsSlacks.inner.numSecondStageVars(scen));	
				modelChanged=true;
			}
		}

	}

	if (originalDimsSlacks.inner.numFirstStageCons() < dimsSlacks.inner.numFirstStageCons()){
		rootSolver.deleteLastFirstStageConsecutiveRows(dimsSlacks.inner.numFirstStageCons()-originalDimsSlacks.inner.numFirstStageCons());	
		modelChanged=true;
	}
	if (originalDimsSlacks.inner.numFirstStageVars() < dimsSlacks.inner.numFirstStageVars()){
		rootSolver.deleteLastFirstStageConsecutiveColumns(dimsSlacks.inner.numFirstStageVars()-originalDimsSlacks.inner.numFirstStageVars());	
		modelChanged=true;
	}
	if (modelChanged){
		BBSMPSSolver::instance()->commitNewColsAndRows();
	}
	currentlyAppliedPlanes.clear();

}


void BBSMPSTree::runParallelSBInitialization(){

	//set up the branching rule
	double beginningTime=MPI_Wtime();
	branchingRuleManager.freeResources();

    BBSMPSParallelPseudoCostBranchingRule *ppcbr= new BBSMPSParallelPseudoCostBranchingRule(100);
    branchingRuleManager.addBranchingRule(ppcbr);

	int mype=BBSMPSSolver::instance()->getMype();
	//cout<<"ABOUT TO START "<<heap.size()<<endl;
	/* While heap not empty and there are still nodes in tree */
	// TODO: Add tolerance on optimality gap, time limit option.
	while (heap.size()<1000) {
		
		int BBSMPSMyPe=BBSMPSSolver::instance()->getSBBMype();
		//cout<<" ABOUT TO ENTER "<<BBSMPSMyPe<<" "<<mype<<" "<<( (counterToLastIter%20==0 && commDone) || heap.size()==0 || (!commDone && heap.size()>16))<<endl;
		PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
		
		/* If heap is empty, update status to Stopped (if possible) and break. */
		if (inTermination) {
			//if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Heap is empty.";
			status.setStatusToStopped();

		/* If solver status is primal feasible, and the heap is empty, then
		the solution must be optimal. */
			if (status.isPrimalFeasible()) {
				status.setStatusToOptimal();
				objLB=objUB;
				if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Optimal solution found.";
			}

		/* If solver status is not primal feasible, then the MILP must be
		infeasible. */
		// TODO: Add test for unboundedness.
			if (status.isLoadedFromFile()) {
				status.setStatusToProvenInfeasible();
				if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "MILP is infeasible.";
			}

			break;
		}
		
		/* Get top-most node and pop it off of heap. */
		BBSMPSNode *currentNode_ptr=*(heap.begin());
		if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Copying node " << currentNode_ptr->getNodeNumber() << " off tree.";
	//	cout<<"COPYING NODE "<<heap.size()<<endl;
		/*int dblSerialSize=0;
		int intSerialSize=0;
		currentNode_ptr->getSerializationSize(intSerialSize,dblSerialSize);
		vector<int> serialVectorInt(intSerialSize);
		vector<double> serialVectorDbl(dblSerialSize);
		currentNode_ptr->serialize(&serialVectorInt[0],&serialVectorDbl[0]);
		cout<<"//--------------"<<endl;
		for (int i=0; i< serialVectorInt.size(); i++)cout<<serialVectorInt[i]<<" ";
			cout<<"//--------------"<<endl;
		for (int i=0; i< serialVectorDbl.size(); i++)cout<<serialVectorDbl[i]<<" ";	
			cout<<"//--------------"<<endl;*/
		heap.erase(heap.begin());

		if (nodesel == BestBound) {
			objLB=currentNode_ptr->getParentObjective();
			if ((objLB - compTol) >= objUB) {
				objLB=objUB;
				if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Can stop if best bound node selection rule";
				status.setStatusToOptimal();
				if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "All nodes can be fathomed! Terminating.";
				currentNode_ptr->eliminate();
				std::set<BBSMPSNode*>::iterator it=heap.begin();
				for (it; it!=heap.end(); it++){
					(*it)->eliminate();
					
				}
				heap.clear();
				continue; 
			}
			
		}

		vector<int> nodeCuttingPlaneUids;
		currentNode_ptr->getAllCuttingUids(nodeCuttingPlaneUids);
		if (nodeCuttingPlaneUids.size()!=currentlyAppliedPlanes.size() || !equal(nodeCuttingPlaneUids.begin(), nodeCuttingPlaneUids.begin() + nodeCuttingPlaneUids.size(), currentlyAppliedPlanes.begin())){
			removeCuts();
			std::vector<BBSMPSCuttingPlane*> cpVector;
			currentNode_ptr->getAllCuttingPlanes(cpVector);
			for (int i=0; i< cpVector.size(); i++){
				cpVector[i]->applyCuttingPlane();
			}
			BBSMPSSolver::instance()->commitNewColsAndRows();
			currentlyAppliedPlanes=nodeCuttingPlaneUids;
		}

		/* Set bounds of LP decision variables from BBSMPSNode */
		if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Setting bounds for LP subproblem.";
		//if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Parent objective of this node "<< currentNode_ptr->getParentObjective();
		//if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Parent pointer of this node "<< (currentNode_ptr->getParentPtr()!=NULL);

		denseBAVector lb(BBSMPSSolver::instance()->getOriginalLB());
		denseBAVector ub(BBSMPSSolver::instance()->getOriginalUB());
		if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "getting branching info.";
		
		currentNode_ptr->getAllBranchingInformation(lb,ub);
		if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Setting bounds.";
		
		rootSolver.setLB(lb);
		rootSolver.setUB(ub);

		/* Set information on basic/nonbasic variables for warm starting */

		BAFlagVector<variableState> ps(BBSMPSSolver::instance()->getOriginalWarmStart());

		currentNode_ptr->reconstructWarmStartState(ps);

		rootSolver.setStates(ps);

		rootSolver.commitStates();
	//	cout<<"RUN SOME PROBLEM "<<heap.size()<<endl;
		/* Solve LP defined by current node*/
		if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Solving LP subproblem.";
		rootSolver.go();
//cout<<"we got here -2 "<<heap.size()<<endl;
		/* Check solver status for infeasibility/optimality */
		solverState lpStatus = rootSolver.getStatus();

		// Only realistic solver states upon completion:
		// ProvenInfeasible, Optimal, ProvenUnbounded
		// Other solver states are intermediate states that should not
		// hold upon return from rootSolver.
		if (0 == mype) outputLPStatus(lpStatus);

		bool isLPinfeasible = (ProvenInfeasible == lpStatus); 
		//if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "isLPinfeasible = " << isLPinfeasible;
		bool isLPunbounded = (ProvenUnbounded == lpStatus);
		//if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "isLPunbounded = " << isLPunbounded ;
		bool isLPoptimal = (Optimal == lpStatus);
		//if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "isLPoptimal = " << isLPoptimal;
		bool isLPother = (!isLPinfeasible && !isLPunbounded && !isLPoptimal);
		//if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "isLPother = " << isLPother;
		assert (!isLPother); // Error if not infeasible/unbounded/optimal

		/* Fathom by infeasibility */
		// If LP solver returns infeasibility, fathom node, go to start of loop
		if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Checking for infeasibility...";
		if (isLPinfeasible) {
			if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Fathoming node " << currentNode_ptr->getNodeNumber() << " by infeasibility.";
			currentNode_ptr->eliminate();
			nodesFathomed++;
		//cout<<"isLPinfeasible "<<heap.size()<<endl;
			
			continue;
		}

		// Otherwise, LP is feasible. LP may be optimal or unbounded.
		// If LP is unbounded, so is the MILP.
		if (isLPunbounded) {
			if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "LP relaxation of node " << currentNode_ptr->getNodeNumber()
				<< " is unbounded.\n"
			<< "Please add additional constraints to "
			<< "bound the MILP.";
		//	cout<<"isLPunbounded "<<heap.size()<<endl;
			return;
		}

		// At this point, LP must be optimal.
		assert (isLPoptimal); // Error if not optimal.
//cout<<"we got here -1 "<<heap.size()<<endl;
		denseBAVector primalSoln(rootSolver.getPrimalSolution());
		
		bool newCuttingPlanes=cuttingPlanesManager.generateCuttingPlanes(currentNode_ptr,primalSoln);

		if(newCuttingPlanes){
			currentNode_ptr->getAllCuttingUids(currentlyAppliedPlanes);
		}
		lpStatus = rootSolver.getStatus();
		primalSoln=denseBAVector(rootSolver.getPrimalSolution());
		isLPoptimal = (Optimal == lpStatus);
		ps=BAFlagVector<variableState>(BBSMPSSolver::instance()->getOriginalWarmStart());
		currentNode_ptr->reconstructWarmStartState(ps);	
		// TODO: Combine the integrality and branching steps later

		/* If LP solution is optimal, can fathom by value dominance. */
		assert (isLPoptimal); // Error if not optimal.
		if (isLPoptimal) {
			// If LP solver returns optimal, then the objective is bounded below.
			// TODO: Change solver status to "Bounded".
			//	setStatusToBounded();

			// Since the branch-and-bound tree is stored as a min-heap, the current
			// node being explored always has the minimal objective function value.
			// If the branch-and-bound tree is ever re-heapified so that it is NOT
			// a min-heap, but has the heap property for some other ordering, then
			// this update cannot occur without taking the min objective function
			// value over all values of the parent objective function for nodes
			// still in the B&B tree (which is expensive if it is not the key used
			// to heapify the min-heap).
			//objLB = currentNode.parentObj;
			//if ((lpObj - compTol) >= objLB) {
			//  if (0 == mype) //cout << "Current best lower bound is " << objLB << endl;
			//  if (0 == mype) //cout << "Updating best lower bound to " << lpObj << endl;
			//  objLB = lpObj;
			//}

			/* Fathom by value dominance */
			// Optimal LP objective function value is lower bound on the objective
			// function value of the LP derived from any node in the subtree
			// of the B&B tree rooted at the current node, so no feasible
			// solution in that subtree can have a lesser objective function
			// value than the current upper bound on the optimal value of
			// the MILP objective function.

			// Get LP objective function value.
			if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Getting LP objective...";
			double lpObj = rootSolver.getObjective();

			currentNode_ptr->setObjective(lpObj);

			if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Checking for value dominance...";

			if ((lpObj - compTol) >= objUB) {
				if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Fathoming node " << currentNode_ptr->getNodeNumber() << " by value dominance.";

				currentNode_ptr->eliminate();
				nodesFathomed++;
				
				continue;
			}
				
			
			
		}
//cout<<"we got here "<<heap.size()<<endl;
		
		/* Get primal solution */
		//if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Getting primal solution...";
		
		denseBAVector originalSpaceSolution;
		originalSpaceSolution.allocate(BBSMPSSolver::instance()->getOriginalBADimensionsSlacks(), BBSMPSSolver::instance()->getBAContext(), PrimalVector); 
		originalSpaceSolution.copyAndShrinkToDims(primalSoln);
		
		/* If primal solution is integral: */
		//  - Update solver status to PrimalFeasible
		//  - Check if upper bound improved
		//  - If so, update current primal solution for that upper bound
		//  - Fathom node, go to start of loop

		if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Checking for integrality of primal solution...";
		if(isLPIntFeas(primalSoln)) {
			if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Node " << currentNode_ptr->getNodeNumber() << " is integer feasible.";
			status.setStatusToPrimalFeasible();

			/* Update upper bound if it's less than current best upper bound, and
			the LP solution is optimal (not unbounded). */
			double newUB = rootSolver.getObjective();
			bool isNewUBbetter = (newUB < (objUB - compTol));
			if (isLPoptimal && isNewUBbetter) {
				if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Updating best upper bound to " << newUB ;
				objUB = rootSolver.getObjective();
				ubPrimalSolution.copyFrom(originalSpaceSolution);
				BBSMPSSolution aux(originalSpaceSolution,newUB, BBSMPSSolver::instance()->getWallTime());
				BBSMPSSolver::instance()->addSolutionToPool(aux);
			}
			currentNode_ptr->eliminate();
			

			nodesBecameInteger++;
			
		}
		else{
		//	cout<<"we got here 2 "<<heap.size()<<endl;
			double lpObj = rootSolver.getObjective();
			const BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
			BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
			BAFlagVector<variableState> states(dimsSlacks, ctx, PrimalVector);
			rootSolver.getStates(states);
	       
			vector<BBSMPSNode*> children;
		//	cout<<"calling branch!!! "<<heap.size()<<endl;
			branchingRuleManager.branch(currentNode_ptr,children,originalSpaceSolution);

			if (children.size()>0){
				for (int i=0; i<children.size();i++){
					generateIncrementalWarmState(children[i], ps, states);
					children[i]->setObjective(lpObj);

					heap.insert(children[i]);
				}

			}
			//removeCuts();
			heuristicsManager.runLPHeuristics(currentNode_ptr,primalSoln);
		

			if(heuristicsManager.willMIPHeuristicsRun(currentNode_ptr,originalSpaceSolution)){
				removeCuts();
				heuristicsManager.runMIPHeuristics(currentNode_ptr,originalSpaceSolution);
				//removeCuts();
			}

		}

		/* Optimality gap termination criteria: if gap between best
		upper bound and best lower bound on objective function is
		less than the optimality gap, stop the solver and return the
		feasible solution corresponding to the best upper bound. */
		
		double gap = fabs(objUB-objLB)*100/(fabs(objUB)+10e-10);
		/*if (gap <= optGapTol) {
			assert(objLB <= objUB); // this statement could be tripped by numerical error
			status.setStatusToOptimal();
			if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Optimality gap reached! Terminating.";
			break;
		}*/
		

		if (BBSMPSSolver::instance()->getSolPoolSize()>0) objUB=min(objUB,BBSMPSSolver::instance()->getSoln(0).getObjValue());

		bbIterationCounter++;
		if (0 == mype && verbosityActivated && bbIterationCounter%500==0) {
			double gap = fabs(objUB-objLB)*100/(fabs(objUB)+10e-10);
			BBSMPS_ALG_LOG_SEV(warning)<<"\n----------------------------------------------------\n"<<
			"Iteration "<<bbIterationCounter<<":LB:"<<objLB<<":UB:"<<objUB<<":GAP:"<<gap<<":Tree Size:"<<heap.size()<<":Time:"<<BBSMPSSolver::instance()->getWallTime()<<"\n"<<
			"----------------------------------------------------";
		}
	}

	
	denseBAVector v1,v2,v3,v4;
	ppcbr->getCostHistory(v1,v2,v3,v4);
	initializationParallelBranchingCommTime=ppcbr->getCommunicationTime();
	branchingRuleManager.freeResources();

     BBSMPSMaxFracBranchingRule *mfbr= new BBSMPSMaxFracBranchingRule(10);
    branchingRuleManager.addBranchingRule(mfbr);

    //BBSMPSPseudoCostBranchingRule *mfbr2= new BBSMPSPseudoCostBranchingRule(100);
   //mfbr2->setCostHistory(v1,v2,v3,v4);
  	 // branchingRuleManager.addBranchingRule(mfbr2);
  	  initializationTime=MPI_Wtime()-beginningTime;

}



		  // TODO: Add way to access single scenario.

		  /* TODO: Refactor MIP solver status into its own class to handle
		     the transition logic. */
		  /* TODO: Add more nuanced version of solver statuses. */
		  /* Methods defining allowable transitions for status */
		  /* status transition state diagram */
		  // NOTE: "Bounded" is not currently a solver state;
		  // TODO: Add Bounded as a solver state.
		  // LoadedFromFile -> {Bounded, PrimalFeasible, Optimal,
		  //                    ProvenInfeasible, Unbounded, Stopped}
		  // Bounded -> {PrimalFeasible, Optimal, ProvenInfeasible, Stopped}
		  // PrimalFeasible -> {Optimal, Stopped}

		  // void setStatusToBounded() {
		  //   bool isInReachableState = (LoadedFromFile == status);
		  //   if (isInReachableState) {
		  //     status = Bounded;
		  //   }
		  // }


void BBSMPSTree::branchAndBound() {

	int BBSMPSMyPe=BBSMPSSolver::instance()->getSBBMype();
   if (BBSMPSMyPe!=0){
   	objLB=COIN_DBL_MAX;
   		std::set<BBSMPSNode*>::iterator it=heap.begin();
		for (it; it!=heap.end(); it++){
			(*it)->eliminate();
			
		}
		heap.clear();
   }
	

	int mype=BBSMPSSolver::instance()->getMype();
	
	bool commDone=false;
	int counterToLastIter=1;
	/* While heap not empty and there are still nodes in tree */
	// TODO: Add tolerance on optimality gap, time limit option.
	while (true) {
		
		
		//cout<<" ABOUT TO ENTER "<<BBSMPSMyPe<<" "<<mype<<" "<<( (counterToLastIter%20==0 && commDone) || heap.size()==0 || (!commDone && heap.size()>16))<<endl;
		if(!communicationActivated){
			checkSequentialTerminationConditions();
		}
		else {
			bool launch=shouldWePerformCommunication( counterToLastIter,  commDone);

			

			if (launch ){
				commDone=true;
				//BBSMPS_ALG_LOG_SEV(warning)<<bbIterationCounter<<" THREAD "<<mype<<" "<<BBSMPSMyPe<<" got into communication ";
   
				communicate(); 
			}
		//	else{
		//		BBSMPS_ALG_LOG_SEV(warning)<<bbIterationCounter<<" THREAD "<<mype<<" "<<BBSMPSMyPe<<" DID NOT!!!!! got into communication ";
			//} 
		}
		counterToLastIter++;
		PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
		
		/* If heap is empty, update status to Stopped (if possible) and break. */
		if (inTermination) {
			//if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Heap is empty.";
			status.setStatusToStopped();

		/* If solver status is primal feasible, and the heap is empty, then
		the solution must be optimal. */
			if (status.isPrimalFeasible()) {
				status.setStatusToOptimal();
				objLB=objUB;
				if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Optimal solution found.";
			}

		/* If solver status is not primal feasible, then the MILP must be
		infeasible. */
		// TODO: Add test for unboundedness.
			if (status.isLoadedFromFile()) {
				status.setStatusToProvenInfeasible();
				if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "MILP is infeasible.";
			}

			break;
		}
		
		if (heap.size()==0) continue;

		/* Get top-most node and pop it off of heap. */
		BBSMPSNode *currentNode_ptr=*(heap.begin());
		if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Copying node " << currentNode_ptr->getNodeNumber() << " off tree.";
//BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
		//MPI_Barrier(ctx.comm());
		//BBSMPS_ALG_LOG_SEV(warning)<<"WORKING ON NODE "<<currentNode_ptr->getNodeNumber();
		/*int dblSerialSize=0;
		int intSerialSize=0;
		currentNode_ptr->getSerializationSize(intSerialSize,dblSerialSize);
		vector<int> serialVectorInt(intSerialSize);
		vector<double> serialVectorDbl(dblSerialSize);
		currentNode_ptr->serialize(&serialVectorInt[0],&serialVectorDbl[0]);
		cout<<"//--------------"<<endl;
		for (int i=0; i< serialVectorInt.size(); i++)cout<<serialVectorInt[i]<<" ";
			cout<<"//--------------"<<endl;
		for (int i=0; i< serialVectorDbl.size(); i++)cout<<serialVectorDbl[i]<<" ";	
			cout<<"//--------------"<<endl;*/
		heap.erase(heap.begin());

		if (nodesel == BestBound) {
			objLB=currentNode_ptr->getParentObjective();
			if ((objLB - compTol) >= objUB) {
				objLB=objUB;
				if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Can stop if best bound node selection rule";
				status.setStatusToOptimal();
				if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "All nodes can be fathomed! Terminating.";
				currentNode_ptr->eliminate();
				std::set<BBSMPSNode*>::iterator it=heap.begin();
				for (it; it!=heap.end(); it++){
					(*it)->eliminate();
					
				}
				heap.clear();
				continue; 
			}
			
		}

		vector<int> nodeCuttingPlaneUids;
		currentNode_ptr->getAllCuttingUids(nodeCuttingPlaneUids);
		if (nodeCuttingPlaneUids.size()!=currentlyAppliedPlanes.size() || !equal(nodeCuttingPlaneUids.begin(), nodeCuttingPlaneUids.begin() + nodeCuttingPlaneUids.size(), currentlyAppliedPlanes.begin())){
			removeCuts();
			std::vector<BBSMPSCuttingPlane*> cpVector;
			currentNode_ptr->getAllCuttingPlanes(cpVector);
			for (int i=0; i< cpVector.size(); i++){
				cpVector[i]->applyCuttingPlane();
			}
			BBSMPSSolver::instance()->commitNewColsAndRows();
			currentlyAppliedPlanes=nodeCuttingPlaneUids;
		}

		/* Set bounds of LP decision variables from BBSMPSNode */
		if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Setting bounds for LP subproblem.";
		//if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Parent objective of this node "<< currentNode_ptr->getParentObjective();
		//if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Parent pointer of this node "<< (currentNode_ptr->getParentPtr()!=NULL);

		denseBAVector lb(BBSMPSSolver::instance()->getOriginalLB());
		denseBAVector ub(BBSMPSSolver::instance()->getOriginalUB());
		if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "getting branching info.";
		
		currentNode_ptr->getAllBranchingInformation(lb,ub);
		if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Setting bounds.";
		
		rootSolver.setLB(lb);
		rootSolver.setUB(ub);

		/* Set information on basic/nonbasic variables for warm starting */

		BAFlagVector<variableState> ps(BBSMPSSolver::instance()->getOriginalWarmStart());

		currentNode_ptr->reconstructWarmStartState(ps);

		rootSolver.setStates(ps);

		rootSolver.commitStates();
		
		/* Solve LP defined by current node*/
		if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Solving LP subproblem.";
		rootSolver.go();

		/* Check solver status for infeasibility/optimality */
		solverState lpStatus = rootSolver.getStatus();

		// Only realistic solver states upon completion:
		// ProvenInfeasible, Optimal, ProvenUnbounded
		// Other solver states are intermediate states that should not
		// hold upon return from rootSolver.
		if (0 == mype) outputLPStatus(lpStatus);

		bool isLPinfeasible = (ProvenInfeasible == lpStatus); 
		//if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "isLPinfeasible = " << isLPinfeasible;
		bool isLPunbounded = (ProvenUnbounded == lpStatus);
		//if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "isLPunbounded = " << isLPunbounded ;
		bool isLPoptimal = (Optimal == lpStatus);
		//if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "isLPoptimal = " << isLPoptimal;
		bool isLPother = (!isLPinfeasible && !isLPunbounded && !isLPoptimal);
		//if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "isLPother = " << isLPother;
		assert (!isLPother); // Error if not infeasible/unbounded/optimal

		/* Fathom by infeasibility */
		// If LP solver returns infeasibility, fathom node, go to start of loop
		if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Checking for infeasibility...";
		if (isLPinfeasible) {
			if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Fathoming node " << currentNode_ptr->getNodeNumber() << " by infeasibility.";
			currentNode_ptr->eliminate();
			nodesFathomed++;
	
			
			continue;
		}

		// Otherwise, LP is feasible. LP may be optimal or unbounded.
		// If LP is unbounded, so is the MILP.
		if (isLPunbounded) {
			if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "LP relaxation of node " << currentNode_ptr->getNodeNumber()
				<< " is unbounded.\n"
			<< "Please add additional constraints to "
			<< "bound the MILP.";
			
			return;
		}

		// At this point, LP must be optimal.
		assert (isLPoptimal); // Error if not optimal.

		denseBAVector primalSoln(rootSolver.getPrimalSolution());
		
		bool newCuttingPlanes=cuttingPlanesManager.generateCuttingPlanes(currentNode_ptr,primalSoln);

		if(newCuttingPlanes){
			currentNode_ptr->getAllCuttingUids(currentlyAppliedPlanes);
		}
		lpStatus = rootSolver.getStatus();
		primalSoln=denseBAVector(rootSolver.getPrimalSolution());
		isLPoptimal = (Optimal == lpStatus);
		ps=BAFlagVector<variableState>(BBSMPSSolver::instance()->getOriginalWarmStart());
		currentNode_ptr->reconstructWarmStartState(ps);	
		// TODO: Combine the integrality and branching steps later

		/* If LP solution is optimal, can fathom by value dominance. */
		assert (isLPoptimal); // Error if not optimal.
		if (isLPoptimal) {
			// If LP solver returns optimal, then the objective is bounded below.
			// TODO: Change solver status to "Bounded".
			//	setStatusToBounded();

			// Since the branch-and-bound tree is stored as a min-heap, the current
			// node being explored always has the minimal objective function value.
			// If the branch-and-bound tree is ever re-heapified so that it is NOT
			// a min-heap, but has the heap property for some other ordering, then
			// this update cannot occur without taking the min objective function
			// value over all values of the parent objective function for nodes
			// still in the B&B tree (which is expensive if it is not the key used
			// to heapify the min-heap).
			//objLB = currentNode.parentObj;
			//if ((lpObj - compTol) >= objLB) {
			//  if (0 == mype) //cout << "Current best lower bound is " << objLB << endl;
			//  if (0 == mype) //cout << "Updating best lower bound to " << lpObj << endl;
			//  objLB = lpObj;
			//}

			/* Fathom by value dominance */
			// Optimal LP objective function value is lower bound on the objective
			// function value of the LP derived from any node in the subtree
			// of the B&B tree rooted at the current node, so no feasible
			// solution in that subtree can have a lesser objective function
			// value than the current upper bound on the optimal value of
			// the MILP objective function.

			// Get LP objective function value.
			if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Getting LP objective...";
			double lpObj = rootSolver.getObjective();

			currentNode_ptr->setObjective(lpObj);

			if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Checking for value dominance...";

			if ((lpObj - compTol) >= objUB) {
				if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Fathoming node " << currentNode_ptr->getNodeNumber() << " by value dominance.";

				currentNode_ptr->eliminate();
				nodesFathomed++;
				
				continue;
			}
				
			
			
		}

	
		/* Get primal solution */
		//if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Getting primal solution...";
		
		denseBAVector originalSpaceSolution;
		originalSpaceSolution.allocate(BBSMPSSolver::instance()->getOriginalBADimensionsSlacks(), BBSMPSSolver::instance()->getBAContext(), PrimalVector); 
		originalSpaceSolution.copyAndShrinkToDims(primalSoln);
		
		/* If primal solution is integral: */
		//  - Update solver status to PrimalFeasible
		//  - Check if upper bound improved
		//  - If so, update current primal solution for that upper bound
		//  - Fathom node, go to start of loop

		if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Checking for integrality of primal solution...";
		if(isLPIntFeas(primalSoln)) {
			if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Node " << currentNode_ptr->getNodeNumber() << " is integer feasible.";
			status.setStatusToPrimalFeasible();

			/* Update upper bound if it's less than current best upper bound, and
			the LP solution is optimal (not unbounded). */
			double newUB = rootSolver.getObjective();
			bool isNewUBbetter = (newUB < (objUB - compTol));
			if (isLPoptimal && isNewUBbetter) {
				if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Updating best upper bound to " << newUB ;
				objUB = rootSolver.getObjective();
				ubPrimalSolution.copyFrom(originalSpaceSolution);
				BBSMPSSolution aux(originalSpaceSolution,newUB, BBSMPSSolver::instance()->getWallTime());
				BBSMPSSolver::instance()->addSolutionToPool(aux);
			}
			currentNode_ptr->eliminate();
			

			nodesBecameInteger++;
			
		}
		else{
			double lpObj = rootSolver.getObjective();
			const BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
			BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
			BAFlagVector<variableState> states(dimsSlacks, ctx, PrimalVector);
			rootSolver.getStates(states);
	       
			vector<BBSMPSNode*> children;

			branchingRuleManager.branch(currentNode_ptr,children,originalSpaceSolution);

			if (children.size()>0){
				for (int i=0; i<children.size();i++){
					generateIncrementalWarmState(children[i], ps, states);
					children[i]->setObjective(lpObj);

					heap.insert(children[i]);
				}

			}
			//removeCuts();
			heuristicsManager.runLPHeuristics(currentNode_ptr,primalSoln);
		

			if(heuristicsManager.willMIPHeuristicsRun(currentNode_ptr,originalSpaceSolution)){
				removeCuts();
				heuristicsManager.runMIPHeuristics(currentNode_ptr,originalSpaceSolution);
				//removeCuts();
			}

		}

		/* Optimality gap termination criteria: if gap between best
		upper bound and best lower bound on objective function is
		less than the optimality gap, stop the solver and return the
		feasible solution corresponding to the best upper bound. */
		
		double gap = fabs(objUB-objLB)*100/(fabs(objUB)+10e-10);
		/*if (gap <= optGapTol) {
			assert(objLB <= objUB); // this statement could be tripped by numerical error
			status.setStatusToOptimal();
			if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Optimality gap reached! Terminating.";
			break;
		}*/
		

		if (BBSMPSSolver::instance()->getSolPoolSize()>0) objUB=min(objUB,BBSMPSSolver::instance()->getSoln(0).getObjValue());

		bbIterationCounter++;
		if (0 == mype && verbosityActivated && bbIterationCounter%200==0) {
			double gap = fabs(objUB-objLB)*100/(fabs(objUB)+10e-10);
			BBSMPS_ALG_LOG_SEV(warning)<<"\n----------------------------------------------------\n"<<
			"Iteration "<<bbIterationCounter<<":LB:"<<objLB<<":UB:"<<objUB<<":GAP:"<<gap<<":Tree Size:"<<heap.size()<<":Time:"<<BBSMPSSolver::instance()->getWallTime()<<"\n"<<
			"----------------------------------------------------";
		}
	}
	//int BBSMPSMyPe=BBSMPSSolver::instance()->getSBBMype();
  		
//cout<<"processor "<<BBSMPSMyPe<<" got out!!!!!!!"<<endl;
//if (0 == mype) BBSMPS_ALG_LOG_SEV(summary) << "Objective function value = " << objUB ;
//if (0 == mype) BBSMPS_ALG_LOG_SEV(summary) << "Objective function LB = " << objLB ;
	if (0 == mype && verbosityActivated) {

		
		BAContext BBSMPSContext=BBSMPSSolver::instance()->getSBBContext();
		int BBSMPSProcs=BBSMPSContext.nprocs();
		double generalBestUB;
		double generalBestLB;
		double idleCommTimeSum;
		double CommTimeSum;
		int totalNodes;
		double totalInitTime;
		double totalInitCommTime;


		MPI_Allreduce(&objUB,&generalBestUB, 1, MPI_DOUBLE,MPI_MIN,BBSMPSContext.comm());
		MPI_Allreduce(&objLB,&generalBestLB, 1, MPI_DOUBLE,MPI_MIN,BBSMPSContext.comm());
		MPI_Allreduce(&totalIdleTime,&idleCommTimeSum, 1, MPI_DOUBLE,MPI_SUM,BBSMPSContext.comm());
		MPI_Allreduce(&totalCommTime,&CommTimeSum, 1, MPI_DOUBLE,MPI_SUM,BBSMPSContext.comm());
		MPI_Allreduce(&bbIterationCounter,&totalNodes, 1, MPI_INT,MPI_SUM,BBSMPSContext.comm());
		MPI_Allreduce(&initializationTime,&totalInitTime, 1, MPI_DOUBLE,MPI_SUM,BBSMPSContext.comm());
		MPI_Allreduce(&initializationParallelBranchingCommTime,&totalInitCommTime, 1, MPI_DOUBLE,MPI_SUM,BBSMPSContext.comm());

		if (BBSMPSMyPe==0){
			double generalGap = fabs(generalBestUB-generalBestLB)*100/(fabs(generalBestUB)+10e-10);
			BBSMPS_ALG_LOG_SEV(warning)<<"-------------EXPLORATION TERMINATED---------------";
			BBSMPS_ALG_LOG_SEV(warning)<<"General GAP:"<<generalGap<<":LB:"<<generalBestLB<<":UB:"<<generalBestUB<<":Nodes Solved:"<<totalNodes;
			BBSMPS_ALG_LOG_SEV(warning)<<"Number of Communication calls:"<<totalTimesCommCalled<<":Avg. Communication Time:"<<CommTimeSum/BBSMPSProcs<<":Avg. Idle Communication Time:"<<idleCommTimeSum/BBSMPSProcs;
			BBSMPS_ALG_LOG_SEV(warning)<<"Avg. Initialization Time:"<<totalInitTime/BBSMPSProcs<<":Avg. Initialization Communication Time:"<<totalInitCommTime/BBSMPSProcs;
		}
		double gap = fabs(objUB-objLB)*100/(fabs(objUB)+10e-10);
  		
  		for (int i=0;i< BBSMPSProcs; i++){
  			MPI_Barrier(BBSMPSContext.comm());

  			if (BBSMPSMyPe==i){
  				
  				BBSMPS_ALG_LOG_SEV(warning)<<"---------------------SOLVER "<<BBSMPSMyPe<<"---------------------";
				BBSMPS_ALG_LOG_SEV(warning)<<"Iteration "<<bbIterationCounter<<":LB:"<<objLB<<":UB:"<<objUB<<":GAP:"<<gap<<":Tree Size:"<<heap.size();
				BBSMPS_ALG_LOG_SEV(warning)<<"Nodes Fathomed:"<<nodesFathomed<<":Nodes with integer Solution:"<<nodesBecameInteger;
				BBSMPS_ALG_LOG_SEV(warning)<<"LP Relaxation Value:"<<LPRelaxationValue<<":LP Relaxation Time:"<<LPRelaxationTime<<":Preprocessing Time:"<<PreProcessingTime<<":Total Time:"<<BBSMPSSolver::instance()->getWallTime();
				BBSMPS_ALG_LOG_SEV(warning)<<"Number of Communication calls:"<<totalTimesCommCalled<<":Communication Time:"<<totalCommTime<<":Communication Idle Time"<<totalIdleTime;
				BBSMPS_ALG_LOG_SEV(warning)<<"Parallel Strong Branch Init Time:"<<initializationTime<<":Parallel Strong Branching Comm. Time:"<<initializationParallelBranchingCommTime;
				

				BBSMPS_ALG_LOG_SEV(warning)<<"----------------------------------------------------";
				heuristicsManager.printStatistics();
				branchingRuleManager.printStatistics();
				cuttingPlanesManager.printStatistics();
				BBSMPSSolver::instance()->printSolutionStatistics(objLB);
				BBSMPSSolver::instance()->printPresolveStatistics( );
  			}
		}
		MPI_Barrier(BBSMPSContext.comm());
	}

	double t = BBSMPSSolver::instance()->getWallTime();
	if (0 == mype && verbosityActivated) {
		BBSMPS_APP_LOG_SEV(warning)<<boost::format("Branch and Bound took %f seconds") % t;
	}
}


void BBSMPSTree::setTimeLimit(int _tiLim){
	tiLim=_tiLim;
}

void BBSMPSTree::setNodeLimit(int _nodeLim){
	nodeLim=_nodeLim;
}



  void BBSMPSTree::loadSimpleHeuristics(){
  	
 	
	/*BBSMPSHeuristicRounding *hr= new BBSMPSHeuristicRounding(0,1,"SimpleRounding");
   heuristicsManager.addLPHeuristic(hr);

    BBSMPSHeuristicFixAndDive *hr2= new BBSMPSHeuristicFixAndDive(0,1,"FixAndDive");
	heuristicsManager.addLPHeuristic(hr2);

  	BBSMPSHeuristicLockRounding *hr3= new BBSMPSHeuristicLockRounding(0,1,"LockRounding");
  	heuristicsManager.addLPHeuristic(hr3);

  	
*/

  	BBSMPSHeuristicLockRounding *hr6= new BBSMPSHeuristicLockRounding(0,3,"HeuristicLockRounding");

  heuristicsManager.addLPHeuristic(hr6);

    BBSMPSHeuristicFixAndDiveLocks *hr5= new BBSMPSHeuristicFixAndDiveLocks(0,3,"FixAndDiveLocks");

   heuristicsManager.addLPHeuristic(hr5);  
  }

  void BBSMPSTree::loadCuttingPlanes(){
	BBSMPSCuttingPlaneGenerator01KP *plane=new BBSMPSCuttingPlaneGenerator01KP("01kp");
    cuttingPlanesManager.addCuttingPlaneGenerator(plane);
  }
  void BBSMPSTree::loadMIPHeuristics(){
 //  	BBSMPSHeuristicRINS *hr= new BBSMPSHeuristicRINS(30,50,"RINS",200);
//	BBSMPSHeuristicRENS *hr2= new BBSMPSHeuristicRENS(30,50,"RENS",200);
  /*	BBSMPSHeuristicCrossover *hr3= new BBSMPSHeuristicCrossover(0,1,"Crossover",1);
  	heuristicsManager.addMIPHeuristic(hr3);*/
 // heuristicsManager.addMIPHeuristic(hr);
  // heuristicsManager.addMIPHeuristic(hr2);
  /*
     BBSMPSHeuristicSolutionRINS *hr4= new BBSMPSHeuristicSolutionRINS(0,1,"SolRINS",1);
   	heuristicsManager.addMIPHeuristic(hr4);
 
  
	 BBSMPSHeuristicBestRINSJump *hr5= new BBSMPSHeuristicBestRINSJump(0,1,"BestRinsJump",1);
  	heuristicsManager.addMIPHeuristic(hr5);
    BBSMPSHeuristicSolutionPolishing *hr6= new BBSMPSHeuristicSolutionPolishing(0,1,"BBSMPSHeuristicSolutionPolishing",1);
  	heuristicsManager.addMIPHeuristic(hr6);*/

  }


  void BBSMPSTree::loadLPHeuristic(BBSMPSHeuristic *heur){
  	heuristicsManager.addLPHeuristic(heur);
  }
   void BBSMPSTree::loadMIPHeuristic(BBSMPSHeuristic *heur){
  	heuristicsManager.addMIPHeuristic(heur);
  }


  BBSMPSNode* BBSMPSTree::topOfHeap(){
  	return *(heap.begin());
  }


  void BBSMPSTree::setVerbosity(bool verbose){
  	verbosityActivated=verbose;
  }


void BBSMPSTree::setSolLimit(int _solLim){
	solsDiscoveredInit=BBSMPSSolver::instance()->getSolPoolSize();

   solsDiscoveredLimit=_solLim;

}

void BBSMPSTree::setGAPTolLimit( double _GAPTolLim){
	optGapTol=_GAPTolLim;
}

int BBSMPSTree::communicate(){

	double startCheckpoint=MPI_Wtime();
	 BAContext &BBSMPSContext=BBSMPSSolver::instance()->getSBBContext();
	int BBSMPSProcs=BBSMPSContext.nprocs();
	int BBSMPSMyPe=BBSMPSSolver::instance()->getSBBMype();

	int mype=BBSMPSSolver::instance()->getMype();


int nSolsExchanged=BBSMPSProcs;
	BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
	vector<double> allLBsGathered(BBSMPSProcs);
	double LBToSend=objLB;

	int totalSumOfNodes = heap.size();
	MPI_Allreduce(MPI_IN_PLACE,&totalSumOfNodes, 1, MPI_INT,MPI_SUM,BBSMPSContext.comm());
	if (totalSumOfNodes<BBSMPSProcs)return 0;
	if (heap.size()==0)LBToSend=COIN_DBL_MAX;

	//BBSMPS_ALG_LOG_SEV(warning)<<"got out of barrier ";
   
	MPI_Allgather(&LBToSend,1,MPI_DOUBLE, &allLBsGathered[0],1,MPI_DOUBLE,BBSMPSContext.comm());
	
	/*if(mype==0){
			for (int i=0; i<allLBsGathered.size();i++){
				cout<<"["<<i<<","<<allLBsGathered[i]<<"]";
			}
			cout<<endl;
	}*/


	double bestLB=COIN_DBL_MAX;
	for (int i=0; i<BBSMPSProcs; i++) if (bestLB>allLBsGathered[i]) bestLB=allLBsGathered[i];
	int doINeedComm=0;
	if (abs(allLBsGathered[BBSMPSMyPe]-bestLB)>=commTol)doINeedComm=1;
	int nProcsThatNeedComm=0;
	for (int i=0; i<BBSMPSProcs; i++)if (abs(allLBsGathered[i]-bestLB)>=commTol)nProcsThatNeedComm++;


	int totalNodesToExchange=nProcsThatNeedComm*nSolsExchanged;

//    cout<<"total nodes to exchange are "<<totalNodesToExchange<<" "<<bestLB<<endl;
	//Gather the top n x n node bounds
	vector<double> localBoundBuff(totalNodesToExchange,COIN_DBL_MAX);
	vector<double> globalBoundBuff(totalNodesToExchange*BBSMPSProcs,COIN_DBL_MAX);

	std::set<BBSMPSNode*>::iterator it=heap.begin();
	for (int i=0; i<totalNodesToExchange && it!=heap.end(); i++){
		localBoundBuff[i]=(*it)->getObjective();
		it++;
		if (it==heap.end())localBoundBuff[i]=COIN_DBL_MAX;
	}
/*	if(mype==0){
		for (int i=0;i< BBSMPSProcs; i++){
  			MPI_Barrier(BBSMPSContext.comm());

  			if (BBSMPSMyPe==i){
  				
				cout<<"before allgather "<<BBSMPSMyPe<<" ";
				for (int i=0; i<localBoundBuff.size();i++){
					cout<<"["<<i<<","<<localBoundBuff[i]<<"]";
				}
				cout<<endl;
			}
		}
	}*/
	//All to all exchange
	MPI_Allgather(&localBoundBuff[0],totalNodesToExchange,MPI_DOUBLE, &globalBoundBuff[0],totalNodesToExchange,MPI_DOUBLE,BBSMPSContext.comm());
	totalIdleTime+=(MPI_Wtime()-startCheckpoint);
	double check1=(MPI_Wtime()-startCheckpoint);
	double check2=MPI_Wtime();
	/*if(mype==0){
			for (int i=0; i<totalNodesToExchange*nProcsThatNeedComm;i++){
				cout<<"["<<i<<","<<globalBoundBuff[i]<<"]";
			}
			cout<<endl;
	}*/
	//Create vector of pairs with bound, proc number and rank within owner
	vector< pair <double, pair<int, int> > > globalBoundVector(totalNodesToExchange*BBSMPSProcs);
	for (int p=0; p<BBSMPSProcs; p++){
		for(int bds=0; bds< totalNodesToExchange; bds ++){
			globalBoundVector[p*totalNodesToExchange+bds].first=globalBoundBuff[p*totalNodesToExchange+bds];
			globalBoundVector[p*totalNodesToExchange+bds].second.first=p;
			globalBoundVector[p*totalNodesToExchange+bds].second.second=bds;
		}
	}
	//sort vector of pairs
	sort(globalBoundVector.begin(), globalBoundVector.end(),boundVectorComparator);
	

/*	if(mype==0 && BBSMPSMyPe==0){
			cout<<"the sort looks like this "<<BBSMPSMyPe<<" ";
			for (int i=0; i<globalBoundVector.size();i++){
				cout<<"["<<globalBoundVector[i].first<<","<<globalBoundVector[i].second.first<<","<<globalBoundVector[i].second.second<<"]";
			}
			cout<<endl;
	}*/
	vector<int> nodesOwned(BBSMPSProcs,0);
	for (int nod=0;nod<totalNodesToExchange;nod++){
		nodesOwned[globalBoundVector[nod].second.first]++;
	}
	vector<int> nodesOwnedDisplacement(BBSMPSProcs,0);

	for(int i=1; i< nodesOwnedDisplacement.size(); i++){
		nodesOwnedDisplacement[i]=nodesOwnedDisplacement[i-1]+nodesOwned[i-1];
	}
	//cout<<" total of nodes owned!!! "<<nodesOwned[BBSMPSMyPe]<<" "<<BBSMPSMyPe<<endl;
	//Make a vector to exchange the sizes of the nodes to exchange
	vector<int> nodeIntSizes (nodesOwned[BBSMPSMyPe],0);
	vector<int> nodeDoubleSizes (nodesOwned[BBSMPSMyPe],0);
	int ptr=0;
	it=heap.begin();
	for (int nod=0;nod<totalNodesToExchange;nod++){
		if(globalBoundVector[nod].second.first==BBSMPSMyPe){
			if (globalBoundVector[nod].first<COIN_DBL_MAX){
				int intVectorSize;
				int dblVectorSize;
				(*it)->getSerializationSize(intVectorSize, dblVectorSize);
				nodeIntSizes[ptr]=intVectorSize;
				nodeDoubleSizes[ptr]=dblVectorSize;
				ptr++;
				it++;
			}
			else{
				nodeIntSizes[ptr]=0;
				nodeDoubleSizes[ptr]=0;
			}
		}
	}
	check2=(MPI_Wtime()-check2);
	double check3=(MPI_Wtime());
	vector<int> globalNodeIntSizes (totalNodesToExchange,0);
	vector<int> globalNodeDoubleSizes (totalNodesToExchange,0);
/*	if(mype==0){
			cout<<"before allgather sizes "<<BBSMPSMyPe<<" ";
			for (int i=0; i<nodeIntSizes.size();i++){
				cout<<"["<<i<<","<<nodeIntSizes[i]<<","<<nodeDoubleSizes[i]<<"]";
			}
			cout<<endl;
	}*/
	//Exchange vector of sizes
	double comm=MPI_Wtime();
	MPI_Allgatherv(&nodeIntSizes[0],nodesOwned[BBSMPSMyPe],MPI_INT,&globalNodeIntSizes[0],&nodesOwned[0],&nodesOwnedDisplacement[0],MPI_INT,BBSMPSContext.comm());
	MPI_Allgatherv(&nodeDoubleSizes[0],nodesOwned[BBSMPSMyPe],MPI_INT,&globalNodeDoubleSizes[0],&nodesOwned[0],&nodesOwnedDisplacement[0],MPI_INT,BBSMPSContext.comm());
	
	check3=(MPI_Wtime()-check3);
	double check4=(MPI_Wtime());
	 comm=MPI_Wtime()-comm;
/*	if(mype==0){
			for (int i=0; i<globalNodeDoubleSizes.size();i++){
				cout<<"["<<i<<","<<globalNodeIntSizes[i]<<","<<globalNodeDoubleSizes[i]<<"]";
			}
			cout<<endl;
	}
*/
	int greatIntBufferSize=0;
	int greatDblBufferSize=0;
	for (int i=0; i< totalNodesToExchange; i++){
		greatIntBufferSize+=globalNodeIntSizes[i];
		greatDblBufferSize+=globalNodeDoubleSizes[i];
	} 

	//Allocate great buffer
	vector<int> greatIntBuffer(greatIntBufferSize,0);
	vector<double> greatDblBuffer(greatDblBufferSize,0);

	//Serialize the nodes that belong to me, delete them from the tree

	
	vector<int> intRecvCounts (BBSMPSProcs,0);
	vector<int> dblRecvCounts (BBSMPSProcs,0);


	vector<int> AcumulatedGlobalNodeIntSizes (totalNodesToExchange,0);
	vector<int> AcumulatedGlobalNodeDoubleSizes (totalNodesToExchange,0);
	for (int i=1; i<AcumulatedGlobalNodeIntSizes.size();i++){
		AcumulatedGlobalNodeIntSizes[i]=AcumulatedGlobalNodeIntSizes[i-1]+globalNodeIntSizes[i-1];
		AcumulatedGlobalNodeDoubleSizes[i]=AcumulatedGlobalNodeDoubleSizes[i-1]+globalNodeDoubleSizes[i-1];
	//	if (BBSMPSMyPe==0) cout<<" accumulated int positions "<<i<< " is "<<AcumulatedGlobalNodeIntSizes[i-1]<<" + "<<globalNodeIntSizes[i-1]<<endl;
	}

	for (int nod=0;nod<totalNodesToExchange;nod++){

		if(globalBoundVector[nod].second.first==BBSMPSMyPe && globalBoundVector[nod].first<COIN_DBL_MAX){

			it=heap.begin();
			
		//	(*it)->printNode();
			int rankWithinOwner=globalBoundVector[nod].second.second;
			int intPtr=AcumulatedGlobalNodeIntSizes[nodesOwnedDisplacement[BBSMPSMyPe]+rankWithinOwner];
			int dblPtr=AcumulatedGlobalNodeDoubleSizes[nodesOwnedDisplacement[BBSMPSMyPe]+rankWithinOwner];
		//	cout<<"PACKING!"<<(*it)->getNodeNumber()<<" "<<intPtr<<" "<<dblPtr<<"!!";
			(*it)->serialize(&greatIntBuffer[intPtr],&greatDblBuffer[dblPtr]);
		//	cout<<" the tree size used to be "<<heap.size();
			heap.erase(heap.begin());
		//	cout<<" and it now is "<<heap.size()<<endl;
			
		}
	//	cout<<BBSMPSMyPe<<" "<<mype<<" "<<nod<<" summing "<<globalNodeIntSizes[nod]<<" "<<globalNodeDoubleSizes[nod]<<endl;
		
		

	}
	for (int proc=0; proc<BBSMPSProcs; proc++){
		for (int nod=0; nod< nodesOwned[proc]; nod++){
			intRecvCounts[proc]+=globalNodeIntSizes[nodesOwnedDisplacement[proc]+nod];
			dblRecvCounts[proc]+=globalNodeDoubleSizes[nodesOwnedDisplacement[proc]+nod];
		}
	}
	check4=(MPI_Wtime()-check4);
	double check5=(MPI_Wtime());
	
	vector<int> intDspls (BBSMPSProcs,0);
	vector<int> dblDspls (BBSMPSProcs,0);

	for(int i=1; i< BBSMPSProcs; i++){
		intDspls[i]=intDspls[i-1]+intRecvCounts[i-1];
		dblDspls[i]=dblDspls[i-1]+dblRecvCounts[i-1];
	}

	/*if (mype==0 ){
		cout<<"INT REcv COUNTS ";
		for(int i=0; i<intRecvCounts.size();i++){
			cout<<intRecvCounts[i]<<" ";

		}
		cout<<endl;
		cout<<"DBL REcv COUNTS ";
		for(int i=0; i<dblRecvCounts.size();i++){
			cout<<dblRecvCounts[i]<<" ";

		}
		cout<<endl;
		cout<<"int dsps COUNTS ";
		for(int i=0; i<intDspls.size();i++){
			cout<<intDspls[i]<<" ";

		}
		cout<<endl;

		cout<<"dbl dsps  counts "<<endl;
		for(int i=0; i<dblDspls.size();i++){
			cout<<dblDspls[i]<<" ";
			
		}
	}*/
 
	//cout<<" great buffer sizes "<<greatIntBuffer.size()<<" "<<greatDblBuffer.size()<<" "<<endl;
	//Allgatherv to the great buffer
	double comm1=MPI_Wtime();
	MPI_Allgatherv(MPI_IN_PLACE,0, MPI_DATATYPE_NULL,&greatIntBuffer[0],&intRecvCounts[0],&intDspls[0],MPI_INT,BBSMPSContext.comm());
	MPI_Allgatherv(MPI_IN_PLACE,0, MPI_DATATYPE_NULL,&greatDblBuffer[0],&dblRecvCounts[0],&dblDspls[0],MPI_DOUBLE,BBSMPSContext.comm());
	check5=(MPI_Wtime()-check5);
	double check6=(MPI_Wtime());
	comm1=MPI_Wtime()-comm1;
	/*if (mype==0 && BBSMPSMyPe==1){
		cout<<"THE GREAT VECTOR INT PART"<<endl;
		for(int i=0; i<greatIntBuffer.size();i++){
			cout<<greatIntBuffer[i]<<" ";

		}
		cout<<"THE GREAT VECTOR DBL PART"<<endl;
		for(int i=0; i<greatDblBuffer.size();i++){
			cout<<greatDblBuffer[i]<<" ";
			
		}
	}*/

		

	//	cout<<" proc "<<mype<<" "<<BBSMPSMyPe<<" got here!!!!!!!!!!"<<endl;
	//For every node that I own, add to the tree.
	if(doINeedComm){


		int myStart=0;
		for (int i=0; i<BBSMPSMyPe;i++) if (abs(allLBsGathered[i]-bestLB)>=commTol)myStart++;
		//cout<<"MYSTART " <<myStart<<endl;
		for (int nod=myStart;nod<totalNodesToExchange;nod+=nProcsThatNeedComm){
			if (globalBoundVector[nod].first<COIN_DBL_MAX){
				int procOwner=globalBoundVector[nod].second.first;
				int rankWithinOwner=globalBoundVector[nod].second.second;
			//	cout<<"trying to gather node no "<<rankWithinOwner<<" of owner "<<procOwner<<" the address that we are trying to access is "<<nodesOwnedDisplacement[procOwner]+rankWithinOwner<<endl;
				int intPtr=AcumulatedGlobalNodeIntSizes[nodesOwnedDisplacement[procOwner]+rankWithinOwner];
				int dblPtr=AcumulatedGlobalNodeDoubleSizes[nodesOwnedDisplacement[procOwner]+rankWithinOwner];
		//		cout<<" proc "<<mype<<" "<<BBSMPSMyPe<<" init a node from "<<intPtr<<" "<<dblPtr<<" INDEX TRYING TO ACCESS "<<nodesOwnedDisplacement[procOwner]+rankWithinOwner<<endl;
			//	cout<<"UNPACKING!!"<<intPtr<<" "<<dblPtr<<"!!";
				BBSMPSNode *node = new BBSMPSNode(&greatIntBuffer[intPtr],&greatDblBuffer[dblPtr]);	
				
				//cout<<" got the obj of th enide LP to value "<<node->getObjective()<<" the node at the top is node "<<node->getNodeNumber()<<endl;
		//		(node)->printNode();
				heap.insert(node);
				//delete node;
			}
		}
	}

	int heapSizeHa=heap.size();
	
	MPI_Allreduce(MPI_IN_PLACE,&heapSizeHa, 1, MPI_INT,MPI_MAX,ctx.comm());
	
	if (heapSizeHa!=heap.size() ) {
		BBSMPS_ALG_LOG_SEV(warning)<<"ATENCIO!!!! DIFFERENCENT TREE SIZES!!! "<<heapSizeHa<<" "<<heap.size();
	}

	//cout<<" the size of the heap is "<<heap.size()<<endl;
	BBSMPSNode *currentNode_ptr=*(heap.begin());
	//Update LB/UB if necessary
	objLB=currentNode_ptr->getObjective();
	//cout<<" set the obj LP to value "<<objLB<<" the node at the top is node "<<currentNode_ptr->getNodeNumber()<<endl;
	check6=(MPI_Wtime()-check6);
	cout<<"TIMES ch1:"<<check1<<" "<<check2<<" "<<check3<<" "<<check4<<" "<<check5<<" "<<check6<<" "<<comm<<" "<<comm1<<endl;
	totalCommTime+=(MPI_Wtime()-startCheckpoint);
    totalTimesCommCalled++;
}


bool BBSMPSTree::shouldWePerformCommunication(int &counterToLastIter, bool &commDone){

	double startCheckpoint=MPI_Wtime();
	BAContext BBSMPSContext=BBSMPSSolver::instance()->getSBBContext();
	int BBSMPSProcs=BBSMPSContext.nprocs();
	int BBSMPSMyPe=BBSMPSSolver::instance()->getSBBMype();
	int mype=BBSMPSSolver::instance()->getMype();
	BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
//	BBSMPS_ALG_LOG_SEV(warning)<<"processor "<<BBSMPSMyPe<<" "<<mype<<" got into should we run!!";
	MPI_Status status;
	if (communicationPending){

		if (!commDone){//taking care of the first iteration
			MPI_Wait(request0,&status);

			MPI_Wait(request1,&status);

			MPI_Wait(request2,&status);

			MPI_Wait(request3,&status);

			MPI_Wait(request4,&status);
		}

		/*if(counterToLastIter%(iterationsBetweenCommunication*4)==(iterationsBetweenCommunication*4-1)){
			cout<<"PROCESSOR "<<BBSMPSMyPe<<" got into wait mode at iteration "<<bbIterationCounter<<" ct to last iter "<<counterToLastIter<<endl;
			double startCheckpoint2=MPI_Wtime();
	
			MPI_Wait(request0,&status);

			MPI_Wait(request1,&status);

			MPI_Wait(request2,&status);

			MPI_Wait(request3,&status);

			MPI_Wait(request4,&status);
			totalIdleTime+=(MPI_Wtime()-startCheckpoint2);
		}*/

		int flag=0;
		int flagAux=0;
		
		MPI_Test(request0,&flagAux,&status);
		flag+=flagAux;
		MPI_Test(request1,&flagAux,&status);
		flag+=flagAux;
		MPI_Test(request2,&flagAux,&status);
		flag+=flagAux;
		MPI_Test(request3,&flagAux,&status);
		flag+=flagAux;
		MPI_Test(request4,&flagAux,&status);
		flag+=flagAux;
	//	BBSMPS_ALG_LOG_SEV(warning)<<"and the sum is "<<flag;
		if (flag==5){
			int weFoundMPIReq=1;
			MPI_Allreduce(MPI_IN_PLACE,&weFoundMPIReq, 1, MPI_INT,MPI_MAX,ctx.comm());
		}
		else{
			int weFoundMPIReq=0;
			MPI_Allreduce(MPI_IN_PLACE,&weFoundMPIReq, 1, MPI_INT,MPI_MAX,ctx.comm());
			if (weFoundMPIReq==1){
				double startCheckpoint2=MPI_Wtime();
	
				MPI_Wait(request0,&status);

				MPI_Wait(request1,&status);

				MPI_Wait(request2,&status);

				MPI_Wait(request3,&status);

				MPI_Wait(request4,&status);
				totalIdleTime+=(MPI_Wtime()-startCheckpoint2);
				flag=5;
			}
		}
		if (flag==5){
	//		BBSMPS_ALG_LOG_SEV(warning)<<"ENTERING ON FIRST COMMUNICATION";
			communicationPending=false;

			objUB=min(objUB,bufferedUB);
		//	BBSMPS_ALG_LOG_SEV(warning)<<"we got out of communication "<<bufferedUB<<" "<<bufferedItCounter<<" "<<bufferedNodesLeft<<" "<<bufferedBestLB<<" "<<bufferedWorstLB;
			if ((BBSMPSSolver::instance()->getWallTime())>tiLim){

				if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Time Limit reached.";
				totalCommTime+=(MPI_Wtime()-startCheckpoint);
				inTermination=true;
				
				return false;
			}
			if (BBSMPSSolver::instance()->getSolPoolSize()-solsDiscoveredInit>=solsDiscoveredLimit){
				if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Solution Limit reached.";
				inTermination=true;
				totalCommTime+=(MPI_Wtime()-startCheckpoint);
				return false;
			}

			if (bufferedItCounter>nodeLim){
				if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Node Limit reached.";
				inTermination=true;
				totalCommTime+=(MPI_Wtime()-startCheckpoint);
				return false;
			}

			if (bufferedNodesLeft==0){
				if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Heap is empty.";
				inTermination=true;
				totalCommTime+=(MPI_Wtime()-startCheckpoint);
				return false;
			}
			counterToLastIter=0;
			double generalGap = fabs(objUB-bufferedBestLB)*100/(fabs(objUB)+10e-10);
	//		BBSMPS_ALG_LOG_SEV(warning)<<" COMPARISON "<<bufferedBestLB-bufferedWorstLB<<" "<<(abs(bufferedBestLB-bufferedWorstLB)>=commTol)<<" "<<((bufferedWorstLB - compTol) >= objUB)<<" "<<(generalGap<0.01);
			if ((bufferedNodesLeft>=BBSMPSProcs) && (abs(bufferedBestLB-bufferedWorstLB)>=commTol || (bufferedWorstLB - compTol) >= objUB  || generalGap<0.05)) {//communication needed
		//		BBSMPS_ALG_LOG_SEV(warning)<<"PROCESSOR "<<BBSMPSMyPe<<" called communication at iteration "<<bbIterationCounter<<" ct to last iter "<<counterToLastIter;
				totalTimesCommCalled++;
				totalCommTime+=(MPI_Wtime()-startCheckpoint);
				iterationsBetweenCommunication/=1.5;
				iterationsBetweenCommunication=max(iterationsBetweenCommunication,MIN_COMM_ITERS);
				if (generalGap<0.05)iterationsBetweenCommunication=RAMPDOWN_COMM_ITERS;
				cout<<"ITERATIONS SET TO "<<iterationsBetweenCommunication<<" "<<bufferedBestLB-bufferedWorstLB<<" "<<((bufferedWorstLB - compTol) >= objUB)<<" "<<(bufferedWorstLB - compTol)<<" "<<objUB<<" "<<generalGap<<endl;
				return true;
			}
			else{
				
				iterationsBetweenCommunication*=1.5;
				iterationsBetweenCommunication=min(iterationsBetweenCommunication,MAX_COMM_ITERS);
				if (generalGap<0.05)iterationsBetweenCommunication=RAMPDOWN_COMM_ITERS;
				cout<<"ITERATIONS SET TO "<<iterationsBetweenCommunication<<" "<<bufferedBestLB-bufferedWorstLB<<" "<<((bufferedWorstLB - compTol) >= objUB)<<" "<<(bufferedWorstLB - compTol)<<" "<<objUB<<" "<<generalGap<<endl;
			}
		}

		totalCommTime+=(MPI_Wtime()-startCheckpoint);
		return false;
	}
//	cout<<"processor is here "<<counterToLastIter<<" "<<commDone<<" "<<(counterToLastIter%100==0 && commDone)<<" "<<(heap.size()==0 || (!commDone && heap.size()>16))<<"FINAL RESULT "<<!((counterToLastIter%10==0 && commDone) || heap.size()==0 || (!commDone && heap.size()>16))<<endl;
	if ( !((counterToLastIter%iterationsBetweenCommunication==0 && commDone) || heap.size()==0 || (!commDone && heap.size()>BBSMPSProcs))){
	//	cout<<"processor "<<BBSMPSMyPe<<" we are not running this turn"<<endl;
		totalCommTime+=(MPI_Wtime()-startCheckpoint);
		return false;
	}
		
	
	bufferedUB=objUB;
	bufferedItCounter=bbIterationCounter;
	bufferedWorstLB=objLB;
	bufferedBestLB=objLB;
	bufferedNodesLeft=heap.size();
//	cout<<"PROCESSOR "<<BBSMPSMyPe<<" LAUNCHED ALLREDUCE communication at iteration "<<bbIterationCounter<<" ct to last iter "<<counterToLastIter<<endl;
	//cout<<"setting up a new allreduce "<<bufferedUB<<" "<<bufferedItCounter<<" "<<bufferedNodesLeft<<" "<<bufferedBestLB<<" "<<bufferedWorstLB<<endl;
    MPI_Iallreduce(MPI_IN_PLACE,&bufferedUB, 1, MPI_DOUBLE,MPI_MIN,BBSMPSContext.comm(),request0);
    MPI_Iallreduce(MPI_IN_PLACE,&bufferedItCounter, 1, MPI_INT,MPI_SUM,BBSMPSContext.comm(),request1);
    MPI_Iallreduce(MPI_IN_PLACE,&bufferedWorstLB, 1, MPI_DOUBLE,MPI_MAX,BBSMPSContext.comm(),request2);
	MPI_Iallreduce(MPI_IN_PLACE,&bufferedBestLB, 1, MPI_DOUBLE,MPI_MIN,BBSMPSContext.comm(),request3);
	MPI_Iallreduce(MPI_IN_PLACE,&bufferedNodesLeft, 1, MPI_INT,MPI_SUM,BBSMPSContext.comm(),request4);

	
	
	communicationPending=true;
	

	totalCommTime+=(MPI_Wtime()-startCheckpoint);
	return false;

	
}

