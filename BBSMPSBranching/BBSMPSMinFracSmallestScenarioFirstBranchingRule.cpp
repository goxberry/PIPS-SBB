#include "BBSMPSMinFracSmallestScenarioFirstBranchingRule.hpp"

using namespace std;
int BBSMPSMinFracSmallestScenarioFirstBranchingRule::getFirstStageMinIntInfeasCol( const denseBAVector& primalSoln,  SMPSInput& input) {
	int col;

    // Return first index of integer variable with fractional value
	for (col = 0; col < input.nFirstStageVars(); col++)
	{

		bool isColInteger = input.isFirstStageColInteger(col);
		bool isValInteger = isIntFeas(primalSoln.getFirstStageVec()[col], intTol); //TODO:: Tolerance is hardcoded for now

		// If the (col)th 1st stage primal variable is integer,
		// but has a fractional value, return idx
		if(isColInteger && !isValInteger) return col;
	}

    // Otherwise, 1st stage is integer feasible: return -1;
	return -1;
}


int BBSMPSMinFracSmallestScenarioFirstBranchingRule::getFirstStageMinFracPartCol( const denseBAVector& primalSoln,  SMPSInput& input){
	int col;

	double minFracPart = COIN_DBL_MAX;
	int minFracPartCol = -1;

    // Return index of integer variable with largest fractional part
	for (col = 0; col < input.nFirstStageVars(); col++)
	{
		bool isColInteger = input.isFirstStageColInteger(col);
		double colFracPart = fracPart(primalSoln.getFirstStageVec()[col]);

		if(isColInteger && !isIntFeas(colFracPart, intTol) && (colFracPart < minFracPart) ) {
			minFracPartCol = col;
			minFracPart = colFracPart;
		}
	}

	
	return minFracPartCol;

}


void BBSMPSMinFracSmallestScenarioFirstBranchingRule::branchOnFirstStage(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes, const denseBAVector& primalSoln,  SMPSInput& input) {

	int mype=BBSMPSSolver::instance()->getMype();
    /* Branching Rule */
    // For now, get minimal index of an integer infeasible variable
    //int branchCol = getFirstStageMinIntInfeasCol(primalSoln);

    // Get index of maximum fractional part.
	int branchCol = getFirstStageMinFracPartCol(primalSoln,input);
    assert(branchCol > -1); // Should always be true if not integer feasible

    if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Branching on first stage variable "<< branchCol <<".";

	//Create both branching infos
    std::vector<BBSMPSBranchingInfo> bInfosLeftKid;
    std::vector<BBSMPSBranchingInfo> bInfosRightKid;
    bInfosLeftKid.push_back( BBSMPSBranchingInfo(branchCol, ceil(primalSoln.getFirstStageVec()[branchCol]), 'L', 1));
    bInfosRightKid.push_back( BBSMPSBranchingInfo(branchCol, floor(primalSoln.getFirstStageVec()[branchCol]), 'U', 1));
	//Create both children

    BBSMPSNode *leftKidNode= new BBSMPSNode(node, bInfosLeftKid,BBSMPSSolver::instance()->getSBBMype());
    BBSMPSNode *rightKidNode= new BBSMPSNode(node, bInfosRightKid,BBSMPSSolver::instance()->getSBBMype());
    childNodes.push_back(leftKidNode);
    childNodes.push_back(rightKidNode);

}

int BBSMPSMinFracSmallestScenarioFirstBranchingRule::getSecondStageMinFracPartCol( const denseBAVector& primalSoln,  SMPSInput& input, int scen){
	int col;

	double minFracPart = COIN_DBL_MAX;
	int minFracPartCol = -1;

    // Return index of integer variable with largest fractional part
	for (col = 0; col < input.nSecondStageVars(scen); col++)
	{
		bool isColInteger = input.isSecondStageColInteger(scen, col);
		double colFracPart = fracPart(primalSoln.getSecondStageVec(scen)[col]);

		if(isColInteger  && !isIntFeas(colFracPart, intTol)&& (colFracPart < minFracPart) ) {
			minFracPartCol = col;
			minFracPart = colFracPart;
		}
	}

	
	return minFracPartCol;

}


int BBSMPSMinFracSmallestScenarioFirstBranchingRule::getSecondStageMinIntInfeasCol( const denseBAVector& primalSoln, int scen,   SMPSInput& input) {
	
	int col;
	for (col = 0; col < input.nSecondStageVars(scen); col++)
	{
		bool isColInteger = input.isSecondStageColInteger(scen, col);
		bool isValInteger = isIntFeas(primalSoln.getSecondStageVec(scen)[col], intTol);

		// If the (col)th 2nd stage primal variable of the (scen)th	
		// scenario is integer, but has fractional value, return idx
		if (isColInteger && !isValInteger) return col;
	}

    // Otherwise, return -1;
	return -1;
}


int BBSMPSMinFracSmallestScenarioFirstBranchingRule::isSecondStageIntFeas( const denseBAVector& primalSoln, int scen, SMPSInput & input) {
	return (getSecondStageMinIntInfeasCol(primalSoln, scen, input) == -1);
}

int BBSMPSMinFracSmallestScenarioFirstBranchingRule::getNSecondStageFracVars( const denseBAVector& primalSoln,  SMPSInput& input, int scen){
	int col;

	int counter=0;

    // Return index of integer variable with largest fractional part
	for (col = 0; col < input.nSecondStageVars(scen); col++)
	{
		bool isColInteger = input.isSecondStageColInteger(scen, col);
		bool isValInteger = isIntFeas(primalSoln.getSecondStageVec(scen)[col], intTol);

		counter+=(isColInteger && !isValInteger);
	}

	return counter;

}


void BBSMPSMinFracSmallestScenarioFirstBranchingRule::branchOnSecondStage(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes,  const denseBAVector& primalSoln,  SMPSInput& input,BAContext &ctx, int mype) {


	//Find out who has the best scenario
	int bestLocalScen=-1;
	double bestLocalScenarioSize=COIN_DBL_MAX;

	for (int scen = 0; scen < input.nScenarios(); scen++)
	{
		if(ctx.assignedScenario(scen)) {
			int scenFracVars=getNSecondStageFracVars (primalSoln,input, scen);
			if (scenFracVars!=0 && bestLocalScenarioSize> scenFracVars){
				bestLocalScen=scen;
				bestLocalScenarioSize=scenFracVars;
			}
		}
	}

	comp localComp;
    comp globalComp;

    localComp.result=bestLocalScenarioSize;
    localComp.rank=mype;
    MPI_Allreduce (&localComp, &globalComp, 1, MPI_DOUBLE_INT, MPI_MINLOC, ctx.comm());


    // Then, for that scenario number, get the minimal index of
    // an integer infeasible decision variable, and branch on that column
	std::vector<BBSMPSBranchingInfo> bInfosLeftKid;
	std::vector<BBSMPSBranchingInfo> bInfosRightKid;

	if(globalComp.rank==mype) {
		int localBestCol=getSecondStageMinFracPartCol( primalSoln,input, bestLocalScen);
		
		BBSMPS_ALG_LOG_SEV(info) << "Processor " << mype << " will branch on variable "<<localBestCol<<" of second stage scenario "
	<< bestLocalScen << " "<<primalSoln.getSecondStageVec(bestLocalScen)[localBestCol]<<".";

		bInfosLeftKid.push_back( BBSMPSBranchingInfo(localBestCol, ceil(primalSoln.getSecondStageVec(bestLocalScen)[localBestCol]), 'L', 2,bestLocalScen));
		bInfosRightKid.push_back( BBSMPSBranchingInfo(localBestCol, floor(primalSoln.getSecondStageVec(bestLocalScen)[localBestCol]), 'U', 2,bestLocalScen));
	}



	BBSMPSNode *leftKidNode= new BBSMPSNode(node, bInfosLeftKid,BBSMPSSolver::instance()->getSBBMype());
	BBSMPSNode *rightKidNode= new BBSMPSNode(node, bInfosRightKid,BBSMPSSolver::instance()->getSBBMype());
	childNodes.push_back(leftKidNode);
	childNodes.push_back(rightKidNode);

}
	// NOTE: MPI standard requires passing ints, not bools
int BBSMPSMinFracSmallestScenarioFirstBranchingRule::isFirstStageIntFeas( const denseBAVector& primalSoln,  SMPSInput& input) {
	return (getFirstStageMinIntInfeasCol(primalSoln,input) == -1);
}

bool BBSMPSMinFracSmallestScenarioFirstBranchingRule::branch(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes,  const denseBAVector& primalSoln){
	/* Branching */
	// Decide which stage to branch on:
	// If first stage decision variables not integer feasible,
	// branch on a first stage variable, go to start of loop
	timesCalled++;
	SMPSInput &input= BBSMPSSolver::instance()->getSMPSInput();
	BAContext &ctx= BBSMPSSolver::instance()->getBAContext();
	int mype=BBSMPSSolver::instance()->getMype();

	if(!isFirstStageIntFeas(primalSoln,input)) {
		branchOnFirstStage(node, childNodes, primalSoln,input);
		timesSuccessful++;
		return true;
	}

    // If we get to this point, we know that the first stage decision variables
    // are integer feasible, but the LP solution is not integer feasible, so
    // one of the second stage scenarios must not be integer feasible, and
    // one of the variables in one of those scenarios should be branched on.
	branchOnSecondStage(node,childNodes,primalSoln,input,ctx,mype);
	timesSuccessful+=(childNodes.size()>0);
	return (childNodes.size()>0);

}