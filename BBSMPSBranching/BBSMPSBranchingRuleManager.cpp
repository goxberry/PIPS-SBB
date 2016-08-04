#include "BBSMPSBranchingRuleManager.hpp"

BBSMPSBranchingRuleManager::BBSMPSBranchingRuleManager(){};
BBSMPSBranchingRuleManager::~BBSMPSBranchingRuleManager(){};


void BBSMPSBranchingRuleManager::addBranchingRule(BBSMPSBranchingRule *rule){
	assert(rule!=NULL);
	branchingRuleList.insert(rule);
}

bool BBSMPSBranchingRuleManager::branch(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes, const denseBAVector& primalSoln){
	bool success=false;
	std::multiset<BBSMPSBranchingRule*>::iterator it;
	for (it=branchingRuleList.begin(); it!=branchingRuleList.end() && !success; ++it){
		BBSMPSBranchingRule *br=(*it);
		double ch1=MPI_Wtime();
		success= success || br->branch(node, childNodes,  primalSoln);
		br->setAccumulatedRunTime(br->getAccumulatedRunTime()+(MPI_Wtime()-ch1));
	}
	return success;
}


void BBSMPSBranchingRuleManager::printStatistics(){
	BBSMPS_ALG_LOG_SEV(warning)<<"**************BRANCHING STATISTICS****************";
	std::multiset<BBSMPSBranchingRule*>::iterator it;
	for (it=branchingRuleList.begin(); it!=branchingRuleList.end(); ++it){
		BBSMPSBranchingRule *br=(*it);
		br->printStatistics();
	}
	BBSMPS_ALG_LOG_SEV(warning)<<"**************************************************";
}

void BBSMPSBranchingRuleManager::freeResources(){
	std::multiset<BBSMPSBranchingRule*>::iterator it;
	for (it=branchingRuleList.begin(); it!=branchingRuleList.end(); ++it){
		BBSMPSBranchingRule *br=(*it);
		delete br;
	}
	branchingRuleList.clear();
}