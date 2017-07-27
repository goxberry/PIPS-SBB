#include "BBSMPSSolver.hpp"
using namespace std;

BBSMPSSolver *BBSMPSSolver::solverInstance = 0;


BAContext& BBSMPSSolver::getBAContext(){
  return pipsCtx;
}

int& BBSMPSSolver::getMype(){
  return mype;
}

 int& BBSMPSSolver::getSBBMype(){
  return sBBMype;
 }

BAContext& BBSMPSSolver::getSBBContext(){
  return pipsSbbCtx;
}

SMPSInput& BBSMPSSolver::getSMPSInput(){
  return input;
}

PIPSSInterface& BBSMPSSolver::getPIPSInterface(){
  return (*rootSolver);
}

BADimensions& BBSMPSSolver::getBADimensions(){
  return dims; 
}

const BADimensionsSlacks& BBSMPSSolver::getBADimensionsSlacks(){
  return (*rootSolver).getSlackDims();
}

const BADimensionsSlacks& BBSMPSSolver::getOriginalBADimensionsSlacks(){
  return originalDims;
}

const denseBAVector& BBSMPSSolver::getOriginalLB(){
  return lbModifiableWithCuts;
}

const denseBAVector& BBSMPSSolver::getOriginalUB(){
  return ubModifiableWithCuts;
}

BBSMPSSolver *BBSMPSSolver::instance(){
  return solverInstance;
}

void BBSMPSSolver::setLPRelaxationObjectiveValue(double lpRelObjVal){
  LPRelaxationObjectiveValue=lpRelObjVal;
}

double BBSMPSSolver::getLPRelaxationObjectiveValue(){
  return LPRelaxationObjectiveValue;
}


//TODO: The relaxation could not be initialized. we need proper guards.
denseBAVector& BBSMPSSolver::getLPRelaxation(){
  return LPRelaxation;
}
void BBSMPSSolver::setLPRelaxation(denseBAVector &_LPRelaxation){
  LPRelaxation=_LPRelaxation;
}


BBSMPSSolver *BBSMPSSolver::initialize(const SMPSInput &_input, int nSolvers){
 if (solverInstance) delete solverInstance;
 //Work the details of the contexts here
//Initialize global context
//  cout<<"initting, solver size is "<<nSolvers<<endl;
 BAContext generalCtx(MPI_COMM_WORLD);
 int generalMype=generalCtx.mype();
 int generalSize=generalCtx.nprocs();
 MPI_Group generalGroup;
 MPI_Comm_group(MPI_COMM_WORLD, &generalGroup);

 //Initialize PIPSCONTEXT
 int pipsContextSize=generalSize/nSolvers;
 int pipsRanks[pipsContextSize];
 int myPipsCtxStart=generalMype/pipsContextSize;
 for(int i=0; i<pipsContextSize;i++) pipsRanks[i]=myPipsCtxStart*pipsContextSize+i;
 //cout<<"MYPE "<<generalMype<< " ctxsize "<<pipsContextSize<<" myPipsCtxStart "<<myPipsCtxStart<<" "<<pipsRanks[0]<<" "<<pipsRanks[1]<<" "<<pipsRanks[2]<<" "<<pipsRanks[3]<<endl;
 MPI_Group pipsGroup;
 MPI_Group_incl(generalGroup, pipsContextSize, pipsRanks, &pipsGroup);

 MPI_Comm pipsComm;
 MPI_Comm_create(MPI_COMM_WORLD, pipsGroup, &pipsComm);

 //Initialize PIPSSBBCONTEXT
 int pipsSbbContextSize=nSolvers;
 int pipsSbbRanks[pipsSbbContextSize];
 int myPipsSbbCtxStart=generalMype%pipsContextSize;

 for(int i=0; i<pipsSbbContextSize;i++) pipsSbbRanks[i]=myPipsSbbCtxStart+i*pipsContextSize;
 MPI_Group pipsSbbGroup;
MPI_Group_incl(generalGroup, pipsSbbContextSize, pipsSbbRanks, &pipsSbbGroup);
 MPI_Comm pipsSbbComm;
 //cout<<"MYPE "<<generalMype<< " ctxsize "<<pipsSbbContextSize<<" myPipsbbCtxStart "<<myPipsSbbCtxStart<<" "<<pipsSbbRanks[0]<<" "<<pipsSbbRanks[1]<<" "<<pipsSbbRanks[2]<<endl;
 MPI_Comm_create(MPI_COMM_WORLD, pipsSbbGroup, &pipsSbbComm);

 int pipsMype;
 int pipsSize;
 MPI_Comm_rank ( pipsComm, &pipsMype );
 MPI_Comm_size ( pipsComm, &pipsSize );
 int pipssbbMype;
 int pipssbbSize;
 MPI_Comm_rank ( pipsSbbComm, &pipssbbMype );
 MPI_Comm_size ( pipsSbbComm, &pipssbbSize );
// cout<<"MYPE "<<generalMype<< " pips "<<pipsMype<<" "<<pipsSize<<" "<<pipssbbMype<<" "<<pipssbbSize<<endl;

 //Initialize the solver with the three contexts
 solverInstance = new BBSMPSSolver(_input,pipsComm,pipsSbbComm);
 return solverInstance;
}

BBSMPSSolver::BBSMPSSolver(const SMPSInput &_input,MPI_Comm &pipsComm, MPI_Comm &pipsSbbComm):
ctx(MPI_COMM_WORLD),
pipsCtx(pipsComm),
pipsSbbCtx(pipsSbbComm),
mype(pipsCtx.mype()),
sBBMype(pipsSbbCtx.mype()),
input(_input),
problemData(input,pipsCtx),
//pre(problemData,input),
dims(problemData.dims.inner),
dimsSlacks(dims),
startTimeStamp(MPI_Wtime()){
 // cout<<"we went through initialization"<<endl;
  rootSolver= new PIPSSInterface(problemData, PIPSSInterface::useDual);
  originalDims=BADimensionsSlacks((*rootSolver).getSlackDims());
  lb = (*rootSolver).getLB();
  ub = (*rootSolver).getUB();
   lbModifiableWithCuts = (*rootSolver).getLB();
  ubModifiableWithCuts = (*rootSolver).getUB();
  

}

BBSMPSSolver::~BBSMPSSolver(){
  if (rootSolver) delete rootSolver;
}



double BBSMPSSolver::getWallTime(){
  return MPI_Wtime()-startTimeStamp;
}

bool BBSMPSSolver::isInitialized(){
  return (solverInstance!=NULL);
}

void BBSMPSSolver::deInitialize(){
  if(isInitialized()){
    delete solverInstance;
  }
}

const BAFlagVector<variableState>& BBSMPSSolver::getOriginalWarmStart(){

  return originalWarmStartModifiableWithCuts;
}
void BBSMPSSolver::setOriginalWarmStart(BAFlagVector<variableState>&warmStart){
  originalWarmStart=warmStart;
  originalWarmStartModifiableWithCuts=warmStart;
}

void BBSMPSSolver::printPresolveStatistics(){


  BBSMPS_ALG_LOG_SEV(warning)<<"~~~~~~~~~~~~~~~PRESOLVE STATISTICS~~~~~~~~~~~~~~~~";
   
  //  BBSMPS_ALG_LOG_SEV(warning)<<"Number Of Presolves:"<<pre.numPresolves<<":Number Of Bound Chgs:"<<pre.numBdsChg;
  //BBSMPS_ALG_LOG_SEV(warning)<<"Number Of RHS Chgs:"<<pre.numRhsChg<<":Number Of Coefficient Chgs:"<<pre.numCoeffChg;
   

BBSMPS_ALG_LOG_SEV(warning)<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
}


void BBSMPSSolver::addSolutionToPool(BBSMPSSolution &sol){

  double objUB=COIN_DBL_MAX;
  if (getSolPoolSize()>0)objUB=getSoln(0).getObjValue();
  if (objUB>sol.getObjValue()){
    denseBAVector solVector;
    sol.getSolutionVector(solVector);
    denseBAVector originalSpaceSolVector;
    originalSpaceSolVector.allocate(getOriginalBADimensionsSlacks(), getBAContext(), PrimalVector); 
    originalSpaceSolVector.copyAndShrinkToDims(solVector);

    BBSMPSSolution originalSpaceSolution(originalSpaceSolVector,sol.getObjValue(),sol.getTimeOfDiscovery());

    solutionPool.insert(originalSpaceSolution);
  }
};


void BBSMPSSolver::printSolutionStatistics(double objLB){


    double bestTime=COIN_DBL_MAX;

    for (std::multiset<BBSMPSSolution,solutionComparison>::iterator it=solutionPool.begin(); it!=solutionPool.end(); ++it){
      BBSMPSSolution s = *it;
      if (s.getTimeOfDiscovery()<bestTime)bestTime=s.getTimeOfDiscovery();
    }

      BBSMPS_ALG_LOG_SEV(warning)<<"---------------SOLUTION STATISTICS----------------";
      
      BBSMPS_ALG_LOG_SEV(warning)<<"Solution Pool Size:"<<solutionPool.size()<<":Time To First Solution:"<<bestTime;
  int itCounter=0;
  for (std::multiset<BBSMPSSolution,solutionComparison>::iterator it=solutionPool.begin(); it!=solutionPool.end(); ++it){
      BBSMPSSolution s = *it;
      double solGap = fabs(s.getObjValue()-objLB)*100/(fabs(s.getObjValue())+10e-10);
      BBSMPS_ALG_LOG_SEV(warning)<<"Solution:"<<itCounter<<":Solution Value:"<<s.getObjValue()<<":Time Of Discovery:"<<s.getTimeOfDiscovery()<<":Solution Gap:"<<solGap;
    itCounter++;
    }
  
  BBSMPS_ALG_LOG_SEV(warning)<<"--------------------------------------------------";
}

const BBSMPSSolution &BBSMPSSolver::getSoln(int index){
  assert(index<solutionPool.size());

  std::multiset<BBSMPSSolution,solutionComparison>::iterator it=solutionPool.begin();
  int ctr=0;
  while (ctr<index && ctr< solutionPool.size()){
    it++;
    ctr++;
  }
  return (*it);

}
const BBSMPSSolution &BBSMPSSolver::getSolnBySolNumber(int number){
  
  std::multiset<BBSMPSSolution,solutionComparison>::iterator it=solutionPool.begin();
  int ctr=0;
  while ( ctr< solutionPool.size()){
    if ((*it).getSolNumber()==number)return (*it);
    it++;
    ctr++;
  }
  assert(false);
  return (*solutionPool.begin());

}



int BBSMPSSolver::getSolPoolSize(){
  return solutionPool.size();
}

void BBSMPSSolver::resetSolver(){
   PIPSSInterface *rootS2 = new PIPSSInterface(problemData, PIPSSInterface::useDual);
   delete rootSolver;
   rootSolver=rootS2;
   originalDims=BADimensionsSlacks((*rootSolver).getSlackDims());
  lb = (*rootSolver).getLB();
  ub = (*rootSolver).getUB();
   lbModifiableWithCuts = (*rootSolver).getLB();
  ubModifiableWithCuts = (*rootSolver).getUB();

  BAFlagVector<variableState> newWS;
  newWS.allocate(originalDims, pipsCtx, PrimalVector); 
  newWS.copyFrom(originalWarmStartModifiableWithCuts);
  originalWarmStartModifiableWithCuts=newWS;


}

void BBSMPSSolver::commitNewColsAndRows(){

  // cout<<"before starting "<<(*rootSolver).getLB().getFirstStageVec().length()<<endl;


    (*rootSolver).commitNewColsAndRows();
   
    const BADimensionsSlacks& currentDims =getBADimensionsSlacks();
    BAFlagVector<variableState> newWS;
    newWS.allocate(currentDims, pipsCtx, PrimalVector); 

    int oldVars1= originalDims.inner.numFirstStageVars();
    int newVars1= currentDims.inner.numFirstStageVars();
    int oldCons1= originalDims.inner.numFirstStageCons();
    int newCons1= currentDims.inner.numFirstStageCons();
   // std::cout<<"settings "<<oldVars1<<" "<<oldCons1<<" new "<<newVars1<<" "<<newCons1<<" "<<std::endl;
    //if (oldVars1==newVars1 && oldCons1==newCons1)
     // (*rootSolver).d.diffAttributes(problemData);
    denseFlagVector<variableState> &stageVec = newWS.getFirstStageVec();
    const denseFlagVector<variableState> &stageVec2 = originalWarmStart.getFirstStageVec();
    for (int j = 0; j < oldVars1; j++) stageVec[j]=stageVec2[j];
    for (int j = oldVars1; j < newVars1; j++) stageVec[j]=AtLower;
    for (int j = 0; j < oldCons1; j++) stageVec[newVars1+j]=stageVec2[oldVars1+j];
    for (int j = oldCons1; j < newCons1; j++) {
      //std::cout<<"setting basis 1 at position "<<newVars1+j<<std::endl;
      stageVec[newVars1+j]=Basic;
    }

    for (int scen = 0; scen < input.nScenarios(); scen++) {
      if(pipsCtx.assignedScenario(scen)) {

        denseFlagVector<variableState> &stageVec = newWS.getSecondStageVec(scen);
        const denseFlagVector<variableState> &stageVec2 = originalWarmStart.getSecondStageVec(scen);
        int oldVars2= originalDims.inner.numSecondStageVars(scen);
        int newVars2= currentDims.inner.numSecondStageVars(scen);
        int oldCons2= originalDims.inner.numSecondStageCons(scen);
        int newCons2= currentDims.inner.numSecondStageCons(scen);
     //   std::cout<<"settings "<<oldVars2<<" "<<oldCons2<<" new "<<newVars2<<" "<<newCons2<<" "<<std::endl;
        for (int j = 0; j < oldVars2; j++) stageVec[j]=stageVec2[j];
        for (int j = oldVars2; j < newVars2; j++) stageVec[j]=AtLower;
        for (int j = 0; j < oldCons2; j++) stageVec[newVars2+j]=stageVec2[oldVars2+j];
        for (int j = oldCons2; j < newCons2; j++) {
          stageVec[newVars2+j]=Basic;
      //   std::cout<<"setting basis 2,"<<scen<<" at position "<<newVars2+j<<std::endl;
        }
      }
    }

    originalWarmStartModifiableWithCuts=newWS;
    
    
    //lbModifiableWithCuts.deallocate();
    //ubModifiableWithCuts.deallocate();

    lbModifiableWithCuts.allocate(currentDims,pipsCtx,PrimalVector);
    ubModifiableWithCuts.allocate(currentDims,pipsCtx,PrimalVector);
    lbModifiableWithCuts = (*rootSolver).getLB();
    ubModifiableWithCuts = (*rootSolver).getUB();

    denseVector &vLB1 = lbModifiableWithCuts.getFirstStageVec();
    const denseVector &vLB2 = lb.getFirstStageVec();
    //cout<<"we just allocated bounds again of size "<<vLB1.length()<<" while original has size of "<<vLB2.length()<<endl;
    for (int j = 0; j < oldVars1; j++) vLB1[j]=vLB2[j]; 
    for (int j = 0; j < oldCons1; j++) vLB1[newVars1+j]=vLB2[oldVars1+j];
    denseVector &vUB1 = ubModifiableWithCuts.getFirstStageVec();
    const denseVector &vUB2 = ub.getFirstStageVec();
    for (int j = 0; j < oldVars1; j++) vUB1[j]=vUB2[j]; 
    for (int j = 0; j < oldCons1; j++) vUB1[newVars1+j]=vUB2[oldVars1+j];
    for (int scen = 0; scen < input.nScenarios(); scen++) {
      if(pipsCtx.assignedScenario(scen)) {

        int oldVars2= originalDims.inner.numSecondStageVars(scen);
        int newVars2= currentDims.inner.numSecondStageVars(scen);
        int oldCons2= originalDims.inner.numSecondStageCons(scen);
        int newCons2= currentDims.inner.numSecondStageCons(scen);
   
        denseVector &vLB1 = lbModifiableWithCuts.getSecondStageVec(scen);
        const denseVector &vLB2 = lb.getSecondStageVec(scen);
        for (int j = 0; j < oldVars2; j++) vLB1[j]=vLB2[j];
        for (int j = 0; j < oldCons2; j++) vLB1[newVars2+j]=vLB2[oldVars2+j];

        denseVector &vUB1 = ubModifiableWithCuts.getSecondStageVec(scen);
        const denseVector &vUB2 = ub.getSecondStageVec(scen);
        for (int j = 0; j < oldVars2; j++) vUB1[j]=vUB2[j];
        for (int j = 0; j < oldCons2; j++) vUB1[newVars2+j]=vUB2[oldVars2+j];


      }
    }


    //
    //TODO: there may be some possible leaks here


}

