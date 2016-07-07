#include "BBSMPSNode.hpp"

using namespace std;

int BBSMPSNode::nodeCounter=0;


// Construct node from {lower, upper} bounds, obj val of parent node
BBSMPSNode::BBSMPSNode(double _objectiveValue, const std::vector< std::pair < BAIndex, variableState > > &states, int procOfOrigin):
partialStartState(states),
procNumber(procOfOrigin) {
  parent =NULL;
  objectiveValue=_objectiveValue;
  childrenAlive=0;
  nodeNumber=(++nodeCounter);
  nodeDepth=-1;
  if (nodeCounter==1) nodeDepth=0;

}


  // Add copy constructor so that priority_queue can use it for
  // instantiation because denseBAVector and BAFlagVector do not have
  // assignment operators or copy constructors.
BBSMPSNode::BBSMPSNode(const BBSMPSNode &sourceNode):
partialStartState(sourceNode.partialStartState),
parent(sourceNode.parent),
objectiveValue(sourceNode.objectiveValue),
childrenAlive(sourceNode.childrenAlive),
nodeNumber(sourceNode.nodeNumber),
nodeDepth(sourceNode.nodeDepth),
procNumber(sourceNode.procNumber) {
}
BBSMPSNode::BBSMPSNode(BBSMPSNode* parent_ptr, std::vector<BBSMPSBranchingInfo>& bInfos, int procOfOrigin){

  
  if (parent_ptr!=NULL){
    parent=parent_ptr;
    parent->incrementAliveChildren();
    nodeDepth=parent->getNodeDepth()+1;
  }
  else {
    parent=NULL;
    nodeDepth=0;

  }

  objectiveValue=-INFINITY;
  branchingInfos=std::vector<BBSMPSBranchingInfo>(bInfos);
  childrenAlive=0;
  nodeNumber=(++nodeCounter);
  procNumber=procOfOrigin;
}

BBSMPSNode::BBSMPSNode(int *intVector, double *dblVector){

  objectiveValue=dblVector[0];
  nodeNumber=intVector[0];
  procNumber=intVector[1];
  nodeDepth=intVector[2];
  parent=NULL;
  childrenAlive=0;

  partialStartState.resize(intVector[3]);
  vector<variableState> auxVector(4,Basic);
  auxVector[1]=Basic;
  auxVector[2]=AtLower;
  auxVector[3]=AtUpper;
  int ptr=5;
  for(int i=0; i<partialStartState.size();i++){
    partialStartState[i].first.scen=intVector[ptr];
    partialStartState[i].first.idx=intVector[ptr+1];
    partialStartState[i].second=auxVector[intVector[ptr+2]];
    ptr+=3;
 
  }
   // //Incremental branching information relative to this node
  branchingInfos.resize(intVector[4]);

  int intSerialSize=0;
  int dblSerialSize=0;
  BBSMPSBranchingInfo::getSerializationSize(intSerialSize, dblSerialSize);
  int intPtr=ptr;
  int dblPtr=1;
  for (int i=0; i< branchingInfos.size(); i++){
    branchingInfos[i]=BBSMPSBranchingInfo(&intVector[intPtr],&dblVector[dblPtr]);
    intPtr+=intSerialSize;
    dblPtr+=dblSerialSize;
  }
}

BBSMPSNode::~BBSMPSNode(){
  assert(childrenAlive==0);
  if(parent!=NULL){
    parent->decrementAliveChildren();
  }
  for (int i=0; i<cuttingPlanes.size();i++) {
    if (cuttingPlanes[i])delete cuttingPlanes[i];
  }
}

double BBSMPSNode::getParentObjective() const{
	if (parent!=NULL){
    return parent->getObjective();
  }
  return -INFINITY;
}

double BBSMPSNode::getObjective() const{
  return objectiveValue;
}

void BBSMPSNode::setParentPtr(BBSMPSNode* parent_ptr){
  parent=&(*parent_ptr);
}

BBSMPSNode* BBSMPSNode::getParentPtr(){
  return parent;
}



void BBSMPSNode::setObjective(double lb){
  objectiveValue=lb;
}


void BBSMPSNode::addBranchingInformation(BBSMPSBranchingInfo& bi){
  branchingInfos.push_back(bi);
}


void BBSMPSNode::auxCopyAllBranchingInformation(std::vector<BBSMPSBranchingInfo> &biVector){
  if (parent!=NULL)parent->auxCopyAllBranchingInformation(biVector);
  if (branchingInfos.size()>0)  {
   
    biVector.insert(biVector.end(),branchingInfos.begin(), branchingInfos.end());
  }
  
}

void BBSMPSNode::auxCopyAllCuttingPlaneInformation(std::vector<BBSMPSCuttingPlane*> &cpVector){
  if (parent!=NULL)parent->auxCopyAllCuttingPlaneInformation(cpVector);
   for (int i=0; i< cuttingPlanes.size(); i++){
    cpVector.push_back(cuttingPlanes[i]);
  }  
  
}

void BBSMPSNode::incrementAliveChildren(){
  childrenAlive++;
} 

void BBSMPSNode::decrementAliveChildren(){
  assert(childrenAlive>0);//we should always have children before calling this.

  childrenAlive--;
  if (childrenAlive==0){

    delete this;
  }


}


void BBSMPSNode::eliminate(){

  if (childrenAlive==0){

    delete this;
  }

}

void BBSMPSNode::setIncrementalWarmStartState(std::vector< std::pair < BAIndex, variableState > > &changes){
  partialStartState=changes;
}


void BBSMPSNode::reconstructWarmStartState(BAFlagVector<variableState> &state){

  if (parent!=NULL) parent->reconstructWarmStartState(state);
  for (int i=0; i< partialStartState.size(); i++){
   // if (state.getVec(partialStartState[i].first.scen).length()<=partialStartState[i].first.idx) std::cout<<"Updating "<<partialStartState[i].first.scen<<" "<<partialStartState[i].first.idx<<" with "<<partialStartState[i].second<<" limit "<<state.getVec(partialStartState[i].first.scen).length()<<endl;
    state.getVec(partialStartState[i].first.scen)[partialStartState[i].first.idx]=partialStartState[i].second;
  }
 // cout<<"we got out "<<endl;
}




int BBSMPSNode::getBranchingInfoSize(){
  return branchingInfos.size();
}

int BBSMPSNode::getPartialStateInfoSize(){
  return partialStartState.size();
}

void BBSMPSNode::getAllBranchingInformation(std::vector<BBSMPSBranchingInfo> &biVector){
  BBSMPSNode* n_ptr = parent;
  int branchingInfoSize=getBranchingInfoSize();
  while(n_ptr!=NULL){
    branchingInfoSize+=n_ptr->getBranchingInfoSize();
    n_ptr=n_ptr->parent;
  }
  biVector.reserve(branchingInfoSize);
  auxCopyAllBranchingInformation(biVector);
}

void BBSMPSNode::getAllBranchingInformation(denseBAVector &lb,denseBAVector &ub){

  if (parent!= NULL) parent->getAllBranchingInformation(lb,ub);
  for (int i=0; i<branchingInfos.size(); i++) branchingInfos[i].applyBranchingInfo(lb,ub);

}

int BBSMPSNode::getNodeNumber()const {
  return nodeNumber;
}

void  BBSMPSNode::setNodeDepth(int depth){
  nodeDepth=depth;
}

int BBSMPSNode::getNodeDepth(){
  return nodeDepth;
}

void BBSMPSNode::addCuttingPlane(BBSMPSCuttingPlane *cp){
  cuttingPlanes.push_back(cp);
  cuttingPlaneUids.push_back(cuttingPlanes[cuttingPlanes.size()-1]->getPlaneUid());
}

int BBSMPSNode::getNumberOfCuttingPlanes(){
  return cuttingPlanes.size();
}

void BBSMPSNode::getAllCuttingPlanes(std::vector<BBSMPSCuttingPlane*> &cpVector){
  BBSMPSNode* n_ptr = parent;
  int nCuttingPlanes=getNumberOfCuttingPlanes();
  while(n_ptr!=NULL){
    nCuttingPlanes+=n_ptr->getNumberOfCuttingPlanes();
    n_ptr=n_ptr->parent;

  }
  cpVector.reserve(nCuttingPlanes);
  auxCopyAllCuttingPlaneInformation(cpVector);

}

void BBSMPSNode::getAllCuttingUids(std::vector<int> &uidVector){
  BBSMPSNode* n_ptr = parent;
  int nCuttingPlanes=getNumberOfCuttingPlanes();
  while(n_ptr!=NULL){
    nCuttingPlanes+=n_ptr->getNumberOfCuttingPlanes();
    n_ptr=n_ptr->parent;

  }
  uidVector.reserve(nCuttingPlanes);
  getCurrentNodeCuttingPlaneUids(uidVector);
  n_ptr = parent;
  while(n_ptr!=NULL){
    n_ptr->getCurrentNodeCuttingPlaneUids(uidVector);
    n_ptr=n_ptr->parent;

  }

}

  void BBSMPSNode::getCurrentNodeCuttingPlaneUids(std::vector<int> &uidVector){
    uidVector.insert(uidVector.end(), cuttingPlaneUids.begin(), cuttingPlaneUids.end());
  }

  void BBSMPSNode::getGrandParentCuttingPlanes(std::vector<BBSMPSCuttingPlane*> &cpVector){
    BBSMPSNode* n_ptr = parent;
    if (n_ptr!=NULL) n_ptr=n_ptr->parent;
    int nCuttingPlanes=0;
    while(n_ptr!=NULL){
      nCuttingPlanes+=n_ptr->getNumberOfCuttingPlanes();
      n_ptr=n_ptr->parent;

    }
    if (nCuttingPlanes>0){
      cpVector.reserve(nCuttingPlanes);
      parent->parent->auxCopyAllCuttingPlaneInformation(cpVector);
    }
    

  }

  void BBSMPSNode::getParentNodeCuttingPlanes(std::vector<BBSMPSCuttingPlane*> &cpVector){
    
    if (parent!=NULL) parent->getCurrentNodeCuttingPlanes(cpVector);


  }

   void BBSMPSNode::getCurrentNodeCuttingPlanes(std::vector<BBSMPSCuttingPlane*> &cpVector){
    for (int i=0; i< cuttingPlanes.size(); i++){
      cpVector.push_back(cuttingPlanes[i]);
    }  

  }


void BBSMPSNode::copyCuttingPlanes(BBSMPSNode *node){
      std::vector<BBSMPSCuttingPlane*> allPlaneVector;
      node->getAllCuttingPlanes(allPlaneVector);
      for (int i=0; i< allPlaneVector.size(); i++){
        BBSMPSCuttingPlane *plane= new BBSMPSCuttingPlane(*(allPlaneVector[i]));
        addCuttingPlane(plane); 
      }
  }



  void BBSMPSNode::getSerializationSize(int &intVectorSize, int &dblVectorSize){
    BBSMPSNode* n_ptr = parent;
    int branchingInfoSize=getBranchingInfoSize();
    while(n_ptr!=NULL){
      branchingInfoSize+=n_ptr->getBranchingInfoSize();
      n_ptr=n_ptr->parent;
    }

    int nWSs=getPartialStateInfoSize();
    n_ptr = parent;
    while(n_ptr!=NULL){
      nWSs+=n_ptr->getPartialStateInfoSize();
      n_ptr=n_ptr->parent;
    }
    int intBranchSerialSize=0;
    int dblBranchSerialSize=0;

    BBSMPSBranchingInfo::getSerializationSize(intBranchSerialSize, dblBranchSerialSize);

    intVectorSize=5+branchingInfoSize*intBranchSerialSize+nWSs*3;
    dblVectorSize=1+branchingInfoSize*dblBranchSerialSize;
    cout<<" int and dbl sizes "<<intVectorSize<<" "<<dblVectorSize<<endl;
  }

  void BBSMPSNode::serialize(int *intVector, double *dblVector){
    dblVector[0]=objectiveValue;
    intVector[0]=nodeNumber;
    intVector[1]=procNumber;
    intVector[2]=nodeDepth;

    BBSMPSNode* n_ptr = parent;
    int branchingInfoSize=getBranchingInfoSize();
    while(n_ptr!=NULL){
      branchingInfoSize+=n_ptr->getBranchingInfoSize();
      n_ptr=n_ptr->parent;
    }

    int nWSs=getPartialStateInfoSize();
    n_ptr = parent;
    while(n_ptr!=NULL){
      nWSs+=n_ptr->getPartialStateInfoSize();
      n_ptr=n_ptr->parent;
    }
    intVector[3]=nWSs;
    intVector[4]=branchingInfoSize;
    int ptrAfterWS=serializePartialStateInfo(&intVector[5]);
    serializeBranchingInfo(&intVector[5+ptrAfterWS],&dblVector[1]);
    cout<<"serializing node "<<nodeNumber<<" "<<procNumber<<" "<<nodeDepth<<" "<<nWSs<<" "<<branchingInfoSize<<endl;
  }

  int BBSMPSNode::serializePartialStateInfo(int* intVector){
    int ptr=0;
    if(parent!=NULL)ptr=parent->serializePartialStateInfo(&intVector[ptr]);
    for (int i=0;i<partialStartState.size(); i++){
      intVector[ptr]=partialStartState[i].first.scen;
      intVector[ptr+1]=partialStartState[i].first.idx;
      intVector[ptr+2]=(partialStartState[i].second==Basic)*1+(partialStartState[i].second==AtLower)*2+(partialStartState[i].second==AtUpper)*3;
      ptr+=3;
    }
    return ptr;
  }


    int BBSMPSNode::serializeBranchingInfo(int *intVector, double *dblVector){

      int intSerialSize=0;
      int dblSerialSize=0;
      int nBranchingInfos=0;
      BBSMPSBranchingInfo::getSerializationSize(intSerialSize, dblSerialSize);
      int intPtr=0;
      int dblPtr=0;
      for (int i=0; i< branchingInfos.size(); i++){
        branchingInfos[i].serialize(&intVector[intPtr],&dblVector[dblPtr]);
        intPtr+=intSerialSize;
        dblPtr+=dblSerialSize;
        nBranchingInfos++;
      }

 
    if(parent!=NULL)nBranchingInfos+=parent->serializeBranchingInfo(&intVector[intPtr],&dblVector[dblPtr]);
    return nBranchingInfos;
  }

  // //Class variable used to assign node numbers upon instantiation
  // static int nodeCounter;

  // int nodeNumber;

  // //Depth of the node within the BB tree
  // int nodeDepth;

  // //Pointer to the node's parent
  // BBSMPSNode *parent;

  // //Child reference counter
  // int childrenAlive;

  // //Incremental branching information relative to this node
  // std::vector<BBSMPSBranchingInfo> branchingInfos;


  // // variable states for warm start information; each index is
  // // one of {Basic, AtLower, AtUpper}

  // std::vector< std::pair < BAIndex, variableState > > partialStartState;

  // // objective function of the actual node (not its parent)
  // double objectiveValue; 


  // //Incremental cutting plane information relative to this node
  // std::vector<BBSMPSCuttingPlane*> cuttingPlanes;
  // vector<int> cuttingPlaneUids;
  // // TODO: Local cut objects; Global cuts are stored in the B&B tree.

  // void auxCopyAllBranchingInformation(std::vector<BBSMPSBranchingInfo> &biVector);

  // void auxCopyAllCuttingPlaneInformation(std::vector<BBSMPSCuttingPlane*> &cpVector);

  // int getBranchingInfoSize();

  // int getNumberOfCuttingPlanes();


  // BBSMPSNode(); // Disallow default constructor
