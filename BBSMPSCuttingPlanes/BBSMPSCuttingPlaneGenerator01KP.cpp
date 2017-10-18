#include "BBSMPSCuttingPlaneGenerator01KP.hpp"

using namespace std;

bool finalKeySort (pair < int, double > i, pair < int, double > j) { 
    if (i.first==j.first) return i.second<j.second;
    return (i.first<j.first); 
}

bool pairObjectSortInc (pair < BAIndex, double > i, pair < BAIndex, double > j) { 
    if (i.second==j.second) {
        if(i.first.idx==j.first.idx)return i.first.scen<j.first.scen;
        return i.first.idx<j.first.idx;
    }
    return (i.second<j.second); 
}
bool pairObjectSortDec (pair < BAIndex, double > i, pair < BAIndex, double > j) { 

    if (i.second==j.second) {
        if(i.first.idx==j.first.idx)return i.first.scen<j.first.scen;
            return i.first.idx<j.first.idx;
    }
    return (i.second>j.second); 

}




bool detectKnapsackFirstRow(int i){

    PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
    const BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getOriginalBADimensionsSlacks();
    BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
    CoinShallowPackedVector row = rootSolver.retrieveARow(i);

    int nRowElems = row.getNumElements();
    for (int j=0; j< nRowElems; j++){
        int col=row.getIndices()[j];
        double colLB = BBSMPSSolver::instance()->getOriginalLB().getFirstStageVec()[col];
        double colUB = BBSMPSSolver::instance()->getOriginalUB().getFirstStageVec()[col];
        if(!isBinary(colLB, colUB, 1e-6) || !input.isFirstStageColInteger(col)|| row.getElements()[j]<0) return false;
    }
    return true;
}

bool detectKnapsackSecondRow(int i, int scen){
    PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
    const BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getOriginalBADimensionsSlacks();
    BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
    CoinShallowPackedVector rowT = rootSolver.retrieveTRow(i,scen); 
    int nRowTElems = rowT.getNumElements();
    for (int j=0; j< nRowTElems; j++){
        int col=rowT.getIndices()[j];
        double colLB = BBSMPSSolver::instance()->getOriginalLB().getFirstStageVec()[col];
        double colUB = BBSMPSSolver::instance()->getOriginalUB().getFirstStageVec()[col];
        if(!isBinary(colLB, colUB, 1e-6)|| !input.isFirstStageColInteger(col) || rowT.getElements()[j]<0) return false;
    }

    CoinShallowPackedVector rowW = rootSolver.retrieveWRow(i,scen);
    int nRowWElems = rowW.getNumElements();
    for (int j=0; j< nRowWElems; j++){
        int col=rowW.getIndices()[j];
        double colLB = BBSMPSSolver::instance()->getOriginalLB().getSecondStageVec(scen)[col];
        double colUB = BBSMPSSolver::instance()->getOriginalUB().getSecondStageVec(scen)[col];
        if(!isBinary(colLB, colUB, 1e-6)|| !input.isSecondStageColInteger(scen,col)|| rowW.getElements()[j]<0) return false;
    }
    return true;
}




BBSMPSCuttingPlaneGenerator01KP::BBSMPSCuttingPlaneGenerator01KP( const char *_name): BBSMPSCuttingPlaneGenerator(_name){
    
    PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
    const BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getOriginalBADimensionsSlacks();
    BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
    totalKnapsackRows=0;
    knapsackRows.allocate(dimsSlacks,ctx, PrimalVector);
    int mype=BBSMPSSolver::instance()->getMype();
    for (int scen = 0; scen <input.nScenarios(); scen++) {
        if (ctx.assignedScenario(scen)){
            int numCons2= dimsSlacks.inner.numSecondStageCons(scen);
            vector<bool> isKnapsack(numCons2,false);
            int nSecondStageKnaps=0;
            for (int i=0; i<numCons2; i++){
                isKnapsack[i]=detectKnapsackSecondRow(i,scen);
                nSecondStageKnaps+=(isKnapsack[i]);
            }
            CoinIndexedVector &v2 = knapsackRows.getSecondStageVec(scen).v;
            v2.clear();
            if (nSecondStageKnaps>0){
                v2.setNumElements(nSecondStageKnaps);
                double *v2Elts = v2.denseVector();
                int *v2Idx = v2.getIndices();
                int ptr=0;
                for(int i=0; i< numCons2; i++){
                    if(isKnapsack[i]){
                        v2Idx[ptr]=i;
                        v2Elts[v2Idx[ptr]]=1;
                        ptr++;
                    }
                }
            }
            totalKnapsackRows+=nSecondStageKnaps;
        }
    }
    MPI_Allreduce(MPI_IN_PLACE,&totalKnapsackRows,1,MPI_INT,MPI_SUM,ctx.comm());

    int numCons1= dimsSlacks.inner.numFirstStageCons();
    vector<bool> isKnapsack(numCons1,false);
    int nFirstStageKnaps=0;
    for (int i=0; i<numCons1; i++){

        isKnapsack[i]=detectKnapsackFirstRow(i);
        nFirstStageKnaps+=(isKnapsack[i]);

    }

    CoinIndexedVector &v1 = knapsackRows.getFirstStageVec().v;
    v1.clear();
    if (nFirstStageKnaps>0){
        v1.setNumElements(nFirstStageKnaps);
        double *v1Elts = v1.denseVector();
        int *v1Idx = v1.getIndices();
        int ptr=0;
        for(int i=0; i< numCons1; i++){
            if(isKnapsack[i]){
                v1Idx[ptr]=i;
                v1Elts[v1Idx[ptr]]=1;
                ptr++;
            }
        }

    }
    totalKnapsackRows+=nFirstStageKnaps;


    

}
bool generateFirstStage01KP(CoinShallowPackedVector row, denseBAVector &LPRelaxationSolution, double lbRHS, double ubRHS, CoinIndexedVector &coverExpression, double &coverLB, double &coverUB){
    int nRowElems = row.getNumElements();
    vector< pair < BAIndex, double > > ratios(nRowElems);

    double a0;
    if (ubRHS==COIN_DBL_MAX){
        double rowElemsAcum=0;
        for (int i=0; i< nRowElems; i++)rowElemsAcum+=row.getElements()[i];
        a0=rowElemsAcum-lbRHS;
    }
    else{
        a0=ubRHS;
    }
    for (int i = 0; i < nRowElems; i++){
        ratios[i].first.scen=-1;
        ratios[i].first.idx=i;

        ratios[i].second=LPRelaxationSolution.getFirstStageVec()[row.getIndices()[i]] ;
        if (ubRHS==COIN_DBL_MAX) ratios[i].second= 1-ratios[i].second;
    }
    
    sort(ratios.begin(),ratios.end(),pairObjectSortDec);
    double w =0;
    int j;
    int coverSize=0;

    for ( j=0; j < nRowElems;j++){
        int idx=ratios[j].first.idx;
        w+=row.getElements()[idx];
        coverSize++;
        if (w >a0){
            break;
        }
    }
    
    vector<bool> isCovered(coverSize,true);
    vector< pair < BAIndex, double > > secondOrdering(coverSize);
    for (int i = 0; i < coverSize; i++){
        secondOrdering[i].first.scen=-1;
        secondOrdering[i].first.idx=ratios[i].first.idx;
        secondOrdering[i].second=row.getElements()[ratios[i].first.idx];
    }   
    sort(secondOrdering.begin(),secondOrdering.end(),pairObjectSortInc);
    for (int i = 0; i < secondOrdering.size(); i++){
        int idx=secondOrdering[i].first.idx;
        if (w-row.getElements()[idx]>a0){
            w-=row.getElements()[idx];
            coverSize--;
            isCovered[i]=false;
        }
        
    }   


    //MAKE C1 and C2, F and R
    vector<BAIndex> C1;
    vector<BAIndex> C2;
    vector<BAIndex> F;
    vector<BAIndex> R;
    for (int i = 0; i < secondOrdering.size(); i++){
        BAIndex aux;
        aux.scen=-1;
        aux.idx=secondOrdering[i].first.idx;
        if (isCovered[i]){
            
            double solValue=LPRelaxationSolution.getFirstStageVec()[row.getIndices()[secondOrdering[i].first.idx]];
            if (ubRHS==COIN_DBL_MAX) solValue=1-solValue;

            if (solValue==1)C2.push_back(aux);
            else C1.push_back(aux);
        }
        else{
            double solValue=LPRelaxationSolution.getFirstStageVec()[row.getIndices()[secondOrdering[i].first.idx]];
            if (ubRHS==COIN_DBL_MAX) solValue=1-solValue;
            if (solValue==0) R.push_back(aux);
            else F.push_back(aux);
        }
    }
    for (int i = secondOrdering.size(); i < nRowElems; i++){
        BAIndex aux;
        aux.scen=-1;
        aux.idx=ratios[i].first.idx;
        double solValue=LPRelaxationSolution.getFirstStageVec()[aux.idx];
        if (ubRHS==COIN_DBL_MAX) solValue=1-solValue;
        if (solValue==0) R.push_back(aux);
        else F.push_back(aux);
    }

    while(C1.size()<2 && C2.size()>0){
        C1.push_back(C2[0]);
        C2.erase(C2.begin());
    }

    
    vector<double> A(nRowElems*100,0);
    vector<double> alphas(nRowElems,0);
    vector<double> z(nRowElems,0);
    vector<double> elementsOrderedByGroup(nRowElems);
    vector<int> idxsOrderedByGroup(nRowElems);
    
    for (int i=0; i<F.size();i++){
        elementsOrderedByGroup[i]=row.getElements()[F[i].idx];
        idxsOrderedByGroup[i]=row.getIndices()[F[i].idx];
    }
    for (int i=0; i<C2.size();i++){
        elementsOrderedByGroup[i+F.size()]=row.getElements()[C2[i].idx];
        idxsOrderedByGroup[i+F.size()]=row.getIndices()[C2[i].idx];
    }
    for (int i=0; i<R.size();i++){
        elementsOrderedByGroup[i+F.size()+C2.size()]=row.getElements()[R[i].idx];
        idxsOrderedByGroup[i+F.size()+C2.size()]=row.getIndices()[R[i].idx];
    }
    for (int om =1; om<= C1.size(); om++) {

        A[om]=A[om-1]+row.getElements()[C1[om-1].idx];
    }
    int alpha0 = C1.size()-1;
    
    double cumulativeAlphasSum=0;
    double C2ElemsSum=0;
    for (int j=F.size();j<C2.size()+F.size();j++){
        C2ElemsSum+=elementsOrderedByGroup[j];
    }
    //UPLIFTING F VARIABLES 
    for (int i=0; i<F.size();i++){
        if(a0-C2ElemsSum-elementsOrderedByGroup[i]<0)z[i]=0;
        else {
            for (int om=alpha0; om>=0; om--){
                if (A[om]<=a0-C2ElemsSum-elementsOrderedByGroup[i]){
                    z[i]=om;
                    break;
                }
            }
        }
        alphas[i]=alpha0-z[i];
        if (C1.size()+cumulativeAlphasSum+alphas[i]+2> A.size()){
            A.resize((C1.size()+cumulativeAlphasSum+alphas[i]+2)*2,0);
        }
        for (int om=C1.size()+cumulativeAlphasSum+1; om <= C1.size()+cumulativeAlphasSum+alphas[i];om++){
            A[om]=9999999999;
        }
        for (int om=C1.size()+cumulativeAlphasSum+alphas[i]; om>=0; om--){
            if(om<alphas[i]){
                if(A[om]>elementsOrderedByGroup[i]){
                    A[om]=elementsOrderedByGroup[i];
                }
            }
            else{
                if(A[om]>A[om-alphas[i]]+elementsOrderedByGroup[i]){
                    A[om]=A[om-alphas[i]]+elementsOrderedByGroup[i];

                }
            }
        }
        cumulativeAlphasSum+=alphas[i];
    }

    double remainingC2Sum=C2ElemsSum;
    //DOWNLIFTING VARIABLES IN C2
    for (int i=F.size(); i<F.size()+C2.size(); i++){
        for (int om= C1.size()+cumulativeAlphasSum; om>=0; om--){
            if (A[om]<=(a0-remainingC2Sum+elementsOrderedByGroup[i])){
                    z[i]=om;
                    break;
                }
        }
        
        remainingC2Sum-=elementsOrderedByGroup[i];
        alphas[i]=z[i]-alpha0;
        alpha0=alpha0+alphas[i];
        if (C1.size()+cumulativeAlphasSum+alphas[i]+2> A.size()){
            A.resize((C1.size()+cumulativeAlphasSum+alphas[i]+2)*2,0);
        }
        for (int om=C1.size()+cumulativeAlphasSum+1; om <= C1.size()+cumulativeAlphasSum+alphas[i];om++){
            A[om]=9999999999;
        }
        for (int om=C1.size()+cumulativeAlphasSum+alphas[i]; om>=0; om--){
            if(om<alphas[i]){
                if(A[om]>elementsOrderedByGroup[i]){
                    A[om]=elementsOrderedByGroup[i];
                }
            }
            else{
                if(A[om]>A[om-alphas[i]]+elementsOrderedByGroup[i]){
                    A[om]=A[om-alphas[i]]+elementsOrderedByGroup[i];
                }
                
            }
        }
        cumulativeAlphasSum+=alphas[i];
    }
    //UPLIFTING VARIABLES IN R
    for (int i=F.size()+C2.size(); i< F.size()+C2.size()+R.size(); i++){
        for (int om=alpha0; om>=0; om--){
            if (A[om]<=a0-elementsOrderedByGroup[i]){
                z[i]=om;
                break;
            }
        }
      
        alphas[i]=alpha0-z[i];
        for (int om=alpha0; om>=0; om--){
            if(om<alphas[i]){
                if(A[om]>elementsOrderedByGroup[i])A[om]=elementsOrderedByGroup[i];
            }
            else{
                if(A[om]>A[om-alphas[i]]+elementsOrderedByGroup[i]){
                    A[om]=A[om-alphas[i]]+elementsOrderedByGroup[i];
                }
            }
        }
    }


    vector< pair< int, double > > finalAlphas(nRowElems);

    int ptr=0;

    for (int i=0; i<C1.size();i++){
        

        finalAlphas[ptr].second=1;
        finalAlphas[ptr].first=row.getIndices()[C1[i].idx];
        ptr++;
        
        
    }

    for (int i=0; i<F.size();i++){
        if(alphas[i]>0){
            finalAlphas[ptr].second=alphas[i];
            finalAlphas[ptr].first=row.getIndices()[F[i].idx];
            ptr++;
        }
    }
    
    for (int i=0; i<C2.size();i++){
        if(alphas[i+F.size()]>0){
            finalAlphas[ptr].second=alphas[i+F.size()];
            finalAlphas[ptr].first=row.getIndices()[C2[i].idx];
            ptr++;
        }
        
    }
    for (int i=0; i<R.size();i++){
        if(alphas[i+F.size()+C2.size()]>0){
            finalAlphas[ptr].second=alphas[i+F.size()+C2.size()];
            finalAlphas[ptr].first=row.getIndices()[R[i].idx];
            ptr++;
        }
    }
    finalAlphas.resize(ptr);
    sort(finalAlphas.begin(),finalAlphas.end(),finalKeySort);
    coverExpression.setNumElements(ptr);
    double *v1Elts = coverExpression.denseVector();
    int *v1Idx = coverExpression.getIndices();

    double rowElemsAcum=0;
    double exprAcum=0;
    for(int i=0; i< ptr; i++){
        
        v1Idx[i]=finalAlphas[i].first;
        v1Elts[v1Idx[i]]=finalAlphas[i].second;
        rowElemsAcum+=finalAlphas[i].second;
        exprAcum+=LPRelaxationSolution.getFirstStageVec()[v1Idx[i]]*v1Elts[v1Idx[i]];

    }

    if (ubRHS==COIN_DBL_MAX){
        
        //for (int i=0; i< nRowElems; i++)rowElemsAcum+=row.getElements()[i];
        alpha0=rowElemsAcum-alpha0;
        coverLB=alpha0;
        coverUB=COIN_DBL_MAX;
    }
    else{
        coverLB=COIN_DBL_MIN;
        coverUB=alpha0;
        
    }
    if (exprAcum-coverLB>=-10E-6 && exprAcum-coverUB<=10E-6)return false;

        return true;
}

bool generateSecondStage01KP(CoinShallowPackedVector rowT,CoinShallowPackedVector rowW, denseBAVector &LPRelaxationSolution, double lbRHS, double ubRHS, int scen, CoinIndexedVector &coverExpressionT, CoinIndexedVector &coverExpressionW, double &coverLB, double &coverUB){
    int nRowTElems = rowT.getNumElements();
    int nRowWElems = rowW.getNumElements();
    
    vector< pair < BAIndex, double > > ratios(nRowWElems+nRowTElems);

    double a0;
    if (ubRHS==COIN_DBL_MAX){
        double rowElemsAcum=0;
        for (int i=0; i< nRowWElems; i++)rowElemsAcum+=rowW.getElements()[i];
        for (int i=0; i< nRowTElems; i++)rowElemsAcum+=rowT.getElements()[i];
        a0=rowElemsAcum-lbRHS;
    }
    else{
        a0=ubRHS;
    }

    for (int i = 0; i < nRowTElems; i++){
        
        ratios[i].first.scen=-1;
        ratios[i].first.idx=i;
        ratios[i].second=LPRelaxationSolution.getFirstStageVec()[rowT.getIndices()[i]] ;
        if (ubRHS==COIN_DBL_MAX)ratios[i].second=1-ratios[i].second;
    }
    for (int i = 0; i < nRowWElems; i++){
        
        ratios[i+nRowTElems].first.scen=scen;
        ratios[i+nRowTElems].first.idx=i;
        ratios[i+nRowTElems].second=LPRelaxationSolution.getSecondStageVec(scen)[rowW.getIndices()[i]] ;
        if (ubRHS==COIN_DBL_MAX)ratios[i+nRowTElems].second=1-ratios[i+nRowTElems].second;
    }
    
    sort(ratios.begin(),ratios.end(),pairObjectSortDec);
    double w =0;
    int j;
    int coverSize=0;

    for ( j=0; j < nRowWElems+nRowTElems;j++){
        int idx=ratios[j].first.idx; 
        if (ratios[j].first.scen!=-1){
            w+=rowW.getElements()[idx];
        }
        else{
            w+=rowT.getElements()[idx];
        }
        coverSize++;
        if (w >a0){
            break;
        }
    }
    
    vector<bool> isCovered(coverSize,true);
    vector< pair < BAIndex, double > > secondOrdering(coverSize);
    for (int i = 0; i < coverSize; i++){
        secondOrdering[i].first.scen=ratios[i].first.scen;
        secondOrdering[i].first.idx=ratios[i].first.idx;
        if(ratios[i].first.scen!=-1)secondOrdering[i].second=rowW.getElements()[ratios[i].first.idx];
        else secondOrdering[i].second=rowT.getElements()[ratios[i].first.idx];
    }   
    sort(secondOrdering.begin(),secondOrdering.end(),pairObjectSortInc);
    for (int i = 0; i < secondOrdering.size(); i++){
        int idx=secondOrdering[i].first.idx;
        if(secondOrdering[i].first.scen!=-1){
            if (w-rowW.getElements()[idx]>a0){
                w-=rowW.getElements()[idx];
                coverSize--;
                isCovered[i]=false;
            }
        }
        else {
            if (w-rowT.getElements()[idx]>a0){
                w-=rowT.getElements()[idx];
                coverSize--;
                isCovered[i]=false;
            }
        }
    }   
    

    //MAKE C1 and C2, F and R
    vector<BAIndex> C1;
    vector<BAIndex> C2;
    vector<BAIndex> F;
    vector<BAIndex> R;
    for (int i = 0; i < secondOrdering.size(); i++){
        BAIndex aux;
        aux.scen=secondOrdering[i].first.scen;
        aux.idx=secondOrdering[i].first.idx;
        double solValue;
        if (aux.scen==-1)solValue=LPRelaxationSolution.getFirstStageVec()[rowT.getIndices()[secondOrdering[i].first.idx]];
        else solValue=LPRelaxationSolution.getSecondStageVec(scen)[rowW.getIndices()[secondOrdering[i].first.idx]];
        if (ubRHS==COIN_DBL_MAX)solValue=1-solValue;
        if (isCovered[i]){
            if (solValue==1)C2.push_back(aux);
            else C1.push_back(aux);
        }
        else{
            if (solValue==0) R.push_back(aux);
            else F.push_back(aux);
        }
    }
    for (int i = secondOrdering.size(); i < nRowWElems+nRowTElems; i++){
        BAIndex aux;
        aux.scen=ratios[i].first.scen;
        aux.idx=ratios[i].first.idx;
        double solValue;
        if (aux.scen==-1)solValue=LPRelaxationSolution.getFirstStageVec()[rowT.getIndices()[ratios[i].first.idx]];
        else solValue=LPRelaxationSolution.getSecondStageVec(scen)[rowW.getIndices()[ratios[i].first.idx]];
        if (ubRHS==COIN_DBL_MAX)solValue=1-solValue;
        if (solValue==0) R.push_back(aux);
        else F.push_back(aux);
    }

    while(C1.size()<2 && C2.size()>0){
        C1.push_back(C2[0]);
        C2.erase(C2.begin());
    }

    
    vector<double> A(nRowWElems+nRowTElems*100 ,0);
    vector<double> alphas(nRowWElems+nRowTElems,0);
    vector<double> z(nRowWElems+nRowTElems,0);

    vector<double> elementsOrderedByGroup(nRowWElems+nRowTElems);
    vector<int> idxsOrderedByGroup(nRowWElems+nRowTElems);
    vector<int> scensOrderedByGroup(nRowWElems+nRowTElems);
    for (int i=0; i<F.size();i++){
        if (F[i].scen==-1){
            elementsOrderedByGroup[i]=rowT.getElements()[F[i].idx];
            idxsOrderedByGroup[i]=rowT.getIndices()[F[i].idx];
            scensOrderedByGroup[i]=-1;
        }
        else{
            elementsOrderedByGroup[i]=rowW.getElements()[F[i].idx];
            idxsOrderedByGroup[i]=rowW.getIndices()[F[i].idx];
            scensOrderedByGroup[i]=scen;
        }
        
    }

    for (int i=0; i<C2.size();i++){
        if (C2[i].scen==-1){
            elementsOrderedByGroup[i+F.size()]=rowT.getElements()[C2[i].idx];
            idxsOrderedByGroup[i+F.size()]=rowT.getIndices()[C2[i].idx];
            scensOrderedByGroup[i+F.size()]=-1;
        }
        else{
            elementsOrderedByGroup[i+F.size()]=rowW.getElements()[C2[i].idx];
            idxsOrderedByGroup[i+F.size()]=rowW.getIndices()[C2[i].idx];
            scensOrderedByGroup[i+F.size()]=scen;
        }
    }

    for (int i=0; i<R.size();i++){
        if (R[i].scen==-1){
            elementsOrderedByGroup[i+F.size()+C2.size()]=rowT.getElements()[R[i].idx];
            idxsOrderedByGroup[i+F.size()+C2.size()]=rowT.getIndices()[R[i].idx];
            scensOrderedByGroup[i+F.size()+C2.size()]=-1;
        }
        else{
            elementsOrderedByGroup[i+F.size()+C2.size()]=rowW.getElements()[R[i].idx];
            idxsOrderedByGroup[i+F.size()+C2.size()]=rowW.getIndices()[R[i].idx];
            scensOrderedByGroup[i+F.size()+C2.size()]=scen;
        }

        
    }
   
    for (int om =1; om<= C1.size(); om++){
        if (C1[om-1].scen==-1){
            A[om]=A[om-1]+rowT.getElements()[C1[om-1].idx];
        }
        else{
            A[om]=A[om-1]+rowW.getElements()[C1[om-1].idx];
        }

    } 
    int alpha0 = C1.size()-1;
    
    double cumulativeAlphasSum=0;
    double C2ElemsSum=0;
    for (int j=F.size();j<C2.size()+F.size();j++){
        C2ElemsSum+=elementsOrderedByGroup[j];
    }
    //UPLIFTING F VARIABLES
    for (int i=0; i<F.size();i++){
        if(a0-C2ElemsSum-elementsOrderedByGroup[i]<0)z[i]=0;
        else {
            for (int om=alpha0; om>=0; om--){
                if (A[om]<=a0-C2ElemsSum-elementsOrderedByGroup[i]){
                    z[i]=om;
                    break;
                }
            }
        }
        
        alphas[i]=alpha0-z[i];
        if (C1.size()+cumulativeAlphasSum+alphas[i]+2> A.size()){
            A.resize((C1.size()+cumulativeAlphasSum+alphas[i]+2)*2,0);
        }
        
        for (int om=C1.size()+cumulativeAlphasSum+1; om <= C1.size()+cumulativeAlphasSum+alphas[i];om++){
            A[om]=9999999999;
        
        }
        for (int om=C1.size()+cumulativeAlphasSum+alphas[i]; om>=0; om--){
            if(om<alphas[i]){
                if(A[om]>elementsOrderedByGroup[i]){
                    A[om]=elementsOrderedByGroup[i];
                }
            }
            else{
                if(A[om]>A[om-alphas[i]]+elementsOrderedByGroup[i]){
                    
                    A[om]=A[om-alphas[i]]+elementsOrderedByGroup[i];
                }
            }
        }
        cumulativeAlphasSum+=alphas[i];
    }

    double remainingC2Sum=C2ElemsSum;
    //DOWNLIFTING VARIABLES IN C2
    for (int i=F.size(); i<F.size()+C2.size(); i++){
        for (int om= C1.size()+cumulativeAlphasSum; om>=0; om--){
            if (A[om]<=(a0-remainingC2Sum+elementsOrderedByGroup[i])){
                    z[i]=om;
                    break;
                }
        }
        
        remainingC2Sum-=elementsOrderedByGroup[i];
        alphas[i]=z[i]-alpha0;
        alpha0=alpha0+alphas[i];
        if (C1.size()+cumulativeAlphasSum+alphas[i]+2> A.size()){
            A.resize((C1.size()+cumulativeAlphasSum+alphas[i]+2)*2,0);
        }

        for (int om=C1.size()+cumulativeAlphasSum+1; om <= C1.size()+cumulativeAlphasSum+alphas[i];om++){
            A[om]=9999999999;
        
        }
        for (int om=C1.size()+cumulativeAlphasSum+alphas[i]; om>=0; om--){
            if(om<alphas[i]){
                if(A[om]>elementsOrderedByGroup[i]){
                    A[om]=elementsOrderedByGroup[i];
                }
            }
            else{
                if(A[om]>A[om-alphas[i]]+elementsOrderedByGroup[i]){
                    
                    A[om]=A[om-alphas[i]]+elementsOrderedByGroup[i];
                }
                
            }
        }
        cumulativeAlphasSum+=alphas[i];
    }
    //UPLIFTING VARIABLES IN R
    for (int i=F.size()+C2.size(); i< F.size()+C2.size()+R.size(); i++){
        for (int om=alpha0; om>=0; om--){
            if (A[om]<=a0-elementsOrderedByGroup[i]){
                z[i]=om;
                break;
            }
        }
        alphas[i]=alpha0-z[i];
        
        for (int om=alpha0; om>=0; om--){
            
            if(om<alphas[i]){
                if(A[om]>elementsOrderedByGroup[i]){
                    A[om]=elementsOrderedByGroup[i];
                }
                
            }
            else{
               
                    
                if(A[om]>A[om-alphas[i]]+elementsOrderedByGroup[i]){
                    A[om]=A[om-alphas[i]]+elementsOrderedByGroup[i];
                 }
                
            }
        }
    }
    vector< pair< int, double > > finalAlphasW(nRowWElems);
    vector< pair< int, double > > finalAlphasT(nRowTElems);

    int tPtr=0;
    int wPtr=0;
    

    for (int i=0; i<C1.size();i++){
        
        if(C1[i].scen==-1){
            finalAlphasT[tPtr].second=1;
            finalAlphasT[tPtr].first=rowT.getIndices()[C1[i].idx];
            tPtr++;
        }
        else {
            finalAlphasW[wPtr].second=1;
            finalAlphasW[wPtr].first=rowW.getIndices()[C1[i].idx];
            wPtr++;
        }

        
    }


    for (int i=0; i<F.size();i++){
        if (alphas[i]>0){
            if(scensOrderedByGroup[i]==-1){
                finalAlphasT[tPtr].second=alphas[i];
                finalAlphasT[tPtr].first=rowT.getIndices()[F[i].idx];
                tPtr++;
            }
            else {
                finalAlphasW[wPtr].second=alphas[i];
                finalAlphasW[wPtr].first=rowW.getIndices()[F[i].idx];
                wPtr++;
            }

        }
        
    }
    for (int i=0; i<C2.size();i++){
        if (alphas[i+F.size()]>0){
            if(scensOrderedByGroup[i+F.size()]==-1){
                finalAlphasT[tPtr].second=alphas[i+F.size()];
                finalAlphasT[tPtr].first=rowT.getIndices()[C2[i].idx];
                tPtr++;
            }
            else {
                finalAlphasW[wPtr].second=alphas[i+F.size()];
                finalAlphasW[wPtr].first=rowW.getIndices()[C2[i].idx];
                wPtr++;
            } 
        }
    }
    for (int i=0; i<R.size();i++){
        if (alphas[i+F.size()+C2.size()]>0){
            if(scensOrderedByGroup[i+F.size()+C2.size()]==-1){
                finalAlphasT[tPtr].second=alphas[i+F.size()+C2.size()];
                finalAlphasT[tPtr].first=rowT.getIndices()[R[i].idx];
                tPtr++;
            }
            else {
                finalAlphasW[wPtr].second=alphas[i+F.size()+C2.size()];
                finalAlphasW[wPtr].first=rowW.getIndices()[R[i].idx];
                wPtr++;
            }
        }
    }
    finalAlphasT.resize(tPtr);
    finalAlphasW.resize(wPtr);
    sort(finalAlphasT.begin(),finalAlphasT.end(),finalKeySort);
    sort(finalAlphasW.begin(),finalAlphasW.end(),finalKeySort);
    coverExpressionT.setNumElements(tPtr);
    double *v1Elts = coverExpressionT.denseVector();
    int *v1Idx = coverExpressionT.getIndices();
    int ptr=0;
    double rowElemsAcum=0;
    double exprAcum=0;
    for(int i=0; i< tPtr; i++){
            rowElemsAcum+=finalAlphasT[i].second;
            v1Idx[ptr]=finalAlphasT[i].first;
            v1Elts[v1Idx[ptr]]=finalAlphasT[i].second;
            exprAcum+=LPRelaxationSolution.getFirstStageVec()[v1Idx[ptr]]*v1Elts[v1Idx[ptr]];
            ptr++;

        
    }

    coverExpressionW.setNumElements(wPtr);
    v1Elts = coverExpressionW.denseVector();
    v1Idx = coverExpressionW.getIndices();
    ptr=0;
    for(int i=0; i< wPtr; i++){
        rowElemsAcum+=finalAlphasW[i].second;
        v1Idx[ptr]=finalAlphasW[i].first;
        v1Elts[v1Idx[ptr]]=finalAlphasW[i].second;
        exprAcum+=LPRelaxationSolution.getSecondStageVec(scen)[v1Idx[ptr]]*v1Elts[v1Idx[ptr]];
        ptr++;
        
    }

    

    if (ubRHS==COIN_DBL_MAX){
        

        alpha0=rowElemsAcum-alpha0;
        coverLB=alpha0;
        coverUB=COIN_DBL_MAX;
    }
    else{
        coverLB=COIN_DBL_MIN;
        coverUB=alpha0;
        
    }
    if (exprAcum-coverLB>=-10E-6 && exprAcum-coverUB<=10E-6)return false;

    return true;
    
}
bool BBSMPSCuttingPlaneGenerator01KP::generateCuttingPlane(BBSMPSNode* node, denseBAVector &LPRelaxationSolution){

    bool planesAdded=false;
    PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
    const BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getOriginalBADimensionsSlacks();
    BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
    int mype=BBSMPSSolver::instance()->getMype();
    
    CoinIndexedVector &kr1 = knapsackRows.getFirstStageVec().v;

    for (int i=0; i<kr1.getNumElements(); i++){

        int rowIdx = kr1.getIndices()[i];
        //For each row in first stage determine if it is a knapsack.
        CoinShallowPackedVector row = rootSolver.retrieveARow(rowIdx);

        double lbRHS=BBSMPSSolver::instance()->getOriginalLB().getFirstStageVec()[dimsSlacks.inner.numFirstStageVars()+rowIdx];
        double ubRHS=BBSMPSSolver::instance()->getOriginalUB().getFirstStageVec()[dimsSlacks.inner.numFirstStageVars()+rowIdx];

        sparseBAVector cutExpression;
        cutExpression.allocate(dimsSlacks,ctx, PrimalVector);
        CoinIndexedVector &v1 = cutExpression.getFirstStageVec().v;
        v1.clear();
        double coverLB;
        double coverUB;
        int valid=0;
        if (mype==0){
            valid+=generateFirstStage01KP(row,LPRelaxationSolution,lbRHS,ubRHS,v1,coverLB,coverUB);
        }
        MPI_Bcast(&valid,1,MPI_INT,0,ctx.comm());
    
        if(valid){

            int nElems=0;
            if (mype==0)nElems=v1.getNumElements();
            MPI_Bcast(&nElems,1,MPI_INT,0,ctx.comm());
            MPI_Bcast(&coverLB,1,MPI_DOUBLE,0,ctx.comm());
            MPI_Bcast(&coverUB,1,MPI_DOUBLE,0,ctx.comm());
            if(mype!=0)v1.setNumElements(nElems);
            MPI_Bcast(v1.getIndices(),nElems,MPI_INT,0,ctx.comm());
            int lastElem=v1.getIndices()[nElems-1];
            MPI_Bcast(v1.denseVector(),lastElem+1,MPI_DOUBLE,0,ctx.comm());


            CoinIndexedVector &v1222 = cutExpression.getFirstStageVec().v;
            double indexSum=0;
            double cuofSum=0;
            for (int j=0; j< v1222.getNumElements();j++){
                indexSum+=v1222.getIndices()[j];
                cuofSum+=v1222.getIndices()[j]*v1222.denseVector()[v1222.getIndices()[j]];
            }
            BBSMPSCuttingPlane *plane= new BBSMPSCuttingPlane(coverLB,coverUB,cutExpression);
            node->addCuttingPlane(plane);
            planesAdded=true;
            cuttingPlanesGenerated++;
        }

    }
    for (int scen = 0; scen <input.nScenarios(); scen++) {



        
        int numCons2= -1;
        CoinIndexedVector kr2;
        if (ctx.assignedScenario(scen)){
            kr2 = knapsackRows.getSecondStageVec(scen).v;

            numCons2= kr2.getNumElements();
        }
        MPI_Allreduce(MPI_IN_PLACE,&numCons2,1,MPI_INT,MPI_MAX,ctx.comm());
        for (int i=0; i<numCons2; i++){
            
            //Stage 2: generate initial cover
           
            sparseBAVector cutExpression;
            cutExpression.allocate(dimsSlacks,ctx, PrimalVector);

            double coverLB=COIN_DBL_MIN;
            double coverUB=COIN_DBL_MAX;
            int valid=0;
            if(ctx.assignedScenario(scen)) {
                int rowIdx = kr2.getIndices()[i];
                CoinShallowPackedVector rowT = rootSolver.retrieveTRow(rowIdx,scen);
                CoinShallowPackedVector rowW = rootSolver.retrieveWRow(rowIdx,scen);
                double lbRHS=BBSMPSSolver::instance()->getOriginalLB().getSecondStageVec(scen)[dimsSlacks.inner.numSecondStageVars(scen)+rowIdx];
                double ubRHS=BBSMPSSolver::instance()->getOriginalUB().getSecondStageVec(scen)[dimsSlacks.inner.numSecondStageVars(scen)+rowIdx];
                CoinIndexedVector &v1 = cutExpression.getFirstStageVec().v;
                CoinIndexedVector &v2 = cutExpression.getSecondStageVec(scen).v;
                v1.clear();
                v2.clear();
                valid=generateSecondStage01KP(rowT,rowW,LPRelaxationSolution,lbRHS,ubRHS,scen,v1,v2,coverLB,coverUB);
            }
            MPI_Allreduce(MPI_IN_PLACE,&valid,1,MPI_INT,MPI_SUM,ctx.comm());
            if(valid){

                BBSMPSCuttingPlane *plane= new BBSMPSCuttingPlane(coverLB,coverUB,cutExpression);
                node->addCuttingPlane(plane);
                planesAdded=true;
                cuttingPlanesGenerated++;
            }
                
        }
    }
    

    return planesAdded;


}

bool BBSMPSCuttingPlaneGenerator01KP::shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution){
    
    return (totalKnapsackRows>0 && node->getNodeDepth()<=0);

}