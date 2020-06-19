/*
Implementing bn.h
*/


#include "bn.h"
#include "iomanip"
BN::BN()
{

}

BN::~BN()
{
    //dtor
}

void print_undirected_network(vector<NODE>& nodes, DATA& data, ostream& os){
    os << "node,neighbor" << endl;
    for(int nodeA_idx=0; nodeA_idx!=nodes.size(); ++nodeA_idx){
        vector<int> this_neighbors = nodes[nodeA_idx].Get_neighbors();
        for(int neighbor_idx=0; neighbor_idx!=this_neighbors.size(); ++neighbor_idx){
            int nodeB_idx = this_neighbors[neighbor_idx];
            os << data.Get_node_name(nodeA_idx) << "," << data.Get_node_name(nodeB_idx) << endl;
        }
    }
    return;
}

void BN::Initial(DATA& data, double cutoff)
{
    _allNodes.reserve(data.Get_NumberOfNodes());
    _score = 0.0;
    _dependencies.reserve(data.Get_NumberOfNodes());
    _numberOfNeighbors.reserve(data.Get_NumberOfNodes());

    _growing = true;
    _oneNeb = false;
    _changedAcc.reserve(data.Get_NumberOfNodes());

    ofstream mitfile;
    mitfile.open(_name + "/" +"MIT.csv", ios::app);
    mitfile << "nodeA,nodeB,prior,data,posterier,score,is significant\n";

    ofstream cmitfile;
    cmitfile.open(_name + "/" + "CMIT.csv", ios::app);
    cmitfile << "nodeA,nodeB,givenNode,prior,data,posterier,score,is significant\n";

    for(int i = 0 ; i<data.Get_NumberOfNodes();i++)
    {
        NODE helpNode;
        helpNode.Set_nodeID(i);
        helpNode.Set_numberOfParameter(data.Get_param(i));
        helpNode.Get_correspondingCPD().Set_correspondingNode(i);
        helpNode.Get_correspondingCPD().Set_numberOfParameter(data.Get_param(i));
        helpNode.UpdateCPD(data);
        helpNode.UpdateBIC(data);
        helpNode.Get_correspondingConstraint().Set_correspondingNode(i);

        helpNode.Get_correspondingConstraint().Set_cutoff(cutoff);

        helpNode.Get_correspondingConstraint().MIT(data, mitfile);

        helpNode.Set_neighbors(helpNode.Get_correspondingConstraint().Get_dependentNode());
        _dependencies.push_back(helpNode.Get_correspondingConstraint().Get_dependentNode().size());
        _allNodes.push_back(helpNode);

        _score += helpNode.Get_scoreContribution();
        _changedAcc.push_back(true);
    }
    mitfile.close();
    ofstream afterMIT;
    afterMIT.open(_name + "/" + "afterMIT.csv");
    print_undirected_network(_allNodes, data, afterMIT);
    afterMIT.close();
    for(int i = 0 ; i<data.Get_NumberOfNodes();i++)
    {
        _allNodes[i].CMIT(_allNodes, data, cutoff, cmitfile);
    }
    cmitfile.close();
    ofstream afterCMIT;
    afterCMIT.open(_name + "/" + "afterCMIT.csv");
    print_undirected_network(_allNodes, data, afterCMIT);
    afterCMIT.close();
    for(int i = 0 ; i<data.Get_NumberOfNodes();i++)
    {
        _allNodes[i].SymmetryCorrection(_allNodes, data);
    }
    ofstream afterSymmetryCorrection;
    afterSymmetryCorrection.open(_name + "/" + "afterSymmetryCorrection.csv");
    print_undirected_network(_allNodes, data, afterSymmetryCorrection);
    afterSymmetryCorrection.close();
    for(int i = 0 ; i<data.Get_NumberOfNodes();i++)
    {
        _allNodes[i].SortNeighbor();
        _numberOfNeighbors.push_back(_allNodes[i].Get_neighbors().size());
    }

}




void BN::Set_arc(int fromNode, int toNode, bool setparent)
{
    // update the edge set, and update numberOfEdge.
    EDGE temp_edge;
    temp_edge.Set_edge(fromNode,toNode);
    _edges.push_back(temp_edge);

    // update the child/parent/ancester/descendent
    _allNodes[fromNode].UpdateChild(toNode);
    _allNodes[fromNode].UpdateDescendant(toNode,_allNodes[toNode].Get_descendant());
    if(setparent)
    {
        _allNodes[toNode].UpdateParent(fromNode);
    }
    _allNodes[toNode].UpdateAncester(fromNode,_allNodes[fromNode].Get_ancester());

    // pass the ancester of the "fromNode" to the descendents of the "toNode"
    // also pass the descendents of the "toNode" to the ancester of the "fromNode"
    set<int>::iterator it;
    for(it=_allNodes[fromNode].Get_ancester().begin(); it!=_allNodes[fromNode].Get_ancester().end(); it++ )
    {
        _allNodes[(*it)].UpdateDescendant(toNode,_allNodes[toNode].Get_descendant());
    }
    for(it=_allNodes[toNode].Get_descendant().begin(); it!=_allNodes[toNode].Get_descendant().end(); it++ )
    {
        _allNodes[(*it)].UpdateAncester(fromNode,_allNodes[fromNode].Get_ancester());
    }
}



void BN::BICscore()
{
    for(unsigned int i=0; i< _allNodes.size(); i++)
    {
        _score += _allNodes[i].Get_scoreContribution();
    }
}

void BN::Update(DATA& data,double temper, ostream& os)
{
    if(!(_oneNeb))
    {
        vector<int> candidates;
        vector<int> triplet(3);
        triplet[0]=-1;
        triplet[1]=-1;
        triplet[2]=-1;
        candidates = _MinOverOne(_numberOfNeighbors);

        if(candidates.empty()) // no possible v-structures
        {
            _candOneNeb = _SearchNoneZero(_numberOfNeighbors);
            if(_candOneNeb.empty())
            {
                _growing = false;
            }
            else if( _candOneNeb.size()>1 )
            {
                int NodeA = -1 , NodeB = -1;
                int Ind;
                double totalScore = 0.0;
                for(int i = 0; i< _candOneNeb.size();i++)
                {
                    _totalScoreOneNeb += _allNodes[_candOneNeb[i]].Get_scoreOfNeighbors(0);
                    _scoreOneNeb.push_back(_allNodes[_candOneNeb[i]].Get_scoreOfNeighbors(0));
                }
                double randomNum = static_cast<double>( rand() ) / RAND_MAX;
                while(randomNum >= 1.0)
                {
                    randomNum = static_cast<double>( rand() ) / RAND_MAX;
                }
                randomNum *= _totalScoreOneNeb;
                double cumScore = 0.0;
                for(int i = 0; i<_candOneNeb.size(); i++)
                {
                    cumScore += _scoreOneNeb[i];
                    if(randomNum < cumScore)
                    {
                        NodeA = _candOneNeb[i];
                        Ind=i;
                        break;
                    }
                }

                if(NodeA==-1)
                {
                    NodeA = _candOneNeb[_candOneNeb.size()-1];
                    Ind = _candOneNeb.size()-1;
                }

                NodeB = _allNodes[NodeA].Get_neighbors(0);
                _candOneNeb.erase(_candOneNeb.begin()+Ind); // remove node A
                _scoreOneNeb.erase(_scoreOneNeb.begin()+Ind);
                int tmpInd = _IndexSearch(NodeB); // remove node B
                _candOneNeb.erase(_candOneNeb.begin()+tmpInd);
                _scoreOneNeb.erase(_scoreOneNeb.begin()+tmpInd);
                _SelectDirection(NodeA,NodeB,data,temper, os);
                _oneNeb=true;
            }
            else
            {
                int NodeA,NodeB;
                NodeA = _candOneNeb[0];
                NodeB = _allNodes[NodeA].Get_neighbors(0);
                _SelectDirection(NodeA,NodeB,data,temper, os);
                _growing = false;
            }
        }
        else
        {
            if(candidates.size()==1)
            {
                triplet[1] = candidates[0];
            }
            else
            {
                double arandomNum = static_cast<double>( rand() ) / RAND_MAX;
                while(arandomNum >= 1.0)
                {
                    arandomNum = static_cast<double>( rand() ) / RAND_MAX;
                }
                triplet[1] = candidates[ min( static_cast<int>(arandomNum*candidates.size()), static_cast<int>(candidates.size())-1 ) ];
            }

            double randomNum = static_cast<double>( rand() ) / RAND_MAX;
            while(randomNum >= 1.0)
            {
                randomNum = static_cast<double>( rand() ) / RAND_MAX;
            }

            double totalScore = 0.0;
            for(int i=0; i< _allNodes[triplet[1]].Get_neighbors().size();i++)
            {
                totalScore += _allNodes[triplet[1]].Get_scoreOfNeighbors(i);
            }

            randomNum *= totalScore;

            double cumScore = 0.0;
            int tmpSampleInde = -1;
            for(int i=0; i<_allNodes[triplet[1]].Get_neighbors().size();i++)
            {
                cumScore += _allNodes[triplet[1]].Get_scoreOfNeighbors(i);
                if(randomNum < cumScore)
                {
                    tmpSampleInde = i;
                    break;
                }
            }
            if(tmpSampleInde == -1)
            {
                tmpSampleInde =  _allNodes[triplet[1]].Get_neighbors().size()-1;
            }
            triplet[0] = _allNodes[triplet[1]].Get_neighbors(tmpSampleInde);
            _allNodes[triplet[1]].ClearRank(tmpSampleInde,true);
            _numberOfNeighbors[triplet[1]]--;
            _allNodes[triplet[0]].ClearRank(triplet[1],false);
            _numberOfNeighbors[triplet[0]]--;

            vector<int> tmpCommonNeb(_allNodes[triplet[0]].Get_neighbors().size()+_allNodes[triplet[1]].Get_neighbors().size());
            vector<int> :: iterator it;
            it = set_intersection(_allNodes[triplet[0]].Get_neighbors().begin(),_allNodes[triplet[0]].Get_neighbors().end(),
                                _allNodes[triplet[1]].Get_neighbors().begin(),_allNodes[triplet[1]].Get_neighbors().end(),
                                tmpCommonNeb.begin());
            tmpCommonNeb.resize(it-tmpCommonNeb.begin());

            randomNum = static_cast<double>( rand() ) / RAND_MAX;
            while(randomNum >= 1.0)
            {
                randomNum = static_cast<double>( rand() ) / RAND_MAX;
            }
            //determine if the triplet forms a complete graph K_3
            if(tmpCommonNeb.empty()) //incomplete graph
            {
                double anotherTotalScore = 0.0;
                for(int i=0; i< _allNodes[triplet[1]].Get_neighbors().size();i++)
                {
                    anotherTotalScore += _allNodes[triplet[1]].Get_scoreOfNeighbors(i);
                }
                randomNum *= anotherTotalScore;

                cumScore = 0.0;
                tmpSampleInde = -1;
                for(int i=0; i<_allNodes[triplet[1]].Get_neighbors().size();i++)
                {
                    cumScore += _allNodes[triplet[1]].Get_scoreOfNeighbors(i);
                    if(randomNum < cumScore)
                    {
                        tmpSampleInde = i;
                        break;
                    }
                }
                if(tmpSampleInde == -1)
                {
                    tmpSampleInde =  _allNodes[triplet[1]].Get_neighbors().size()-1;
                }
                triplet[2] = _allNodes[triplet[1]].Get_neighbors(tmpSampleInde);
                _allNodes[triplet[1]].ClearRank(tmpSampleInde,true);
                _numberOfNeighbors[triplet[1]]--;
                _allNodes[triplet[2]].ClearRank(triplet[1],false);
                _numberOfNeighbors[triplet[2]]--;
                _TripletSelection(triplet,data,temper,false, os);
            }
            else //complete graph
            {
                triplet[2] = tmpCommonNeb[min(static_cast<int> (randomNum * tmpCommonNeb.size()), static_cast<int>(tmpCommonNeb.size())-1)];

                _allNodes[triplet[0]].ClearRank(triplet[2],false);

                _numberOfNeighbors[triplet[0]]--;
                _allNodes[triplet[1]].ClearRank(triplet[2],false);

                _numberOfNeighbors[triplet[1]]--;
                _allNodes[triplet[2]].ClearRank(triplet[0],false);

                _numberOfNeighbors[triplet[2]]--;
                _allNodes[triplet[2]].ClearRank(triplet[1],false);
                _numberOfNeighbors[triplet[2]]--;

                _TripletSelection(triplet,data,temper,true, os);
            }
        }
    }
    else // only nodes with one neighbor remain
    {
        if(_candOneNeb.empty())
        {
            _growing = false;
        }
        else if( _candOneNeb.size()>2 )
        {
            int NodeA = -1 , NodeB = -1;
            int Ind;
            double totalScore = 0.0;
            for(int i = 0; i< _candOneNeb.size();i++)
            {
                _totalScoreOneNeb += _allNodes[_candOneNeb[i]].Get_scoreOfNeighbors(0);
            }
            double randomNum = static_cast<double>( rand() ) / RAND_MAX;
            while(randomNum >= 1.0)
            {
                randomNum = static_cast<double>( rand() ) / RAND_MAX;
            }
            randomNum *= _totalScoreOneNeb;
            double cumScore = 0.0;
            for(int i = 0; i<_candOneNeb.size(); i++)
            {
                cumScore += _scoreOneNeb[i];
                if(randomNum < cumScore)
                {
                    NodeA = _candOneNeb[i];
                    Ind=i;
                    break;
                }
            }

            if(NodeA==-1)
            {
                NodeA = _candOneNeb[_candOneNeb.size()-1];
                Ind = _candOneNeb.size()-1;
            }

            NodeB = _allNodes[NodeA].Get_neighbors(0);
            _candOneNeb.erase(_candOneNeb.begin()+Ind);// remove node A
            _scoreOneNeb.erase(_scoreOneNeb.begin()+Ind);
            int tmpInd = _IndexSearch(NodeB);// remove node B
            _candOneNeb.erase(_candOneNeb.begin()+tmpInd);
            _scoreOneNeb.erase(_scoreOneNeb.begin()+tmpInd);

            _SelectDirection(NodeA,NodeB,data,temper, os);
        }
        else
        {
            int NodeA,NodeB;
            NodeA = _candOneNeb[0];
            NodeB = _allNodes[NodeA].Get_neighbors(0);

            _SelectDirection(NodeA,NodeB,data,temper, os);
            _growing = false;
        }
    }
}




void BN::HC(DATA& data)
{
    vector<int> sequence;
    int No_nodes = data.Get_NumberOfNodes();
	for(int i=0;i<No_nodes;i++)
	{
	 sequence.push_back(i);
	}

	for(int i=0;i<No_nodes;i++)
	{
		int j = rand()%(i+1);
		sequence[i]=sequence[j];
		sequence[j]=i;
	}

    for(int i=0; i<No_nodes;i++)
    {
        if(_changedAcc[i])
        {
            for(int j=0;j<No_nodes;j++) //sequence
            {
                if(i!= sequence[j])
                {
                    if(!(_ExistArc(sequence[j],i)) && !(_ExistArc(i,sequence[j])) ) //no arc
                    {
                        if(_IsDAG(sequence[j],i))
                        {
                            NODE tmpNode = _allNodes[i];
                            tmpNode.UpdateParent(sequence[j]);
                            tmpNode.UpdateCPD(data);
                            tmpNode.UpdateBIC(data);
                            if(tmpNode.Get_scoreContribution()>_allNodes[i].Get_scoreContribution())
                            {
                                _score -= _allNodes[i].Get_scoreContribution();
                                _allNodes[i] = tmpNode;
                                Set_arc(sequence[j],i,false);
                                _score += tmpNode.Get_scoreContribution();
                            }
                            else
                            {
                                _changedAcc[i] = false;
                            }
                        }
                    }
                }
            }
        }
    }
}


void BN::_SelectDirection(int NodeA, int NodeB, DATA& data,double temper, ostream& os)
{
    bool flagAtoB, flagBtoA;
    flagAtoB = _IsDAG(NodeA,NodeB);
    flagBtoA = _IsDAG(NodeB,NodeA);

    /* some cases:
    1) only A to B valid,
    2) only B to A valid,
    3) both are valid,
    4) neither is valid.
    */
    if(flagAtoB && (!flagBtoA) ) //only A to B valid,
    {
        double tempScore = _score;
        NODE tempNodeA = _allNodes[NodeA];
        NODE tempNodeB = _allNodes[NodeB];
        tempScore -= tempNodeA.Get_scoreContribution() + tempNodeB.Get_scoreContribution();
        tempNodeB.UpdateParent(NodeA);
        tempNodeA.UpdateChild(NodeB);
        tempNodeB.UpdateCPD(data);
        vector<double> nodeB_loglikelihoods = tempNodeB.UpdateBIC(data);
        tempNodeA.UpdateCPD(data);
        vector<double> nodeA_loglikelihoods = tempNodeA.UpdateBIC(data);
        vector<double> orig_nodeB_llk = _allNodes[NodeB].UpdateBIC(data);
        vector<double> orig_nodeA_llk = _allNodes[NodeA].UpdateBIC(data);
        tempScore += tempNodeA.Get_scoreContribution() + tempNodeB.Get_scoreContribution();
        double probSampleAtoB = 0.0;
        if(tempScore > _score)
        {
            probSampleAtoB = 1.0/(1.0+exp((_score-tempScore)/(temper)));
        }
        else
        {
            probSampleAtoB = exp((tempScore-_score)/(temper))/(1.0+exp((tempScore-_score)/(temper)));
        }
        double aRandomNumber = static_cast<double>(rand())/RAND_MAX;
        while(aRandomNumber >= 1.0)
        {
            aRandomNumber = static_cast<double>( rand() ) / RAND_MAX;
        }
        if(aRandomNumber <= probSampleAtoB ) // which means we do accept A to B
        {
            _allNodes[NodeA] = tempNodeA;
            _allNodes[NodeB] = tempNodeB;
            Set_arc(NodeA,NodeB,false);
            _score = tempScore;
        }
        std::ios cout_state(nullptr);
        cout_state.copyfmt(os);
        os << std::showpos << std::fixed << std::setprecision(4) << std::boolalpha;

        os << "Select Direction: " << data.Get_node_name(NodeA) << ", " << data.Get_node_name(NodeB) << endl;
        os << "\tBayesian network prior to this selection" << endl;
        os << "\t\t" << Write_Bnlearn_modelstring(data) << endl;
        vector<EDGE> edges = _edges;

        os << "\tConfiguration, is selected, probability of selection, prior log-likelihood, data BIC, posterior log-likelihood, loop introducing" << endl;
        os << "\t" << data.Get_node_name(NodeA) << " --- " << data.Get_node_name(NodeB) << "\t" << !(aRandomNumber<=probSampleAtoB) << "\t" << 1-probSampleAtoB << "\t";
        os << orig_nodeA_llk[0]+orig_nodeB_llk[0] << " {" << orig_nodeA_llk[0] << "," << orig_nodeB_llk[0] << "}\t";
        os << orig_nodeA_llk[1]+orig_nodeB_llk[1] << " {" << orig_nodeA_llk[1] << "," << orig_nodeB_llk[1] << "}\t";
        os << orig_nodeA_llk[2]+orig_nodeB_llk[2] << " {" << orig_nodeA_llk[2] << "," << orig_nodeB_llk[2] << "}\t";
        os << false << endl;
        os << "\t" << data.Get_node_name(NodeA) << " --> " << data.Get_node_name(NodeB) << "\t" << (aRandomNumber<=probSampleAtoB) << "\t" << probSampleAtoB << "\t";
        os << nodeA_loglikelihoods[0]+nodeB_loglikelihoods[0] << " {" << nodeA_loglikelihoods[0] << "," << nodeB_loglikelihoods[0] << "}\t";
        os << nodeA_loglikelihoods[1]+nodeB_loglikelihoods[1] << " {" << nodeA_loglikelihoods[1] << "," << nodeB_loglikelihoods[1] << "}\t";
        os << nodeA_loglikelihoods[2]+nodeB_loglikelihoods[2] << " {" << nodeA_loglikelihoods[2] << "," << nodeB_loglikelihoods[2] << "}\t";
        os << false << endl;
        os << "\t" << data.Get_node_name(NodeA) << " <-- " << data.Get_node_name(NodeB) << "\tIt will introduce a loop." << endl;
        os.copyfmt(cout_state);
    }
    else if( (!flagAtoB) && flagBtoA ) //only B to A valid,
    {
        double tempScore = _score;
        NODE tempNodeA = _allNodes[NodeA];
        NODE tempNodeB = _allNodes[NodeB];
        tempScore -= tempNodeA.Get_scoreContribution() + tempNodeB.Get_scoreContribution();
        tempNodeA.UpdateParent(NodeB);
        tempNodeA.UpdateCPD(data);
        vector<double> nodeA_loglikelihoods = tempNodeA.UpdateBIC(data);
        tempNodeB.UpdateChild(NodeA);
        tempNodeB.UpdateCPD(data);
        vector<double> nodeB_loglikelihoods = tempNodeB.UpdateBIC(data);
        vector<double> orig_nodeB_llk = _allNodes[NodeB].UpdateBIC(data);
        vector<double> orig_nodeA_llk = _allNodes[NodeA].UpdateBIC(data);
        tempScore += tempNodeA.Get_scoreContribution() + tempNodeB.Get_scoreContribution();
        double probSampleBtoA;
        if(tempScore > _score)
        {
            probSampleBtoA = 1.0/(1.0+exp((_score-tempScore)/(temper)));
        }
        else
        {
            probSampleBtoA = exp((tempScore-_score)/(temper))/(1.0+exp((tempScore-_score)/(temper)));
        }

        double aRandomNumber = static_cast<double>(rand())/RAND_MAX;
        while(aRandomNumber >= 1.0)
        {
            aRandomNumber = static_cast<double>( rand() ) / RAND_MAX;
        }
        if(aRandomNumber <= probSampleBtoA ) // which means we do accept B to A
        {
            _allNodes[NodeA] = tempNodeA;
            _allNodes[NodeB] = tempNodeB;
            Set_arc(NodeB,NodeA,false);
            _score = tempScore;
        }
            std::ios cout_state(nullptr);
            cout_state.copyfmt(os);
            os << std::showpos << std::fixed << std::setprecision(4) << std::boolalpha;

            os << "Select Direction: " << data.Get_node_name(NodeA) << ", " << data.Get_node_name(NodeB) << endl;
            os << "\tBayesian network prior to this selection" << endl;
            os << "\t\t" << Write_Bnlearn_modelstring(data) << endl;
            vector<EDGE> edges = _edges;

            os << "\tConfiguration, is selected, probability of selection, prior log-likelihood, data BIC, posterior log-likelihood, loop introducing" << endl;
            os << "\t" << data.Get_node_name(NodeA) << " --- " << data.Get_node_name(NodeB) << "\t" << !(aRandomNumber<=probSampleBtoA) << "\t" << 1-probSampleBtoA << "\t";
            os << orig_nodeA_llk[0]+orig_nodeB_llk[0] << " {" << orig_nodeA_llk[0] << "," << orig_nodeB_llk[0] << "}\t";
            os << orig_nodeA_llk[1]+orig_nodeB_llk[1] << " {" << orig_nodeA_llk[1] << "," << orig_nodeB_llk[1] << "}\t";
            os << orig_nodeA_llk[2]+orig_nodeB_llk[2] << " {" << orig_nodeA_llk[2] << "," << orig_nodeB_llk[2] << "}\t";
            os << false << endl;
            os << "\t" << data.Get_node_name(NodeA) << " --> " << data.Get_node_name(NodeB) << "\tIt will introduce a loop." << endl;
            os << "\t" << data.Get_node_name(NodeA) << " <-- " << data.Get_node_name(NodeB) << "\t" << (aRandomNumber<=probSampleBtoA) << "\t" << probSampleBtoA << "\t";
            os << nodeA_loglikelihoods[0]+nodeB_loglikelihoods[0] << " {" << nodeA_loglikelihoods[0] << "," << nodeB_loglikelihoods[0] << "}\t";
            os << nodeA_loglikelihoods[1]+nodeB_loglikelihoods[1] << " {" << nodeA_loglikelihoods[1] << "," << nodeB_loglikelihoods[1] << "}\t";
            os << nodeA_loglikelihoods[2]+nodeB_loglikelihoods[2] << " {" << nodeA_loglikelihoods[2] << "," << nodeB_loglikelihoods[2] << "}\t";
            os << false << endl;
            os.copyfmt(cout_state);
    }
    else if( flagAtoB && flagBtoA ) // both are valid,
    {
        // A to B part
        double tempScoreAtoB = _score;
        NODE tempNodeAAtoB = _allNodes[NodeA];
        NODE tempNodeBAtoB = _allNodes[NodeB];
        vector<double> orig_nodeA_llk = _allNodes[NodeA].UpdateBIC(data);
        vector<double> orig_nodeB_llk = _allNodes[NodeB].UpdateBIC(data);

        tempScoreAtoB -= tempNodeAAtoB.Get_scoreContribution() + tempNodeBAtoB.Get_scoreContribution();
        tempNodeBAtoB.UpdateParent(NodeA);
        tempNodeBAtoB.UpdateCPD(data);
        vector<double> B_AtoB_llks = tempNodeBAtoB.UpdateBIC(data);
        tempNodeAAtoB.UpdateChild(NodeB);
        tempNodeAAtoB.UpdateCPD(data);
        vector<double> A_AtoB_llks = tempNodeAAtoB.UpdateBIC(data);
        tempScoreAtoB += tempNodeAAtoB.Get_scoreContribution() + tempNodeBAtoB.Get_scoreContribution();
        
        // B to A part
        double tempScoreBtoA = _score;
        NODE tempNodeABtoA = _allNodes[NodeA];
        NODE tempNodeBBtoA = _allNodes[NodeB];
        tempScoreBtoA -= tempNodeABtoA.Get_scoreContribution() + tempNodeBBtoA.Get_scoreContribution();
        tempNodeABtoA.UpdateParent(NodeB);
        tempNodeABtoA.UpdateCPD(data);
        vector<double> A_BtoA_llks = tempNodeABtoA.UpdateBIC(data);
        tempNodeBBtoA.UpdateChild(NodeA);
        tempNodeBBtoA.UpdateCPD(data);
        vector<double> B_BtoA_llks = tempNodeBBtoA.UpdateBIC(data);
        tempScoreBtoA += tempNodeABtoA.Get_scoreContribution() + tempNodeBBtoA.Get_scoreContribution();
        
        // probabilities
        double probSampleAtoB,probSampleBtoA;
        if(tempScoreAtoB > tempScoreBtoA && tempScoreAtoB >_score)
        {
            probSampleAtoB = exp((tempScoreAtoB-_score)/(temper))/(1.0+exp((tempScoreAtoB-_score)/(temper))+exp((tempScoreBtoA-_score)/(temper)));
            probSampleBtoA = exp((tempScoreBtoA-_score)/(temper))/(1.0+exp((tempScoreAtoB-_score)/(temper))+exp((tempScoreBtoA-_score)/(temper)));
        }
        else if(tempScoreAtoB < tempScoreBtoA && tempScoreAtoB >_score)
        {
            probSampleAtoB = exp((tempScoreAtoB-_score)/(temper))/(1.0+exp((tempScoreAtoB-_score)/(temper))+exp((tempScoreBtoA-_score)/(temper)));
            probSampleBtoA = exp((tempScoreBtoA-_score)/(temper))/(1.0+exp((tempScoreAtoB-_score)/(temper))+exp((tempScoreBtoA-_score)/(temper)));
        }
        else
        {
            probSampleAtoB = exp((tempScoreAtoB-_score)/(temper))/(1.0+exp((tempScoreAtoB-_score)/(temper))+exp((tempScoreBtoA-_score)/(temper)));
            probSampleBtoA = exp((tempScoreBtoA-_score)/(temper))/(1.0+exp((tempScoreAtoB-_score)/(temper))+exp((tempScoreBtoA-_score)/(temper)));
        }
        
        // sample
        double aRandomNumber = static_cast<double>(rand())/RAND_MAX;
        while(aRandomNumber >= 1.0)
        {
            aRandomNumber = static_cast<double>( rand() ) / RAND_MAX;
        }
        int selected_direction = 0;
        if(aRandomNumber <= probSampleAtoB ) // which means we do accept A to B
        {
            selected_direction = 1;
            _allNodes[NodeA] = tempNodeAAtoB;
            _allNodes[NodeB] = tempNodeBAtoB;
            Set_arc(NodeA,NodeB,false);
            _score = tempScoreAtoB;
        }
        else if(aRandomNumber <= probSampleAtoB+probSampleBtoA) // accept B to A
        {
            selected_direction = 2;
            _allNodes[NodeA] = tempNodeABtoA;
            _allNodes[NodeB] = tempNodeBBtoA;
            Set_arc(NodeB,NodeA,false);
            _score = tempScoreBtoA;
        }
        std::ios cout_state(nullptr);
        cout_state.copyfmt(os);
        os << std::showpos << std::fixed << std::setprecision(4) << std::boolalpha;

        os << "Select Direction: " << data.Get_node_name(NodeA) << ", " << data.Get_node_name(NodeB) << endl;
        os << "\tBayesian network prior to this selection" << endl;
        os << "\t\t" << Write_Bnlearn_modelstring(data) << endl;
        vector<EDGE> edges = _edges;

        os << "\tConfiguration, is selected, probability of selection, prior log-likelihood, data BIC, posterior log-likelihood, loop introducing" << endl;
        os << "\t" << data.Get_node_name(NodeA) << " --- " << data.Get_node_name(NodeB) << "\t" << (selected_direction==0) << "\t" << 1-probSampleAtoB-probSampleBtoA << "\t";
        os << orig_nodeA_llk[0]+orig_nodeB_llk[0] << " {" << orig_nodeA_llk[0] << "," << orig_nodeB_llk[0] << "}\t";
        os << orig_nodeA_llk[1]+orig_nodeB_llk[1] << " {" << orig_nodeA_llk[1] << "," << orig_nodeB_llk[1] << "}\t";
        os << orig_nodeA_llk[2]+orig_nodeB_llk[2] << " {" << orig_nodeA_llk[2] << "," << orig_nodeB_llk[2] << "}\t";
        os << false << endl;
        os << "\t" << data.Get_node_name(NodeA) << " --> " << data.Get_node_name(NodeB) << "\t" << (selected_direction==1) << "\t" << probSampleAtoB << "\t";
        os << A_AtoB_llks[0]+B_AtoB_llks[0] << " {" << A_AtoB_llks[0] << "," << B_AtoB_llks[0] << "}\t";
        os << A_AtoB_llks[1]+B_AtoB_llks[1] << " {" << A_AtoB_llks[1] << "," << B_AtoB_llks[1] << "}\t";
        os << A_AtoB_llks[2]+B_AtoB_llks[2] << " {" << A_AtoB_llks[2] << "," << B_AtoB_llks[2] << "}\t";
        os << false << endl;
        os << "\t" << data.Get_node_name(NodeA) << " <-- " << data.Get_node_name(NodeB) << "\t" << (selected_direction==2) << "\t" << probSampleBtoA << "\t";
        os << A_BtoA_llks[0]+B_BtoA_llks[0] << " {" << A_BtoA_llks[0] << "," << B_BtoA_llks[0] << "}\t";
        os << A_BtoA_llks[1]+B_BtoA_llks[1] << " {" << A_BtoA_llks[1] << "," << B_BtoA_llks[1] << "}\t";
        os << A_BtoA_llks[2]+B_BtoA_llks[2] << " {" << A_BtoA_llks[2] << "," << B_BtoA_llks[2] << "}\t";
        os << false << endl;
        os.copyfmt(cout_state);
    }
    else
    {
        //dummy...
    }
}

string get_conf_string(vector<int>& candidateNodes, vector<int>& configuration, bool isCompleteGraph, DATA& data){
    string ret_og = "";
    string ret_afterloop = "";
    if(isCompleteGraph){
        vector<string> arrows = {" --- ", " --> ", " <-- "};

        ret_og += data.Get_node_name(candidateNodes[0]) + arrows[configuration[0]] + data.Get_node_name(candidateNodes[1]);
        ret_og += arrows[configuration[1]] + data.Get_node_name(candidateNodes[2]);
        ret_og += arrows[configuration[2]] + data.Get_node_name(candidateNodes[0]);

    } else {
        vector<string> arrows = {" --- ", " --> ", " <-- "};

        ret_og += data.Get_node_name(candidateNodes[0]) + arrows[configuration[0]] + data.Get_node_name(candidateNodes[1]);
        ret_og += arrows[configuration[1]] + data.Get_node_name(candidateNodes[2]);

    }
    return ret_og;
}

void BN::_TripletSelection(vector<int>& candidateNodes, DATA& data, double temper,bool isCompleteGraph, ostream& os)
{
    // if isCompleteGraph: candidateNodes[0] <--> candidateNodes[1] <--> candidateNodes[2] <--> candidateNodes[0]
    // if not: candidateNodes[0] <--> candidateNodes[1] <--> candidateNodes[2]
    if(isCompleteGraph){
        os << "Triplet Selection: " << data.Get_node_name(candidateNodes[0]) << " <--> " << data.Get_node_name(candidateNodes[1]) << " <--> " << data.Get_node_name(candidateNodes[2]);
        os << " <--> " << data.Get_node_name(candidateNodes[0]) << endl;
    } else {
        os << "Triplet Selection: " << data.Get_node_name(candidateNodes[0]) << " <--> " << data.Get_node_name(candidateNodes[1]) << " <--> " << data.Get_node_name(candidateNodes[2]) << endl;
    }
    
    double originalBIC = 0.0;
    vector<NODE> originalNodes(3);
    double TotalProb = 0.0;
    vector<bool> tmpAddedEdge(3);
    int MaxInd = 3;
    vector< vector<NODE> > candidates;
    vector<double> potentialBIC;
    vector<double> probOfSampling;
    vector<vector<int> > configurations;
    vector<int> tmpConfigurations(3);
    vector<bool> IsValidConfiguration; // A configuration is valid if: it doesn't introduce any loop 
    vector<string> configuration_string;
    vector< vector< vector<double> > > debug_BICs;
    for(int i=0; i<3;i++)
    {
        originalBIC += _allNodes[candidateNodes[i]].Get_scoreContribution();
        originalNodes[i] = _allNodes[candidateNodes[i]];
        tmpAddedEdge[i] = true;

    }

    if(isCompleteGraph)
    {
        MaxInd = 3;
        candidates.reserve(27);
        potentialBIC.reserve(27);
        probOfSampling.reserve(27);
        configurations.reserve(27);
        IsValidConfiguration.reserve(27);
        configuration_string.reserve(27);
        for(int i =0; i< 27;i++)
        {
            tmpConfigurations[0] = i%3;
            tmpConfigurations[1] = (i/3)%3;
            tmpConfigurations[2] = (i/9)%3;
            configurations.push_back(tmpConfigurations);
            candidates.push_back(originalNodes);
            potentialBIC.push_back(0.0);
            probOfSampling.push_back(0.0);
            configuration_string.push_back("NA");
            vector< vector<double> > tmpBICvector = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
            debug_BICs.push_back(tmpBICvector);
            IsValidConfiguration.push_back(true);
        }

    }
    else
    {
        MaxInd = 2;
        candidates.reserve(9);
        potentialBIC.reserve(9);
        probOfSampling.reserve(9);
        configurations.reserve(9);
        IsValidConfiguration.reserve(9);
        configuration_string.reserve(9);
        for(int i =0; i< 9;i++)
        {
            tmpConfigurations[0] = i%3;
            tmpConfigurations[1] = (i/3)%3;
            configurations.push_back(tmpConfigurations);
            candidates.push_back(originalNodes);
            potentialBIC.push_back(0.0);
            probOfSampling.push_back(0.0);
            configuration_string.push_back("NA");
            vector< vector<double> > tmpBICvector = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
            debug_BICs.push_back(tmpBICvector);
            IsValidConfiguration.push_back(true);
        }
    }

    for(int i=0 ; i< pow(3,MaxInd);i++)
    {
        for(int j=0;j<MaxInd;j++)
        {
            switch( configurations[i][j])
            {
                case 0:
                    break;
                case 1: // j -> j+1, fromNode= j%3, and toNode= (j+1)%3
                    if( candidates[i][j%3].Get_ancester().find(candidateNodes[(j+1)%3]) ==  candidates[i][j%3].Get_ancester().end() ) //if it is DAG
                    {

                        //update fromNode: descendent/ child, and if (j+2)%3 in ancester of fromNode add toNode to its descendents as well
                        candidates[i][j%3].UpdateChild(candidateNodes[(j+1)%3]);
                        candidates[i][j%3].UpdateDescendant(candidateNodes[(j+1)%3],candidates[i][(j+1)%3].Get_descendant() );
                        if(candidates[i][j%3].Get_ancester().find(candidateNodes[(j+2)%3]) !=  candidates[i][j%3].Get_ancester().end()) //(j+2)%3 is an ancester of j%3
                        {
                            candidates[i][(j+2)%3].UpdateDescendant(candidateNodes[(j+1)%3],candidates[i][(j+1)%3].Get_descendant() );
                        }
                        // update toNode
                        candidates[i][(j+1)%3].UpdateParent(candidateNodes[(j)%3]);
                        candidates[i][(j+1)%3].UpdateAncester(candidateNodes[(j)%3],candidates[i][(j)%3].Get_ancester() );
                        if(candidates[i][(j+1)%3].Get_descendant().find(candidateNodes[(j+2)%3]) !=  candidates[i][(j+1)%3].Get_descendant().end()) //(j+2)%3 is an ancester of j%3
                        {
                            candidates[i][(j+2)%3].UpdateAncester(candidateNodes[(j)%3],candidates[i][(j)%3].Get_ancester() );
                        }

                    }
                    else
                    {
                        IsValidConfiguration[i] = false;
                    }
                    break;
                case 2: // j<- j+1, fromNode= (j+1)%3, and toNode= j%3
                    if( candidates[i][(j+1)%3].Get_ancester().find(candidateNodes[(j)%3]) ==  candidates[i][(j+1)%3].Get_ancester().end() ) //if it is DAG
                    {
                        //update fromNode: descendent/ child, and if (j+2)%3 in ancester of fromNode add toNode to its descendents as well
                        candidates[i][(j+1)%3].UpdateChild(candidateNodes[(j)%3]);
                        candidates[i][(j+1)%3].UpdateDescendant(candidateNodes[(j)%3],candidates[i][(j)%3].Get_descendant() );
                        if(candidates[i][(j+1)%3].Get_ancester().find(candidateNodes[(j+2)%3]) !=  candidates[i][(j+1)%3].Get_ancester().end()) //(j+2)%3 is an ancester of j%3
                        {
                            candidates[i][(j+2)%3].UpdateDescendant(candidateNodes[(j)%3],candidates[i][(j)%3].Get_descendant() );
                        }
                        // update toNode
                        candidates[i][(j)%3].UpdateParent(candidateNodes[(j+1)%3]);
                        candidates[i][(j)%3].UpdateAncester(candidateNodes[(j+1)%3],candidates[i][(j+1)%3].Get_ancester() );
                        if(candidates[i][(j)%3].Get_descendant().find(candidateNodes[(j+2)%3]) !=  candidates[i][j%3].Get_descendant().end()) //(j+2)%3 is an ancester of j%3
                        {
                            candidates[i][(j+2)%3].UpdateAncester(candidateNodes[(j+1)%3],candidates[i][(j+1)%3].Get_ancester() );
                        }

                    }
                    else
                    {
                        IsValidConfiguration[i] = false;
                    }
                    break;
                default :
                cout<<"invalid type of connections, and you should not see this message unless there is something wrong."<<endl;
                break;
            }
        }
        if(IsValidConfiguration[i]){
            for(int j=0;j<3;j++)
            {
                candidates[i][j].UpdateCPD(data);
                vector<double> this_loglikelihoods = candidates[i][j].UpdateBIC(data);
                debug_BICs[i][j] = this_loglikelihoods;
                potentialBIC[i] += candidates[i][j].Get_scoreContribution();
            }
        }
    }
    
    // find the largest bic.
    double MaxBIC = numeric_limits<double>::lowest();
    
    for(int i=0; i<pow(3,MaxInd);i++)
    {
        if(potentialBIC[i]>MaxBIC && IsValidConfiguration[i])
        {
            MaxBIC = potentialBIC[i];
        }
    }
    
    for(int i=0; i<pow(3,MaxInd);i++)
    {
        if(IsValidConfiguration[i])
        {
            probOfSampling[i] = exp( (potentialBIC[i]-MaxBIC)/temper );
        }
        else
        {
            probOfSampling[i] = 0.0;
        }
        TotalProb += probOfSampling[i];
    }

    for(int i=0; i<pow(3,MaxInd);i++)
    {
        probOfSampling[i] /= TotalProb;
    }
    
    // Sample a random number
    double aRandomNumber = static_cast<double>( rand() ) / RAND_MAX;
    while(aRandomNumber > 1.0)
    {
        aRandomNumber = static_cast<double>( rand() ) / RAND_MAX;
    }
    double cumulativeProb = 0.0;
    int SampledEdge = -1;
    for(int i=0;i<pow(3,MaxInd);i++)
    {
        cumulativeProb += probOfSampling[i];
        if(aRandomNumber < cumulativeProb)
        {
            SampledEdge = i;
            break;
        }
    }

    // this part is to avoid numerical issue.
    if(SampledEdge == -1)
    {
        // TODO: what is this???
        for(int i=probOfSampling.size()-1;i>=0;i--)
        {
            if(probOfSampling[i]!=0.0)
            {
                SampledEdge = pow(3,MaxInd)-1;
            }
        }

    }

    os << "\tBayesian network prior to this selection" << endl;
    os << "\t\t" << Write_Bnlearn_modelstring(data) << endl;

    os << "\tConfiguration, is selected, probability of selection, prior log-likelihood, data BIC, posterior log-likelihood" << endl;
    os << "\t(Numbers in brackets are the individual values for each node: " << data.Get_node_name(candidateNodes[0]) << ", " << data.Get_node_name(candidateNodes[1]) << ", " << data.Get_node_name(candidateNodes[2]) << ")" << endl;
    for(int conf_idx=0; conf_idx!=pow(3, MaxInd); ++conf_idx){
        string conf_string = get_conf_string(candidateNodes, configurations[conf_idx], isCompleteGraph, data);
        bool is_selected = (conf_idx == SampledEdge);
        bool loop_introducing = !IsValidConfiguration[conf_idx];
        double prob_of_selection = probOfSampling[conf_idx];
        vector<double> prior_loglikelihoods = {debug_BICs[conf_idx][0][0], debug_BICs[conf_idx][1][0], debug_BICs[conf_idx][2][0]};
        double prior_loglikelihood = debug_BICs[conf_idx][0][0] + debug_BICs[conf_idx][1][0] + debug_BICs[conf_idx][2][0];
        vector<double> data_BICs = {debug_BICs[conf_idx][0][1], debug_BICs[conf_idx][1][1], debug_BICs[conf_idx][2][1]};
        double data_BIC = debug_BICs[conf_idx][0][1] + debug_BICs[conf_idx][1][1] + debug_BICs[conf_idx][2][1];
        vector<double> posterior_loglikelihoods = {debug_BICs[conf_idx][0][2], debug_BICs[conf_idx][1][2], debug_BICs[conf_idx][2][2]};
        double posterior_loglikelihood = debug_BICs[conf_idx][0][2] + debug_BICs[conf_idx][1][2] + debug_BICs[conf_idx][2][2];

        std::ios cout_state(nullptr);
        cout_state.copyfmt(os);
        os << std::showpos << std::fixed << std::setprecision(4) << std::boolalpha;
        os << "\t\t" << conf_string << "\t" << is_selected << "\t" << prob_of_selection << ", ";
        if(IsValidConfiguration[conf_idx]){
            os << prior_loglikelihood << " {" << prior_loglikelihoods[0] << "," << prior_loglikelihoods[1] << "," << prior_loglikelihoods[2] << "}\t";
            os << data_BIC << " {" << data_BICs[0] << "," << data_BICs[1] << "," << data_BICs[2] << "}\t";
            os << posterior_loglikelihood << " {" << posterior_loglikelihoods[0] << "," << posterior_loglikelihoods[1] << "," << posterior_loglikelihoods[2] << "}\t";
        } else {
            os << "This configuration is not valid. It will introduce a loop.\t";
        }
        os << endl;
        os.copyfmt(cout_state);
    }

    _score -= originalBIC;

    _score += potentialBIC[SampledEdge];

    // copy the nodes from temp to allnodes
    for(int j=0; j<3 ; j++)
    {
       // candidates[SampledEdge][j].Set_copyall(false);
        _allNodes[candidateNodes[j]].UpdateResult(candidates[SampledEdge][j]);
    }

    for(int j=0;j<MaxInd;j++)
        {
            switch( configurations[SampledEdge][j])
            {
                case 0:

                    break;
                case 1: // j -> j+1, fromNode= j%3, and toNode= (j+1)%3
                    if( _IsDAG(candidateNodes[j%3],candidateNodes[(j+1)%3]) ) //if it is DAG
                    {
                        Set_arc(candidateNodes[j%3],candidateNodes[(j+1)%3],true);
                    }
                    break;
                case 2: // j<- j+1, fromNode= (j+1)%3, and toNode= j%3
                    if(_IsDAG(candidateNodes[(j+1)%3],candidateNodes[(j)%3])) //if it is DAG
                    {
                        Set_arc(candidateNodes[(j+1)%3],candidateNodes[(j)%3],true);
                    }
                    break;
                default :
                    cout<<"invalid type of connections, and you should not see this message unless there is something wrong."<<endl;
                    break;
            }
        }

}


bool BN::_IsDAG(int fromNode, int toNode)
{
    // The method here is to see whether the "toNode" is an element of the ancester set of the "fromNode".
    bool DAG = true;
    if ( _allNodes[fromNode].Get_ancester().find( toNode ) != _allNodes[fromNode].Get_ancester().end() ) // if the iterator is different from the end means the "toNode" is found
        //from the ancester set of the "fromNode", then it is not a DAG.
    {
        DAG = false;
    }

    return DAG;
}


vector<int> BN::_SearchNoneZero(vector<int>& target)
{
    vector<int> index;
    for(unsigned int j=0; j<target.size();j++)
    {
        if(target[j] != 0)
        {
            index.push_back(j);
        }
    }
    return index;
}


vector< int > BN::_MinOverOne(vector<int>& target)// The algorithm is linear search.
{

    vector< int > index;
    int minimum;
    int ind;
    // the first loop is intended to detected if there is an at-least-two entry and start the search from the first at-least-two entry.
    for(unsigned int j=0; j<target.size();j++)
    {

        if(target[j]>1)
        {
            index.push_back(j);
            minimum = target[j];
            ind = j;
            break;
        }
        else
        {
            continue;
        }

    }

    // start the linear search.
    if(!index.empty())
    {
        for( unsigned int i=ind+1; i<target.size();i++)
        {
            if(target[i] == minimum)
            {
                index.push_back(i);
            }
            else if(target[i] < minimum && target[i]>1)
            {
                index.clear();// in case a smaller value came, clear the exisiting vector.
                index.push_back(i);
                minimum = target[i];
            }
            else
            {
                continue;
            }
        }
    }
    return index;
}


double BN::_Summation(vector<double>& target)// linear addition..
{
    double result = 0.0;

    for(unsigned int i=0; i< target.size();i++)
    {
        result += target[i];
    }
    return result;
}

int BN::_IndexSearch(int target)
{
    int L = 0;
    int U = _candOneNeb.size()-1;
    int tar = 0;
    while( (_candOneNeb[L]!=target) && (_candOneNeb[U]!=target) )
    {
        if(target<_candOneNeb[(L+U)/2])
        {
            U = (L+U)/2;
        }
        else
        {
            L = (L+U)/2;
        }
    }
    if (_candOneNeb[L]==target)
    {
        tar = L;
    }
    else if(_candOneNeb[U]==target)
    {
        tar = U;
    }
    else
    {
        cerr<<"Index search failed. It should not happen really."<<endl;
    }
    return tar;

}

bool BN::_ExistArc(int fromNode,int toNode)
{
    bool edge = false;

    if(_allNodes[fromNode].Get_child().find(toNode)!=_allNodes[fromNode].Get_child().end())
    {
        edge=true;
    }
    return edge;
}

string BN::Write_Bnlearn_modelstring(DATA& data)
{
    string model_string = "";
    int num_nodes = _allNodes.size();
    for(int to_node=0; to_node!=num_nodes; ++to_node){
        string to_node_name = data.Get_node_name(to_node);
        vector<int> parents = _allNodes[to_node].Get_parent();
        if(parents.empty()) model_string += "[" + to_node_name + "]";
        else {
            model_string += "[" + to_node_name + "|";
            for(int from_node_idx=0; from_node_idx!=parents.size(); ++ from_node_idx){
                int from_node = parents[from_node_idx];
                string from_node_name = data.Get_node_name(from_node);
                model_string += from_node_name;
                if(from_node_idx != parents.size()-1) model_string += ":";
            }
            model_string +="]";
        }
    }
    return model_string;
}