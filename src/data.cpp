/*
Implementing of data.h

*/


#include "data.h"
#include <string.h>

DATA::DATA()
{
    //ctor
    _NumberOfNodes = 0;
    _NumberOfObservations = 0;
    _hasPrior = false;
}

DATA::DATA(int obs, int nodes)
{
    _data.reserve(obs);
    _NumberOfNodes = nodes;
    _NumberOfObservations = obs;
    _param.reserve(nodes);
    
    // Fill the _prior with -1 (acting as no default value for this nodes)
    for(int idx=0; idx < nodes; ++idx){
        vector<double> init_priors(nodes, -1.0);
        _prior_prob.push_back(init_priors);
        _prior_pval.push_back(init_priors);    
    }
}

DATA::~DATA()
{
    //dtor
}

void DATA::ReadData(const char* fileName)
{
    ifstream inFile;
    inFile.open(fileName);
    typedef boost::tokenizer< boost::char_separator< char > > Tokenizer;
    boost::char_separator< char > sep(",");
    string line;

    if(inFile.fail())
    {
        cerr << "Error Opening data File" << endl;
        exit(1);
    }
    else
    {
        while( getline(inFile, line ))
        {
            Tokenizer info(line, sep);
            vector <int> val(_NumberOfNodes);
            int tempInd = 0;
            for (Tokenizer::iterator it = info.begin(); it != info.end(); ++it)
            {
                val[tempInd++] = strtod(it->c_str(),0);

            }
            _data.push_back(val);
        }
    }
    inFile.close();
}


void DATA::ReadParam(const char* fileName)
{
    ifstream inFile;
    inFile.open(fileName);
    int helper;
    if(inFile.fail())
    {
        cerr<< "unable to open the param file"<<endl;
        exit(1);
    }
    else
    {
        for(int i=0; i<_NumberOfNodes; i++)
        {
            inFile>> helper;
            _param.push_back(helper);

        }
    }
    inFile.close();
}


void DATA::ReadPrior(const char* fileName)
{
    ifstream inFile;
    inFile.open(fileName);
    typedef boost::tokenizer< boost::char_separator< char > > Tokenizer;
    boost::char_separator< char > sep(",");
    string line;

    if(inFile.fail())
    {
        cerr << "Error Opening prior File" << endl;
        exit(1);
    }
    else
    {
        while( getline(inFile, line ))
        {
            _hasPrior = true;
            Tokenizer info(line, sep);
            int from_node = stoi(*(info.begin())) - 1;
            int to_node = stoi( *(++info.begin()) ) - 1;
            double p_val = stod( *(++(++info.begin())) );
            double prob = stod( *(++(++(++info.begin()))) );
            cout << from_node << '\t' << to_node << '\t' << p_val << prob << endl;
            _prior_pval[from_node][to_node] = p_val;
            _prior_prob[from_node][to_node] = prob;
        }
    }
    inFile.close();
}
