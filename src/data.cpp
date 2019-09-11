/*
Implementing of data.h

*/


#include "data.h"

DATA::DATA()
{
    //ctor
    _NumberOfNodes = 0;
    _NumberOfObservations = 0;
}

DATA::DATA(int obs, int nodes)
{
    _data.reserve(obs);
    _NumberOfNodes = nodes;
    _NumberOfObservations = obs;
    _param.reserve(nodes);
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
