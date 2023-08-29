#ifndef __ADD_FILES__
#define __ADD_FILES__


#include "TString.h"
#include <vector>

using namespace std;

 vector<TString> addMultFiles(TString commaSeparatedList)
{
        // add multiple files as comma separated file list : "file1,file2,file3"
        vector<TString> fileArray;
        TObjArray* arr = commaSeparatedList.Tokenize(",");
        Int_t n = arr->GetEntries();
        if(n == 0) {
                //Warning("addMultFiles()","No file in input string '%s' found !",commaSeparatedList.Data());
        } else {
                for(Int_t i = 0; i < n; i++) {
                        TString name = ((TObjString*)arr->At(i))->GetString();
                        fileArray.push_back(name);
                        }
                }
        
        delete arr;
	return fileArray;
}


#endif
