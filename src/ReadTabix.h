/*
 * ReadTabix
 * Date: Aug-13-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#ifndef ReadTabix_h
#define ReadTabix_h

#include <string> 
#include <iostream> 
#include <stdlib.h>

#include "utils.h"

extern "C"{
#include "tabix.h"
}

using namespace std;

class ReadTabix{
 private:
    tabix_t * fpTab;
    ti_iter_t iteratorTab;

 public:
    ReadTabix(string file,string indexForFile,string chrName,int start,int end);
    ReadTabix(string file,string indexForFile,string chrName);

    ~ReadTabix();
    bool readLine(string & line);
    void repositionIterator(string chrName,int start,int end);
    void repositionIterator(string chrName);

    string getHeader();
    
};
#endif
