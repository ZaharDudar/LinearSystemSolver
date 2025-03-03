#include <LinearSystemSolver.h>

std::vector<int> getPerms(int r, std::vector<int> in){
    if (r == 1){
        return in;
    }
    std::vector<int> newVec(in.size()*2);
    int newL = newVec.size();
    for(int i=0; i<in.size(); i++){
        newVec[i*2] = in[i];
        newVec[i*2+1] = newL - 1 - in[i];
    }
    return getPerms(r-1, newVec);
} 