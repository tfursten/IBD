#include "disk.h"

void Disk::initialize(double r)
{
    radius = r;
    rSq = r*r;
    halfRsq = rSq/2.0;
    getCellRange();
    vecDim = cellRange.size()-1;
    m_maxX = getMaxX();
    totalArea = M_PI*rSq;
    makeTables();

}

double Disk::integrate(double x1, double x2)
{
    //determine area of portion of the quarter disk in a column of lattice cells.  
    return (0.5 * sqrt(rSq-x2*x2) * x2 + rSq*0.5* asin(x2/radius)) - (0.5 * sqrt(rSq-x1*x1) * x1 + rSq*0.5* asin(x1/radius));
}

void Disk::getCellRange()
{
    //define borders of each lattice cell.
    cellRange.push_back(0.0);
    for(int i=0; i < radius; i++)
    {
        double x = i + 0.5;
        if (x<radius)
            cellRange.push_back(x);
    }
    cellRange.push_back(radius);
}

double Disk::getMaxX(){
    //Find column of cells that intesect the quarter circle at 45 degrees. 
    //Values for cells in any column after this point can be copied from 
    //cells in previous columns in mirrored locations.  
    double mX = pol2xy(M_PI/4.0).first;
    for (int i=vecDim; i>=1; i--){
        if((mX <= cellRange[i]) && (mX > cellRange[i-1]))
            return i;
    }
}

pair<double,double> Disk::pol2xy(double theta)
{
    //convert from polar to xy coordinates
    double x = radius*cos(theta);
    double y = radius*sin(theta);
    return make_pair(x,y);
}

double Disk::circle(double &x)
{
    //given the radius and x coord. solve for y coord. of circle
    return (double)sqrt(radius*radius-x*x);

}

int Disk::getBin(double x){
    for (int i=0; i<vecDim; i++){
        if(x <= cellRange[i+1])
            return i;
    }
}

void Disk::areas(double x1, double x2, int i){
    //find area for each cell in a column starting at the base
    //each full cell is subtracted from the total area.
    double A = integrate(x1,x2);
    for (int y = 0; y<vecDim; y++){
        double y1 = cellRange[y];
        double y2 = cellRange[y+1];
        double a = (x2-x1)*(y2-y1);
        if(a<A){
            //add area to current value (important to add for when columns are separated into two fragments)
            probMap[make_pair(i,y)] += a;
            A -= a;
        }
        else{
        probMap[make_pair(i,y)] += A;
        break;
        }
    }
}

void Disk::getAreas(int i){
    double x1 = cellRange[i];
    double x2 = cellRange[i+1];
    double fx1 = circle(x1);
    double fx2 = circle(x2);

    if(getBin(fx1) != getBin(fx2)){
    //if the arch of the circle crosses through the bottom of a cell, the area 
    //of that column will be calculated as two separate column fragments.  
        x2 = circle(cellRange[getBin(fx1)]);
        areas(x1,x2,i);
        x1 = x2;
        x2 = cellRange[i+1];
        fx1 = circle(x1);
        fx2 = circle(x2);
        }
    areas(x1,x2,i);
}




void Disk::makeTables(){
    //Calculates the area of lattice cells that are covered by a disk centered
    //on a focal lattice cell. The areas are then converted into a probability.  
    //This is simplified by only calculating values for a quarter slice of
    //disk and mirroring values to the other quadrants.  It is further 
    //simplified by only calculating approximately half of the values for 
    //the quarter slice and mirroring those as well.
    for(int i=0; i<m_maxX; i++)
        getAreas(i);
    //mirror values to the other half of the quarter disk
    for(int x=m_maxX; x<vecDim; x++){
        for(int y=0; y<vecDim; y++){
            if (probMap.count(make_pair(y,x)))
                probMap[make_pair(x,y)] = probMap[make_pair(y,x)];
        }
    }
    for(int x=0; x<vecDim; x++){
        for(int y=0; y<vecDim; y++){
            if (probMap.count(make_pair(x,y))){
                //convert area to a probability
                double prob = probMap[make_pair(x,y)]/totalArea;
                //mirror probability to other 3/4 of the disk
                probMap[make_pair(x,y)] = prob;
                probMap[make_pair(x*-1,y*-1)] += prob;
                probMap[make_pair(x*-1,y)] += prob;
                probMap[make_pair(x,y*-1)] += prob;
            }
        }
    }
    makeVectors();
    makeAliasTable();
}

xyCoord Disk::disperse(uint64_t u){
    return coordVec[xyTable(u)];
}

void Disk::makeVectors(){
    for (map<pair<int,int>,double>::iterator it = probMap.begin(); it != probMap.end(); ++it){
        coordVec.push_back(it->first);
        probVec.push_back(it->second);
    }
}

void Disk::makeAliasTable(){
    xyTable.create(probVec.begin(),probVec.end());
}

void Disk::printTables(){
    for (map<pair<int,int>,double>::iterator it = probMap.begin(); it != probMap.end(); ++it)
            cout << "X: " << it->first.first<< " Y: " << it->first.second << " Prob: " << it->second << endl;
}

