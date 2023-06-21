#include <iostream> //header file library, lets us work with input and output objects
#include <string> //string library
#include <cmath> //math library
#include <fstream>// Stream class to both read and write from/to files.
#include <list>
#include <vector> 
#include <sstream> 
#include <algorithm>



using namespace std; //we can use names for objects and variables from the standard library

// from https://java2blog.com/split-string-space-cpp/


void tokenize(std::string const &str, const char delim, 
            std::vector<std::string> &out) 
{ 
    // construct a stream from the string 
    std::stringstream ss(str); 
 
    std::string s; 
    while (std::getline(ss, s, delim)) { 
        out.push_back(s); 
    } 
}


// from https://www.geeksforgeeks.org/extract-integers-string-c/

void extractIntegerWords(string str, vector <int> & lstr)
{
    stringstream ss;
 
    /* Storing the whole string into string stream */
    ss << str;
 
    /* Running loop till the end of the stream */
    string temp;
    int found;
    while (!ss.eof()) {
 
        /* extracting word by word from stream */
        ss >> temp;
 
        /* Checking the given word is integer or not */
        if (stringstream(temp) >> found)
            lstr.push_back(found);
            //cout << found << " ";
            
        /* To save from space at the end of string */
        temp = "";
    }
}

void norm(double x, double y){
    //return  x+y;
}

class Tissue{
    public:

    vector <int> centersloc{}; //creates a vector to add center indices
    vector <vector <double>>  centerscoord; //creates a vector to add center coordinates
    vector <int> centerspointreg{}; //creates a vector to add center point regions
    vector <vector <int>> centersregion{}; //creates a vector to add center regions

    vector <int> verticesloc{}; //creates a vector to add vertex indices
    vector <vector <double>>  verticescoord; //creates a vector to add vertex coordinates
    vector <vector <int>>  neighbourscent{}; //creates a vector to add neighbouring cells to a given vertex
    vector <vector <int>>  neighboursvert{}; //creates a vector to add neighbouring vertices to a given vertex

    vector <int> edgesloc{}; //creates a vector to add edge indices
    vector <vector <int>>  edgeslist{}; //creates a list to add edges

    vector <int> boundarytissue{}; //boundary of tissue
    vector <int> boundarywound{}; //boundary of wound

    void printboundary()
    {
        cout << "Tissue boundary ";
        for(auto i:boundarytissue){
            
            cout << i << " ";
        }

        cout << "\n";
        cout << "Wound boundary ";

        for(auto i:boundarywound){
            
            cout << i << " ";
        }

        cout << "\n";
    };

    void neighCenterV(int v){
        for(auto c:neighbourscent[v]){
            //cout << c << " ";
        }
        cout  << "\n";
    };

    void neighVertexV(int v){
        for(auto j:neighboursvert[v]){
            //cout << j << " ";
        }
        //cout  << "\n";
    };

    double AreaPoly(int c){
        vector <int> regionC = centersregion[c];
        double A = 0.0;
        double a = 0.0;
        for(int i = 0; i < regionC.size()-1; i++){
            double a1 = verticescoord[regionC[i]][0]*verticescoord[regionC[i+1]][1];
            double a2 = verticescoord[regionC[i]][1]*verticescoord[regionC[i+1]][0];
            a += a1-a2;
        }
        //cout << regionC.back() << "\n";
        //cout << regionC.front() << "\n";
        double aend1 = verticescoord[regionC.back()][0]*verticescoord[regionC.front()][1];
        double aend2 = verticescoord[regionC.front()][0]*verticescoord[regionC.back()][1];
    
        a += aend1 - aend2;
        //cout << sqrt(a*a) << "\n";
        A = sqrt(a*a)/2;
        return A;
        //cout << A << "\n";
    };


    double PerimeterPoly(int c){

        vector <int> regionC = centersregion[c];
        double P = 0.0;
        double p = 0.0;
        for(int i = 0; i < regionC.size()-1; i++){
            double a1 = verticescoord[regionC[i]][0] - verticescoord[regionC[i+1]][0];
            double a2 = verticescoord[regionC[i]][1] - verticescoord[regionC[i+1]][1];
            p += sqrt(a1*a1+a2*a2);
        }
        //cout << regionC.back() << "\n";
        //cout << regionC.front() << "\n";
        double aend1 = verticescoord[regionC.back()][0]-verticescoord[regionC.front()][0];
        double aend2 = verticescoord[regionC.front()][1]-verticescoord[regionC.back()][1];
    
        p += sqrt(aend1*aend1 + aend2*aend2);
        //cout << sqrt(a*a) << "\n";
        P = p;
        //cout << P << "\n";

        return P;

    };

    double LengthV(int v1, int v2){
        double x1 = verticescoord[v1][0] - verticescoord[v2][0];
        double y1 = verticescoord[v1][1] - verticescoord[v2][1];
        double l12 = sqrt(x1*x1+y1*y1);
        //cout << l12 << "\n";
        return l12;
    };

    int Ncells(){
        return centersloc.size();
    }

    int Nvertices(){
        return verticesloc.size();
    }

    vector <double> InitialArea(){
        int N = Ncells();
        vector <double> AreaVec{};
        for(int i = 0; i < N; i++){
            if(centersregion[i].size() > 2){
                AreaVec.push_back(AreaPoly(i));
            }
        }
        return AreaVec;

    };

    double meanArea(){
        int N = Ncells();
        double A0;
        double a = 0.0;
        for(int i = 0; i < N; i++){
            if(centersregion[i].size() > 2){
                a += AreaPoly(i);
            }
        }
        A0 = a/N;
        //cout << A0 << "\n";

        return A0;

    };




};

    double AreaOut(int c, vector <vector <int>> centersregion, vector <vector <double>> verticescoord){
        vector <int> regionC = centersregion[c];
        double A = 0.0;
        double a = 0.0;
        for(int i = 0; i < regionC.size()-1; i++){
            double a1 = verticescoord[regionC[i]][0]*verticescoord[regionC[i+1]][1];
            double a2 = verticescoord[regionC[i]][1]*verticescoord[regionC[i+1]][0];
            a += a1-a2;
        }
        //cout << regionC.back() << "\n";
        //cout << regionC.front() << "\n";
        double aend1 = verticescoord[regionC.back()][0]*verticescoord[regionC.front()][1];
        double aend2 = verticescoord[regionC.front()][0]*verticescoord[regionC.back()][1];
    
        a += aend1 - aend2;
        //cout << sqrt(a*a) << "\n";
        A = sqrt(a*a)/2;
        return A;
        //cout << A << "\n";
    };

        double AreaOutWound(int c, vector <vector <int>> centersregion, vector <vector <double>> verticescoord){
        vector <int> regionC = centersregion[c];
        double A = 0.0;
        double a = 0.0;
        for(int i = 0; i < regionC.size()-1; i++){
            double a1 = verticescoord[regionC[i]][0]*verticescoord[regionC[i+1]][1];
            cout << verticescoord[regionC[i]][0]<< verticescoord[regionC[i]][1] << "\n";
            double a2 = verticescoord[regionC[i]][1]*verticescoord[regionC[i+1]][0];
            a += a1-a2;
        }
        //cout << regionC.back() << "\n";
        //cout << regionC.front() << "\n";
        double aend1 = verticescoord[regionC.back()][0]*verticescoord[regionC.front()][1];
        double aend2 = verticescoord[regionC.front()][0]*verticescoord[regionC.back()][1];
    
        a += aend1 - aend2;
        //cout << sqrt(a*a) << "\n";
        A = sqrt(a*a)/2;
        return A;
        //cout << A << "\n";
    };


    double PerimeterOut(int c, vector <vector <int>> centersregion, vector <vector <double>> verticescoord){

        vector <int> regionC = centersregion[c];
        double P = 0.0;
        double p = 0.0;
        for(int i = 0; i < regionC.size()-1; i++){
            double a1 = verticescoord[regionC[i]][0] - verticescoord[regionC[i+1]][0];
            double a2 = verticescoord[regionC[i]][1] - verticescoord[regionC[i+1]][1];
            p += sqrt(a1*a1+a2*a2);
        }
        //cout << regionC.back() << "\n";
        //cout << regionC.front() << "\n";
        double aend1 = verticescoord[regionC.back()][0]-verticescoord[regionC.front()][0];
        double aend2 = verticescoord[regionC.front()][1]-verticescoord[regionC.back()][1];
    
        p += sqrt(aend1*aend1 + aend2*aend2);
        //cout << sqrt(a*a) << "\n";
        P = p;
        //cout << P << "\n";

        return P;

    };

double energy(vector <vector <int>> NeighC, vector <vector <int>> NeighV, vector <vector <int>> R, vector <vector <double>> vertices, int v,
double Kv, double Gv, double Lv, double Lwv, vector <double> A0, vector <int> woundbound, vector <int> boundary, int woundloc){

    vector <int> Nc = NeighC[v];
    vector <int> Nv = NeighV[v];
    double E = 0.0;

    for(auto c:Nc){
        double Pc = PerimeterOut(c, R, vertices);
        double Ac = AreaOut(c, R, vertices);
        double Ea = (Ac-A0[c]);
        E += 0.5*Kv*Ea*Ea + 0.5*Gv*Pc*Pc; 
    }

    vector <double> vertex1 = vertices[v];

    for(auto j:Nv){
        if(count(boundary.begin(),boundary.end(),v)==0){
            if(count(boundary.begin(),boundary.end(),j)==0){
                vector <double> vertex2 = vertices[j];  
                double xij = vertex1[0] - vertex2[0];
                double yij = vertex1[1] - vertex2[1];
                double lij = sqrt(xij*xij + yij*yij); 
                if(count(woundbound.begin(),woundbound.end(),v)){
                    if(count(woundbound.begin(),woundbound.end(),j)){
                        E+= Lwv*lij;
                    }
                }
                else{
                    E += -Lv*lij;
                }

        }

        }

    }

    return E;

};

//Calculate the force acting on a vertex through the finite gradient method
vector <double> force_vtx(vector <vector <int>> NeighC, vector <vector <int>> NeighV, vector <vector <int>> R, vector <vector <double>> vertices, int v,
double Kv, double Gv, double Lv, double Lwv, vector <double> A0, double dx, vector <int> woundbound, vector <int> boundary, int woundloc){

    //Create dummy vertex networks
    vector <vector <double>> vertices1x = vertices;
    vector <vector <double>> vertices2x = vertices;
    vector <vector <double>> vertices1y = vertices;
    vector <vector <double>> vertices2y = vertices;

    vertices1x[v][0] -= dx;
    vertices2x[v][0] += dx;
    vertices1y[v][1] -= dx;
    vertices2y[v][1] += dx;

    double E1x = energy(NeighC, NeighV,R, vertices1x,v, Kv, Gv,Lv, Lwv, A0, woundbound, boundary, woundloc);
    double E2x = energy(NeighC, NeighV,R, vertices2x,v, Kv, Gv,Lv, Lwv, A0, woundbound, boundary, woundloc);
    double E1y = energy(NeighC, NeighV,R, vertices1y,v, Kv, Gv,Lv, Lwv, A0, woundbound, boundary, woundloc);
    double E2y = energy(NeighC, NeighV,R, vertices2y,v, Kv, Gv,Lv, Lwv, A0, woundbound, boundary, woundloc);

    vector <double> f_v = {-0.5*(E2x-E1x)/dx,-0.5*(E2y-E1y)/dx};
    return f_v;
 
};


//Build an array of forces acting on all vertices of the network
vector <vector <double>> force_vtx_array(int Nvertices, double Kv, double Gv, double Lv, 
double Lwv, double r0, vector <double> A0, vector <vector <double>> cells, 
vector <vector <double>> vertices, int woundl, double dx, vector <int> boundary,
 vector <int> woundbound, vector <vector <int>> NeighC, vector <vector <int>> NeighV, vector <vector <int>> R){
    vector <vector <double>> Fv{};
    vector <double> f_v;
    for(int v=0; v < Nvertices; v++){
        if(count(boundary.begin(),boundary.end(),v)==0){
            f_v = force_vtx(NeighC, NeighV,R, vertices, v, Kv, Gv, Lv, Lwv, A0, dx, woundbound, boundary,woundl);
            for(int i = 0; i < NeighC[v].size(); i++){
                for(int j = i+1; j < NeighC[v].size(); j++){
                    vector <double> ci = cells[NeighC[v][i]];
                    vector <double> cj = cells[NeighC[v][j]];

                    double cxij = cj[0]-ci[0];
                    double cyij = cj[1]-ci[1];

                    double rij = sqrt(cxij*cxij + cyij*cyij);
                    
                    if(rij <= sqrt(r0)/2){
                        f_v[0] += 0.1*(rij-sqrt(r0)/2)*cxij/rij;
                        f_v[1] += 0.1*(rij-sqrt(r0)/2)*cyij/rij;
                    } 
                };
            };
        }
        else {
            f_v = {0.0, 0.0};
        };
        Fv.push_back(f_v);
    };
    return Fv;

};



int main(){

    //Read all the files to make a tisse

    Tissue tissue1;

 
   fstream centers;
   

   centers.open("tissue225_3/centers.txt",ios::in); //open a file to perform read operation using file object
   if (centers.is_open()){ //checking whether the file is open
      string tp;
      while(getline(centers, tp)){ //read data from file object and put it into string.
       
        const char delim = ';';
        std::vector<std::string> out; 
        tokenize(tp, delim, out);
        tissue1.centersloc.push_back(stoi(out.at(0)));
        tissue1.centerscoord.push_back({stod(out.at(1)),stod(out.at(2))});
        tissue1.centerspointreg.push_back(stoi(out.at(3)));
        string strout = out.at(4);
        strout.erase(strout.begin());
        vector <int> outlist{};
        extractIntegerWords(strout,outlist);
        //cout << out.at(4) << "\n";
        tissue1.centersregion.push_back(outlist);
        
      }



      centers.close(); //close the file object.
   }

   fstream vertices;

   vertices.open("tissue225_3/vertices.txt",ios::in); //open a file to perform read operation using file object
   if (vertices.is_open()){ //checking whether the file is open
      string tp;
      while(getline(vertices, tp)){ //read data from file object and put it into string.
       
        const char delim = ';';
        std::vector<std::string> out; 
        tokenize(tp, delim, out);
        tissue1.verticesloc.push_back(stoi(out.at(0)));

        tissue1.verticescoord.push_back({stod(out.at(1)),stod(out.at(2))});

        string strout1 = out.at(3);
        strout1.erase(strout1.begin());
        vector <int> outlist1{};
        extractIntegerWords(strout1,outlist1);
        //cout << out.at(4) << "\n";
        tissue1.neighbourscent.push_back(outlist1);

        string strout2 = out.at(4);
        strout2.erase(strout2.begin());
        vector <int> outlist2{};
        extractIntegerWords(strout2,outlist2);
        //cout << out.at(4) << "\n";
        tissue1.neighboursvert.push_back(outlist2);
        
      }

      vertices.close(); //close the file object.
   }



   fstream edges;

   edges.open("tissue225_3/edges.txt",ios::in); //open a file to perform read operation using file object
   if (edges.is_open()){ //checking whether the file is open
      string tp;
      while(getline(edges, tp)){ //read data from file object and put it into string.
       
        const char delim = ';';
        vector <string> out; 
        tokenize(tp, delim, out);
        tissue1.edgesloc.push_back(stoi(out.at(0)));

        string strout = out.at(1);
        strout.erase(strout.begin());
        vector <int> outlist{};
        //cout << outlist.empty() << endl;
        extractIntegerWords(strout,outlist);
        //cout << outlist.reserve() << endl;


        tissue1.edgeslist.push_back(outlist);
        
      }

      edges.close(); //close the file object.
   }

   fstream boundaries;

   vector <string> boundarieslist{};
   boundaries.open("tissue225_3/boundaries.txt",ios::in);
   if (boundaries.is_open()){ //checking whether the file is open
      string tp;
      while(getline(boundaries, tp)){ //read data from file object and put it into string.
        boundarieslist.push_back(tp);
      }

    string bt = boundarieslist[0];
    string bw = boundarieslist[1];

    bt.erase(bt.begin());
    bw.erase(bw.begin());

    extractIntegerWords(bt,tissue1.boundarytissue);
    extractIntegerWords(bw,tissue1.boundarywound);

    boundaries.close();

    

   }

   tissue1.meanArea();
   tissue1.LengthV(62,63);

   //Define some relevant constants

   const double K = 1.0;
   const double G = 1.0;
   const double mu = 1.0;
   const double DeltaL = 10;
   const double Tsim = 10;
   const vector <double> A0 = tissue1.InitialArea();

   double L = 5.0;
   double Lw = 4.0;
   int woundl = 70;
    


   //Define time and space steps


   double h = DeltaL/(2*sqrt(tissue1.Ncells()))*1e-4;
   double dt = mu/(K*tissue1.meanArea())*1e-2;

   cout << "space step -  " << h << " time step - " << dt << "\n";


   int Nsim = Tsim/dt;

   std::ofstream woundout;
   string stringArea = to_string(round(tissue1.meanArea()*1000)/1000);
   stringArea.erase(stringArea.size()-3);

   string stringL = to_string(round(L*100)/100);
   stringL.erase(stringL.size()-4);
   string stringLw =to_string(round(Lw*100)/100);
   stringLw.erase(stringLw.size()-4);

   woundout.open("tissue225_3/woundA"+stringArea+"L-"+stringL+"Lw-"+stringLw+".txt");

   if (woundout.is_open()){

   cout << Nsim << "\n";
   
   //woundout << to_string(0) << " " << round(AreaOut(woundl,tissue1.centersregion, tissue1.verticescoord)*1000)/1000 << " " << round(PerimeterOut(woundl,tissue1.centersregion, tissue1.verticescoord)*1000)/1000 << "\n";
   //cout << 62 << " " << tissue1.verticescoord[62][0] << " " << tissue1.verticescoord[62][1] << "\n";
   cout << tissue1.Nvertices() << "\n";
   for(int t = 0; t < 1; t++){
    vector <vector <double>> Fvertex = force_vtx_array(tissue1.Nvertices(),K,G,L,Lw,sqrt(tissue1.meanArea()),A0
    ,tissue1.centerscoord, tissue1.verticescoord,woundl, h,tissue1.boundarytissue, tissue1.boundarywound, tissue1.neighbourscent,
    tissue1.neighboursvert,tissue1.centersregion);
    for (int v = 0; v < tissue1.Nvertices(); v++){
        vector <double> Fv = Fvertex[v];

        //cout << v << " " << tissue1.verticescoord[v][0] << " " << tissue1.verticescoord[v][1] << "\n";

        tissue1.verticescoord[v][0] = tissue1.verticescoord[v][0]+mu*Fv[0]*dt;
        tissue1.verticescoord[v][1] = tissue1.verticescoord[v][1]+mu*Fv[1]*dt; 
        woundout << to_string(v) << " " << round(mu*Fv[0]*1000)/1000 << " " << round(mu*Fv[1]*1000)/1000 << "\n";
        //cout << v << " " << tissue1.verticescoord[v][0] << " " << tissue1.verticescoord[v][1] << "\n";

        //cout << v << " " << mu*Fv[0]*100*dt << " " << mu*Fv[1]*100*dt << "\n";
    }
    //woundout << to_string(t+1) << " " << round(AreaOut(woundl,tissue1.centersregion, tissue1.verticescoord)*1000)/1000 << " " << round(PerimeterOut(woundl,tissue1.centersregion, tissue1.verticescoord)*1000)/1000 << "\n";

    }
    //woundout.close();

   }
   else cout << "Unable to open file";
   //cout << 62 << " " << tissue1.verticescoord[62][0] << " " << tissue1.verticescoord[62][1] << "\n";








    return 0;
}
