#include <iostream> //header file library, lets us work with input and output objects
#include <string> //string library
#include <cmath> //math library
#include <fstream>// Stream class to both read and write from/to files.
#include <list>
#include <vector> 
#include <sstream> 




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

int main(){

    //Read all the files to make a tisse

 
   fstream centers;

   //lists associated with centers
   vector <int> centersloc{}; //creates a vector to add center indices
   vector <vector <double>>  centerscoord; //creates a vector to add center coordinates
   vector <int> centerspointreg{}; //creates a vector to add center point regions
   vector <vector <int>> centersregion{}; //creates a vector to add center regions

   centers.open("tissue225_3/centers.txt",ios::in); //open a file to perform read operation using file object
   if (centers.is_open()){ //checking whether the file is open
      string tp;
      while(getline(centers, tp)){ //read data from file object and put it into string.
       
        const char delim = ';';
        std::vector<std::string> out; 
        tokenize(tp, delim, out);
        centersloc.push_back(stoi(out.at(0)));
        centerscoord.push_back({stod(out.at(1)),stod(out.at(2))});
        centerspointreg.push_back(stoi(out.at(3)));
        string strout = out.at(4);
        strout.erase(strout.begin());
        vector <int> outlist{};
        extractIntegerWords(strout,outlist);
        //cout << out.at(4) << "\n";
        centersregion.push_back(outlist);
        
      }



      centers.close(); //close the file object.
   }

   fstream vertices;

   //lists associated with vertices
   vector <int> verticesloc{}; //creates a vector to add center indices
   vector <vector <double>>  verticescoord; //creates a vector to add center coordinates
   vector <vector <int>>  neighbourscent{}; //creates a vector to add center point regions
   vector <vector <int>>  neighboursvert{}; //creates a vector to add center regions

   vertices.open("tissue225_3/vertices.txt",ios::in); //open a file to perform read operation using file object
   if (vertices.is_open()){ //checking whether the file is open
      string tp;
      while(getline(vertices, tp)){ //read data from file object and put it into string.
       
        const char delim = ';';
        std::vector<std::string> out; 
        tokenize(tp, delim, out);
        verticesloc.push_back(stoi(out.at(0)));

        verticescoord.push_back({stod(out.at(1)),stod(out.at(2))});

        string strout1 = out.at(3);
        strout1.erase(strout1.begin());
        vector <int> outlist1{};
        extractIntegerWords(strout1,outlist1);
        //cout << out.at(4) << "\n";
        neighbourscent.push_back(outlist1);

        string strout2 = out.at(4);
        strout2.erase(strout2.begin());
        vector <int> outlist2{};
        extractIntegerWords(strout2,outlist2);
        //cout << out.at(4) << "\n";
        neighboursvert.push_back(outlist2);
        
      }

      vertices.close(); //close the file object.
   }



   fstream edges;

   //lists associated with edges
   vector <int> edgesloc{}; //creates a vector to add center indices
   vector <vector <int>>  edgeslist{}; //creates a list to add edges

   edges.open("tissue225_3/edges.txt",ios::in); //open a file to perform read operation using file object
   if (edges.is_open()){ //checking whether the file is open
      string tp;
      while(getline(edges, tp)){ //read data from file object and put it into string.
       
        const char delim = ';';
        vector <string> out; 
        tokenize(tp, delim, out);
        edgesloc.push_back(stoi(out.at(0)));

        string strout = out.at(1);
        strout.erase(strout.begin());
        vector <int> outlist{};
        //cout << outlist.empty() << endl;
        extractIntegerWords(strout,outlist);
        //cout << outlist.reserve() << endl;


        edgeslist.push_back(outlist);
        
      }

      edges.close(); //close the file object.
   }

   fstream boundaries;

   vector <int> boundarytissue{};
   vector <int> boundarywound{};
   vector <string> boundarieslist{};
   boundaries.open("tissue225_3/boundaries.txt",ios::in);
   if (edges.is_open()){ //checking whether the file is open
      string tp;
      while(getline(boundaries, tp)){ //read data from file object and put it into string.
        boundarieslist.push_back(tp);
      }

    string bt = boundarieslist[0];
    string bw = boundarieslist[1];

    bt.erase(bt.begin());
    bw.erase(bw.begin());

    extractIntegerWords(bt,boundarytissue);
    extractIntegerWords(bw,boundarywound);

    boundaries.close();

   }




    return 0;
}

   /*int x,y,myNum;
    myNum = 15;
    double myFloat = 5.99;
    char myLetter = 'D';
    string myText = "Hello";
    bool myBoolean = true;

    cout << "Hello World! \n";
    cout << "I'm learning c++ \n";
    cout << myFloat << endl;
    cout << "I am " << myNum << " years old. "; 
    
    
    / - Division
    % - Modulus - returns the remainder 
    ++ - increment
    -- - Decrement
    ^ - Power
    myText.length() - length of myText string
    myText.append(newText) - concatenates two strings together

    index starts with zero

if (condition1) {
  // block of code to be executed if condition1 is true
} else if (condition2) {
  // block of code to be executed if the condition1 is false and condition2 is true
} else {
  // block of code to be executed if the condition1 is false and condition2 is false
}

switch(expression) {
  case x:
    // code block
    break;
  case y:
    // code block
    break;
  default:
    // code block
}

while (condition) {
  // code block to be executed
}


do {
  // code block to be executed
}
while (condition);

 for (statement 1; statement 2; statement 3) {
  // code block to be executed
}
Statement 1 is executed (one time) before the execution of the code block.

Statement 2 defines the condition for executing the code block.

Statement 3 is executed (every time) after the code block has been executed.

    */
   //double const TSimul = 10;
   //double const K_run = 1;
   //double const G_run = 1;
   //int N;