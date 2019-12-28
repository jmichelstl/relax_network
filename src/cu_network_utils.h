#ifndef CU_NET_UTILS_H
#define CU_NET_UTILS_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <unordered_map>

#define PRIME1 (int) 1021297
#define PRIME2 (int) 1021651 //Primes for computing hash functions
#define FLOAT_TOL 1E-6 //Tolerance for deeming floating point numbers different
#define BIG_SLOPE 1E10
#define INF 1E12

using std::vector;
using std::string;
using std::ifstream;
using std::ofstream;
using std::cout;
using std::cin;

struct Point{
    Point() {
        this->tol = 0.0000001;
    }

    Point(double xval, double yval) : x(xval), y(yval){
        this->tol = 0.0000001;
    }

    Point(double xval, double yval, int i) : x(xval), y(yval), index(i){
        this->tol = 0.0000001;
    }

    double x,y;
    double tol;
    int index;

    bool operator == (const Point& point) const{
        return(abs(point.x - x) < tol && abs(point.y - y) < tol);
    }

    friend bool operator < (const Point &p1, const Point &p2){
        if(p2.y > p1.y + p1.tol) return true;
        else if(p1.y > p2.y + p2.tol) return false;
        else if(p2.x > p1.x + p1.tol) return true;
        else return false;
    }
};

struct Edge{
    Edge() {}

    Edge(Point start, Point end) : p1(start), p2(end){
        top = p1.y > p2.y ? p1.y : p2.y;
        bottom = p1.y < p2.y ? p1.y : p2.y;
        left = p1.x < p2.x ? p1.x : p2.x;
        right = p1.x > p2.x ? p1.x : p2.x;
    }

    Point p1, p2;
    double top, bottom, left, right;

    void reset(Point newp1, Point newp2){
        p1 = newp1;
        p2 = newp2;

        top = p1.y > p2.y ? p1.y : p2.y;
        bottom = p1.y < p2.y ? p1.y : p2.y;
        left = p1.x < p2.x ? p1.x : p2.x;
        right = p1.x > p2.x ? p1.x : p1.x;
    }

    bool operator == (const Edge& edge) const {
        return((edge.p1==p1 && edge.p2==p2)||(edge.p2==p1 && edge.p1==p2));
    }

    bool operator < (const Edge& edge) const {
        if(top > edge.top) return false;
        else if(top == edge.top && right > edge.right) return false;
        else return true;
    }

    double length(){
        return sqrt((p2.x - p1.x)*(p2.x - p1.x) + (p2.y - p1.y)*(p2.y - p1.y));
    }

    Point midpoint(){
        return Point((p1.x + p2.x) / 2, (p1.y + p2.y) / 2);
    }
};

namespace std {

    template<> struct hash<Point>{
        typedef size_t result_type;
        typedef Point argument_type;

        size_t operator() (const Point& p) const;
    };

    size_t hash<Point>::operator() (const Point& p) const {
        return (size_t) ((p.x+p.y)*(p.x+p.y)*(p.x*p.y+p.x*p.x+p.y*p.y*PRIME1))^PRIME2;
    }

    template<> struct hash<Edge>{
        typedef size_t result_type;
        typedef Edge argument_type;

        size_t operator() (const Edge& e) const;
    };

    size_t hash<Edge>::operator() (const Edge& e) const{
        size_t h1, h2;
        Point p1 = e.p1, p2 = e.p2;

        h1 = (size_t) ((p1.x+p1.y)*(p1.x+p1.y)*(p1.x*p1.y+p1.x*p1.x+p1.y*p1.y*PRIME1))^PRIME2;
        h2 = (size_t) ((p2.x+p2.y)*(p2.x+p2.y)*(p2.x*p2.y+p2.x*p2.x+p2.y*p2.y*PRIME1))^PRIME2;

        return h1*h1 + h2*h2;
    }
}


bool yesno(string message){
    string response;

    cout << message << "(y/n): ";

    while(true){
        getline(cin, response);
        if(!response.compare("y")) return true;
        else if(!response.compare("n")) return false;
        else printf("Respond with \"y\" or \"n\": ");
    }
}

vector<string> split(string input, char delim){
    vector<string> result;
    size_t start, iter;
    bool reading_token = false;

    for(iter = 0; iter < input.size(); iter++){
        if(delim != input[iter]){
            if(!reading_token){
                reading_token = true;
                start = iter;
            }
        }

        else{
            if(reading_token){
                reading_token = false;
                result.push_back(input.substr(start, iter - start));
            }
        }
    }

    if(reading_token) result.push_back(input.substr(start, iter - start));


    return result;
}

vector<double> parse_doubles(vector<string> numstring){
    double nextnum;
    vector<double> result;

    for(string s : numstring){
        if(sscanf(s.c_str(), "%lf", &nextnum)) result.push_back(nextnum);
    }

    return result;
}

vector<double> getdoubles(string prompt){
    vector<double> result;
    string response;

    printf("%s", prompt.c_str());
    getline(cin, response);

    return parse_doubles(split(response, ' '));
}

string enter_decline(string message){
    string response;
    cout << message << ", or return to decline: ";
    getline(cin, response);
    return response;
}

void open_dat_file(string prompt, ifstream& file_stream){
    string response;

    do{
        cout << prompt;
        getline(cin, response);
        file_stream.open(response);
        if(! file_stream.is_open()){
            if(! yesno("The file could not be read. Try again?")) break;
        }
    }while(! file_stream.is_open());
}

void open_output_file(string prompt, ofstream& file_stream){
    string response;

    do{
        cout << prompt;
        getline(cin, response);
        file_stream.open(response, ofstream::out);
        if(! file_stream.is_open()){
            if(! yesno("The file could not be read. Try again?")) break;
        }
    }while(! file_stream.is_open());
}

#endif
