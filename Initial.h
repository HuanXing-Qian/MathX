#pragma once

#include <iostream>
#include <stack>
#include <queue>
#include <vector>
#include <unordered_map>
#include <string>
#include <random>
#include <complex>
#include<fstream>
#include<chrono>
#include<iomanip>
#include <cstdlib>
#include <map>
#include <regex>
#include <set>
#ifdef _WIN32
    #include <windows.h>  // 用于 Windows 平台
#else
    #include <unistd.h>   // 用于 Unix 平台
#endif

#define rep(i,a,b,l) for(auto i=(a);(l)>0?i<=(b):i>=(b);i+=(l))
#define RESET   "\033[0m"  
#define BLACK   "\033[30m"  
#define RED     "\033[31m"   
#define GREEN   "\033[32m"
#define YELLOW  "\033[33m"    
#define BLUE    "\033[34m"     
#define MAGENTA "\033[35m"  
#define CYAN    "\033[36m"
#define WHITE   "\033[37m"
#define BOLDBLACK   "\033[1m\033[30m" 
#define BOLDRED     "\033[1m\033[31m"
#define BOLDGREEN   "\033[1m\033[32m"
#define BOLDYELLOW  "\033[1m\033[33m"
#define BOLDBLUE    "\033[1m\033[34m"
#define BOLDMAGENTA "\033[1m\033[35m"
#define BOLDCYAN    "\033[1m\033[36m"
#define BOLDWHITE   "\033[1m\033[37m"
using namespace std;
const int INF = 0x3f3f3f3f;
using ll = long long;
using ld = long double;
using u64 = uint64_t;
using cld = complex<ld>;

constexpr ld PI = 3.14159265358979323846264338327950288L;

ll precision = 5;
int lagg = 1;
int baoliu = 4;
bool timing = true;
bool rpn = false;
bool high_precision = false;
bool science = false;
bool complex_mode = false;

const string func[100] = {
    "sqrt",
    "sin",
    "cos",
    "tan",
    "log",
    "ln",
    "exp",
    "exp2",
    "log10",
    "log2",
    "arcsin",
    "arccos",
    "arctan",
    "sinh",
    "cosh",
    "tanh",
    "ceil",
    "floor",
    "round",
    "abs",
    "erf",
    "Gamma",
    "sinc",
    "lngamma",
    "atr",
    "rta",
    "prime",
    "isprime",
    "min",
    "max",
    "arctan2"
    "gcd",
    "lcm",
    "hypot",
    "C",
    "rand",
};
int yuan[100];

const unordered_map<string, string> constants = {
    {"pi", "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899"},
    {"e", "2.71828182845904523536028747135266249775724709369995957496696762772407663035354759"},
    {"gamma", "0.57721566490153286060651209008240243104215933593992359880576723488486772677766467"},
    {"phi", "1.61803398874989484820458683436563811772030917980576286213544862270526046281890245"},
    {"sqrt2", "1.41421356237309504880168872420969807856967187537694807317667973799073247846210704"}
};

vector<string> little_knowledge;

set<int> knowledge_vis;

unordered_map<string, vector<pair<string, vector<string> > > > ledge;

unordered_map<string, string> variables;

unordered_map<string,cld>complex_variables;

struct node {
    string in;
    string result;
};

extern vector<node> history;

unordered_map<string, pair<string,string> > deffunc;

vector<pair<pair<string,ll>,pair<ll,int>>>Achievement;

vector<ll> primes;

inline void clear_screen() {
#ifdef _WIN32
    system("cls");
#else
    system("clear");
#endif
}

inline int prio(const string& c) {
    if (c == "+" || c == "-") return 1;
    if (c == "*" || c == "/") return 2;
    if (c == "^" || c == "%") return 3;
    if (find(begin(func), end(func), c) != end(func)) return 4;
    if (constants.find(c) != constants.end()) return 5;
    if (deffunc.find(c)!=deffunc.end())return 6;
    return 0;
}

inline string ns(string expr) {
    string result;
    int len = expr.length();
    rep(i,0,len-1,1) {
        char c = expr[i];
        if (c == ',') result += ' ';
        else result += c;

        if ((isalnum(c) || c == '.') && i < len - 1) {
            char next = expr[i + 1];
            if (string("+-*/^%()").find(next) != string::npos) result += ' ';
        }

        if (string("+*/^%()").find(c) != string::npos && i < len - 1) {
            if (expr[i + 1] != ' ') result += ' ';
        }

        if (c == '-' && i < len - 1) {
            bool isNegativeSign = i == 0 || string("+-*/^%( ").find(expr[i - 1]) != string::npos;
            int j = i + 1;
            string lin;
            while (j < len && isalpha(expr[j])) lin += expr[j++];
            if (prio(lin)==4) {
                result.pop_back();  // 移除已添加的 '-'
                result += "0 - ";  // 添加 "0-"
            } else {
                if (prio(lin) == 5) isNegativeSign = true;
                if (!isNegativeSign && expr[i + 1] != ' ') result += ' ';
            }
        }
    }

    string final;
    len = result.length();
    for (int i = 0; i < len; ) {
        if (result[i] == '/' && i + 1 < len) {
            int j = i + 1;
            while (j < len && result[j] == ' ') j++;

            if (j < len && (isdigit(result[j]) || result[j] == '.')) {
                int num_start = j;
                while (j < len && (isdigit(result[j]) || result[j] == '.')) j++;
                while (j < len && result[j] == ' ') j++;

                if (j < len && isalpha(result[j])) {
                    final += " / ( ";
                    final += result.substr(num_start, j - num_start);
                    final += " * ";
                    int ident_start = j;
                    while (j < len && isalpha(result[j])) j++;
                    final += result.substr(ident_start, j - ident_start);
                    final += " )";

                    i = j;
                    continue;
                }
            }
        }
        if (i + 1 < len && (isdigit(result[i]) || result[i] == '.') &&
            isalpha(result[i + 1])) {

            final += result[i];
            final += " * ";
            i++;
            continue;
            }

        if (i + 2 < len && result[i] == ')' &&
            isalpha(result[i + 2])) {
            final += result[i];
            final += " * ";
            i++;
            continue;
            }

        if (i + 2 < len && result[i] == ')' && result[i + 2] == '(') {
            final += result[i];
            final += " * ";
            i++;
            continue;
        }
		
        final += result[i];
        i++;
    }

    return final;
}

inline void out(const string& str, const string& color) {
    cout << color << str << RESET;
}

inline void inledge() {
    ledge.clear();
    string path;
    if (lagg == 1)path = "knowledge_en.txt";
    else if (lagg == 2)path = "knowledge_zh.txt";
    else if (lagg == 3)path = "knowledge_ru.txt";
    ifstream file(path);
    if (file) {
        int n;
        file >> n;
        while (n--) {
            string domain;
            int func_count;
            file >> domain >> func_count;
            while (func_count--) {
                string func_name, line;
                int line_count;
                file >> func_name >> line_count;
                file.ignore();
                if (!line_count)exit(0);
                vector<string> desc_lines;
                while(line_count--) {
                    getline(file, line);
                    desc_lines.push_back(line);
                }
                ledge[domain].emplace_back(func_name, desc_lines);
            }
        }
        file.close();
    }
}

inline void in_easter_egg()
{
    ifstream file("easter_egg.txt");
    if (file)
    {
        string line;
        while (getline(file,line))little_knowledge.push_back(line);
        file.close();
    }else cout<<"我的小知识文件呢口牙？"<<"\n";
}

inline void inAchievement()
{
    ifstream file("Achievement.txt");
    if (file)
    {
        string name;
        ll num,mu;
        int is;
        while (file>>name>>num>>mu>>is)
        {
            pair<string,ll> lin = {name,num};
            pair<ll,int> lin2 = {mu,is};
            Achievement.push_back({lin,lin2});
        }
    }else
    {
        cout<<"连这个也没了？？？"<<endl;
        return;
    }
    file.close();
    ifstream fe("Shi.txt");
    if (fe)
    {
        int x;
        while (fe>>x)knowledge_vis.insert(x);
    }
    fe.close();
}

inline void is_achievement() {
    for (auto& v : Achievement) { 
        if (v.first.second == v.second.first && v.second.second != 1) {
            cout << "获得成就“" << v.first.first << "”！\n";
            v.second.second = 1;
        }
    }

    //ofstream file("Achievement.txt");
    //if (file.is_open()) { 
    //    for (const auto& v : Achievement)file << v.first.first << " " << v.first.second << " " << v.second.first << " " << v.second.second << "\n";
    //    file.close();
    //} else cout << "无法打开文件以保存成就。\n";
}

bool parseComplex(const string& str, ld& real, ld& imag) {
    regex complex_re(R"(([-+]?\d*\.?\d+)([-+]\d*\.?\d*)i)");
    smatch match;
    if (regex_match(str, match, complex_re)) {
        real = stod(match[1]);
        imag = stod(match[2]);
        return true;
    }
    return false;
}

inline void Ready() {
    inledge();
    in_easter_egg();
    inAchievement();
    int x = 1;
    rep(i, 0, 99, 1) {
        if (func[i] == "min")x = 2;
        if (!func[i].empty())yuan[i] = x;
    }
    Achievement[0].first.second=1;
    variables["ans"]="0";
    complex_variables["ans"]=cld(0);
    primes.push_back(2);
    cout << "========================================\n"
         << "|            MathX 2.0 Beta            |\n"
         << "|      Enter [help] for commands       |\n"
         << "|          Example: 2sin(pi/4)         |\n"
         << "========================================\n";
}
