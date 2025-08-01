// Developers: ShiYuze 
#pragma GCC optimize("O3,unroll-loops")
#include<bits/stdc++.h>
#include <Windows.h>
using namespace std;
#define rep(i,a,b,l) for(auto i=(a);(l)>0?i<=(b):i>=(b);i+=(l))
#define RESET   "\033[0m"       // RESET
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */
const int INF = 0x3f3f3f3f;
using ll = long long;
using ld = long double;

ll preci = 10;
int language = 1;
bool timing = 1;
bool RPN = 0;
bool high_precision = 0;

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
	"min",
	"max",
	"arctan2"
	"gcd",
	"lcm",
	"hypot" 
};
int yuan[100];

const unordered_map<string, string> constants = {
	{"pi", "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899"},
	{"e", "2.71828182845904523536028747135266249775724709369995957496696762772407663035354759"},
	{"gamma", "0.57721566490153286060651209008240243104215933593992359880576723488486772677766467"},
	{"phi", "1.61803398874989484820458683436563811772030917980576286213544862270526046281890245"},
	{"sqrt2", "1.41421356237309504880168872420969807856967187537694807317667973799073247846210704"}
};

unordered_map<string, vector<pair<string, vector<string> > > > ledge;

unordered_map<string, string> variables;

struct node {
	string in;
	string result;
};
vector<node> history;

int prio(string c) {
	if (c == "+" || c == "-") return 1;
	if (c == "*" || c == "/") return 2;
	if (c == "^" || c == "%") return 3;
	if (find(begin(func), end(func), c) != end(func)) return 4;
	if (constants.find(c) != constants.end()) return 5;
	return 0;
}

string NS(string expr) {
	string result;
	int len = expr.length();

	for (int i = 0; i < len; i++) {
		char c = expr[i];
		if (c == ',')result += ' ';
		else result += c;

		if ((isalnum(c) || c == '.') && i < len - 1) {
			char next = expr[i + 1];
			if (string("+-*/^()").find(next) != string::npos) {
				result += ' ';
			}
		}

		if (string("+*/^()").find(c) != string::npos && i < len - 1) {
			if (expr[i + 1] != ' ') {
				result += ' ';
			}
		}
		if (c == '-' && i < len - 1) {
			bool isNegativeSign = (i == 0) || (string("+-*/^( ").find(expr[i - 1]) != string::npos);
			int j = i + 1;
			string lin;
			while (isalpha(expr[j]) && j < len)lin += expr[j++];
			if (prio(lin) == 5)isNegativeSign = 1;
			if (!isNegativeSign && expr[i + 1] != ' ') {
				result += ' ';
			}
		}
	}

	return result;
}

void out(string str, string color) {
	cout << color << str << RESET;
}

ld calc(string str) {
	str = NS(str);
	stack<string> st;
	queue<string> q;
	ll i = 0;
	while (i < str.size()) {
		string token;
		while (i < str.size() && str[i] != ' ') {
			token += str[i];
			i++;
		}

		if (token.empty()) {
			i++;
			continue;
		}

		if (token == "setprecision") {
			token = "";
			i++;
			while (i < str.size() && str[i] != ' ') {
				token += str[i];
				i++;
			}
			if (token.empty()) {
				out("The usage of this function is to add a non-negative integer after the function name.\n", RED);
				return 0;
			}
			ld sum = stold(token);
			if (int(sum) != sum) {
				out("The usage of this function is to add a non-negative integer after the function name.\n", RED);
				return 0;
			}
			if (sum > 80) {
				out("Error: more than long double.\n", RED);
				return 0;
			}
			preci = sum;
			return 0;
		}

		if (variables.find(token) != variables.end()) q.push(variables[token]);
		else if (isdigit(token[0]) || (token[0] == '-' && isdigit(token[1])))q.push(token);
		else if (prio(token) == 4) st.push(token);
		else if (token == "^" || token == "%") {
			while (!st.empty() && st.top() != "(" && prio(token) < prio(st.top())) {
				q.push(st.top());
				st.pop();
			}
			st.push(token);
		} else if (token == "+" || token == "-" || token == "*" || token == "/") {
			while (!st.empty() && st.top() != "(" && prio(token) <= prio(st.top())) {
				q.push(st.top());
				st.pop();
			}
			st.push(token);
		} else if (prio(token) == 5) q.push(constants.at(token));
		else if (token[0] == '-' && prio(token.substr(1)) == 5) {
			string const_name = token.substr(1);
			q.push("-" + constants.at(const_name));
		} else if (token == "(") st.push(token);
		else if (token == ")") {
			while (!st.empty() && st.top() != "(") {
				q.push(st.top());
				st.pop();
			}
			if (!st.empty() && st.top() == "(") st.pop();
			if (!st.empty() && prio(st.top()) == 4) {
				q.push(st.top());
				st.pop();
			}
		} else {
			out("Error\n", RED);
			return 0;
		}
		i++;
	}

	while (!st.empty()) {
		q.push(st.top());
		st.pop();
	}

	if (RPN) {
		cout << "RPN: ";
		auto newq = q;
		while (!newq.empty()) {
			cout << newq.front() << ' ';
			newq.pop();
		}
		cout << "\n";
	}

	stack<ld> stk;
	while (!q.empty()) {
		string token = q.front();
		q.pop();
		if (isdigit(token[0]) || (token[0] == '-' && isdigit(token[1]))) stk.push(stold(token));
		else if (prio(token) < 4) {
			if (stk.size() < 2) {
				out("Error: Insufficient " + token + " operators.\n", RED);
				return 0;
			}
			ld r = stk.top();
			stk.pop();
			ld l = stk.top();
			stk.pop();
			ld result;
			if (token == "+") result = l + r;
			else if (token == "-") result = l - r;
			else if (token == "*") result = l * r;
			else if (token == "/") {
				if (!r) {
					out("inf\n", YELLOW);
					return 0;
				}
				result = l / r;
			} else if (token == "^") result = pow(l, r);
			else if (token == "%") result = fmod(l,r);
			stk.push(result);
		} else if (prio(token) == 4) {
			if (stk.empty()) {
				out("No " + token + " function.\n", RED);
				return 0;
			}
			ld res;
			int index = distance(begin(func), find(begin(func), end(func), token));
			if(yuan[index]==1){
				ld x = stk.top();
				stk.pop();
				if (token == "sin") res = sin(x);
				else if (token == "cos") res = cos(x);
				else if (token == "tan") res = tan(x);
				else if (token == "sqrt") res = sqrt(x);
				else if (token == "exp") res = exp(x);
				else if (token == "log" || token == "ln") res = log(x);
				else if (token == "exp2") res = exp2(x);
				else if (token == "log10") res = log10(x);
				else if (token == "log2") res = log2(x);
				else if (token == "arcsin") res = asin(x);
				else if (token == "arccos") res = acos(x);
				else if (token == "arctan") res = atan(x);
				else if (token == "sinh") res = sinh(x);
				else if (token == "cosh") res = cosh(x);
				else if (token == "tanh") res = tanh(x);
				else if (token == "ceil") res = ceil(x);
				else if (token == "floor") res = floor(x);
				else if (token == "round") res = round(x);
				else if (token == "abs") res = fabs(x);
				else if (token == "erf") res = erf(x);
				else if (token == "Gamma") res = tgamma(x);
				else if (token == "sinc") res = (x == 0) ? 1.0 : sin(x) / x;
				else if (token == "lngamma") res = lgamma(x);
				else if (token == "atr") res = x/180*stold(constants.at("pi"));
				else if (token == "rta") res = x/stold(constants.at("pi"))*180;
			}
			else if (yuan[index]==2) {
				ld x = stk.top();
				stk.pop();
				ld y = 0;
				if (!stk.empty())y = stk.top();
				if (!stk.empty()) {
					if (token == "min")res = min(x, y);
					else if (token == "max")res = max(x, y);
					else if (token == "arctan2")res = atan2(x, y);
					else if (token == "gcd"){
						if(int(x)!=x||int(y)!=y){
							out("Error: gcd requires integer inputs\n", RED);
        					return 0;
						} 
						res = __gcd(int(x),int(y));
					}else if (token == "lcm"){
						if(int(x)!=x||int(y)!=y){
							out("Error: lcm requires integer inputs\n", RED);
        					return 0;
						} 
						res = int(x*y)/__gcd(int(x),int(y));
					}else if (token == "hypot")res = sqrt(x*x+y*y);
					
					stk.pop();
				} else {
					out("Insufficient parameters.\n", RED);
					return 0;
				}
			}
			stk.push(res);
		}
	}
	return stk.top();
}

namespace high_precision_func {
	string add(string a,string b);
	string sub(string a, string b);
    string trim_zeros(const string &s) {
        auto pos = s.find_first_not_of('0');
        return (pos == string::npos) ? "0" : s.substr(pos);
    }
    string add(string a, string b) {
        bool a_neg = (a[0] == '-');
        bool b_neg = (b[0] == '-');
        if (a_neg && b_neg) return "-" + add(a.substr(1), b.substr(1));
        else if (a_neg) return sub(b, a.substr(1));
        else if (b_neg) return sub(a, b.substr(1));
        
        string res;
        int carry = 0;
        int i = a.size() - 1, j = b.size() - 1;
        while (i >= 0 || j >= 0 || carry) {
            int x = (i >= 0) ? a[i--] - '0' : 0;
            int y = (j >= 0) ? b[j--] - '0' : 0;
            int sum = x + y + carry;
            carry = sum / 10;
            res.push_back(sum % 10 + '0');
        }
        reverse(res.begin(), res.end());
        return trim_zeros(res);
    }
    string sub(string a, string b) {
        bool a_neg = (a[0] == '-');
        bool b_neg = (b[0] == '-');
        if (a_neg && b_neg) return sub(b.substr(1), a.substr(1)); 
        else if (a_neg) return "-" + add(a.substr(1), b);
        else if (b_neg) return add(a, b.substr(1)); 
        
        bool is_negative = false;
        if (a.size() < b.size() || (a.size() == b.size() && a < b)) {
            swap(a, b);
            is_negative = true;
        }

        string res;
        int borrow = 0;
        int i = a.size() - 1, j = b.size() - 1;
        while (i >= 0) {
            int x = a[i--] - '0' - borrow;
            int y = (j >= 0) ? b[j--] - '0' : 0;
            borrow = 0;
            if (x < y) {
                x += 10;
                borrow = 1;
            }
            res.push_back(x - y + '0');
        }
        reverse(res.begin(), res.end());
        res = trim_zeros(res);
        return is_negative ? "-" + res : res;
    }

    string div(string a, string b) {
        if (b == "0") throw runtime_error("Division by zero");
        bool neg = (a[0] == '-') ^ (b[0] == '-');
        string a_abs = (a[0] == '-') ? a.substr(1) : a;
        string b_abs = (b[0] == '-') ? b.substr(1) : b;

        if (a_abs == "0") return "0";
        if (b_abs == "1") return neg ? "-" + a_abs : a_abs;

        string quotient;
        string current;
        for (char c : a_abs) {
            current.push_back(c);
            current = trim_zeros(current);
            int count = 0;
            while (current.size() > b_abs.size() || 
                  (current.size() == b_abs.size() && current >= b_abs)) {
                current = sub(current, b_abs);
                count++;
            }
            quotient.push_back(count + '0');
        }
        quotient = trim_zeros(quotient);
        return neg ? "-" + quotient : quotient;
    }
    
    const double MM_PI = acos(-1.0);

	void parallel_fft(vector<complex<double>>& a, bool invert) {
	    int n = a.size();
	    if (n == 1) return;
	    vector<complex<double>> a0(n / 2), a1(n / 2);
	    for (int i = 0; i < n; i += 2) {
	        a0[i / 2] = a[i];
	        a1[i / 2] = a[i + 1];
	    }
	
	    parallel_fft(a0, invert);
	    parallel_fft(a1, invert);
	
	    double ang = 2 * MM_PI / n * (invert ? -1 : 1);
	    complex<double> w(1), wn(cos(ang), sin(ang));
	
	    for (int i = 0; i < n / 2; ++i) {
	        a[i] = a0[i] + w * a1[i];
	        a[i + n / 2] = a0[i] - w * a1[i];
	        if (invert) {
	            a[i] /= 2;
	            a[i + n / 2] /= 2;
	        }
	        w *= wn;
	    }
	}
	string mul(string num1, string num2) {
	    bool negative = false;
	    if (num1[0] == '-') {
	        negative = !negative;
	        num1 = num1.substr(1);
	    }
	    if (num2[0] == '-') {
	        negative = !negative;
	        num2 = num2.substr(1);
	    }
	
	    num1.erase(0, num1.find_first_not_of('0'));
	    num2.erase(0, num2.find_first_not_of('0'));
	    
	    if (num1.empty() || num2.empty()) return "0";
	    
	
	    int n = 1;
	    while (n < num1.size() + num2.size()) n <<= 1;
	    
	
	    vector<complex<double>> fa(n), fb(n);
	    for (int i = 0; i < num1.size(); ++i) fa[i] = num1[num1.size() - 1 - i] - '0';
	    
	    for (int i = 0; i < num2.size(); ++i) fb[i] = num2[num2.size() - 1 - i] - '0';
	    parallel_fft(fa, false);
	    parallel_fft(fb, false);
	    for (int i = 0; i < n; ++i) {
	        fa[i] *= fb[i];
	    }
	
	    parallel_fft(fa, true);
	
	    vector<int> result(n);
	    int carry = 0;
	    for (int i = 0; i < n; ++i) {
	        int value = (int)round(fa[i].real()) + carry;
	        result[i] = value % 10;
	        carry = value / 10;
	    }
	
	    int pos = n - 1;
	    while (pos > 0 && result[pos] == 0) {
	        --pos;
	    }
	
	    string res;
	    if (negative && !(pos == 0 && result[0] == 0)) {
	        res += '-';
	    }
	    for (int i = pos; i >= 0; --i) {
	        res += to_string(result[i]);
	    }
	
	    return res;
	}
    
    string pow(string base, string exponent) {
        string res = "1";
	    while (exponent != "0") {
	        if ((exponent.back() - '0') % 2 == 1) res = mul(res, base);
	        base = mul(base, base);
	        exponent = div(exponent, "2");
	    }
	    return res;
    }
}

string high_precision_calc(string str){
	str=NS(str);
	stack<string> st;
	queue<string> q;
	ll i=0;
	while (i < str.size()) {
		string token;
		while (i < str.size() && str[i] != ' ') {
			token += str[i];
			i++;
		}

		if (token.empty()) {
			i++;
			continue;
		}
		if (isdigit(token[0]) || (token[0] == '-' && isdigit(token[1])))q.push(token);
		else if (token == "^") {
			while (!st.empty() && st.top() != "(" && prio(token) < prio(st.top())) {
				q.push(st.top());
				st.pop();
			}
			st.push(token);
		}else if (token == "+" || token == "-" || token == "*" || token == "/") {
			while (!st.empty() && st.top() != "(" && prio(token) <= prio(st.top())) {
				q.push(st.top());
				st.pop();
			}
			st.push(token);
		}else if (token == "(") st.push(token);
		else if (token == ")") {
			while (!st.empty() && st.top() != "(") {
				q.push(st.top());
				st.pop();
			}
			if (!st.empty() && st.top() == "(") st.pop();
			if (!st.empty() && prio(st.top()) == 4) {
				q.push(st.top());
				st.pop();
			}
		} else {
			out("Error\n", RED);
			return 0;
		}
		i++;
	}
	
	while (!st.empty()) {
		q.push(st.top());
		st.pop();
	}
	
	stack<string> stk;
	
	while (!q.empty()) {
		string token = q.front();
		q.pop();
		if (isdigit(token[0]) || (token[0] == '-' && isdigit(token[1]))) stk.push(token);
		else if (prio(token) < 4) {
			if (stk.size() < 2) {
				out("Error: Insufficient " + token + " operators.\n", RED);
				return 0;
			}
			string r = stk.top();
			stk.pop();
			string l = stk.top();
			stk.pop();
			string result;
			if (token == "+") result = high_precision_func::add(l,r);
			else if (token == "-") result = high_precision_func::sub(l,r);
			else if (token == "*") result = high_precision_func::mul(l,r);
			else if (token == "/") {
				if (r=="0") {
					out("inf\n", YELLOW);
					return 0;
				}
				result = high_precision_func::div(l,r);
			} else if (token == "^") result = high_precision_func::pow(l, r);
			else{
				out("Error",RED);
				return 0;
			}
		}
	}
	return stk.top();
}

void inledge(){
	ledge.clear();
	string path; 
	if(language==1)path="knowledge_en.txt";
	else if(language==2)path="knowledge_zh.txt";
	else if(language==3)path="knowledge_ru.txt";
	ifstream file(path);
	if(file){
		int n;
		file>>n;
		while(n--){
			string domain;
			int func_count;
			file >> domain >> func_count;
			while (func_count--) {
				string func_name,line;
				int line_count; 
			    file >> func_name >> line_count;
			    file.ignore();
			    if(!line_count)exit(0);
			    vector<string> desc_lines;
			    for (int i = 0; i < line_count; i++) {
			        getline(file, line);
			        desc_lines.push_back(line);
			    }
			    ledge[domain].push_back({func_name, desc_lines});
			}
		}
		file.close();
	}
}

bool chao(string str){
	if (str.empty()) return 0;
	if (str == "exit") exit(0);
	if (str == "clear history") {
		out("Finished.\n", GREEN);
		history.clear();
		return 0;
	}

	if (str == "clear local history") {
		ofstream outFile("History.txt", ios::trunc);
		if (!outFile) {
			out("Error: Unable to clear history file.\n", RED);
			return 0;
		}
		outFile.close();
		out("Local history file (History.txt) has been cleared.\n", BLUE);
		return 0;
	}

	if (str.rfind("save history", 0) == 0) {
		string path;
		size_t pos = str.find_first_not_of(" \t", 13);

		if (pos == string::npos) path = "History.txt";
		else {
			path = str.substr(pos);
			if ((path.front() == '"' && path.back() == '"') ||
			    (path.front() == '\'' && path.back() == '\'')) {
				path = path.substr(1, path.size() - 2);
			}
		}
		ofstream outFile(path, ios::app);

		if (!outFile) {
			out("Error: Unable to open file '" + path + "' for writing.\n", RED);
			return 0;
		}

		time_t now = time(0);
		char* dt = ctime(&now);
		outFile << "\n=== History saved on: " << dt;

		for (size_t i = 0; i < history.size(); ++i) {
			outFile << i + 1 << ". " << history[i].in << " = "
			        << fixed << setprecision(preci) << history[i].result << endl;
		}

		outFile.close();
		out("History saved to '" + path + "' (" + to_string(history.size()) + " records)\n", GREEN);
		return 0;
	}

	if (str.rfind("load history", 0) == 0) {
		string path;
		size_t pos = str.find_first_not_of(" \t", 13);
		if (pos == string::npos) path = "History.txt";
		else {
			path = str.substr(pos);
			if ((path.front() == '"' && path.back() == '"') ||
			    (path.front() == '\'' && path.back() == '\'')) {
				path = path.substr(1, path.size() - 2);
			}
		}
		ifstream inFile(path);
		if (!inFile) {
			out("No history file found.\n", RED);
			return 0;
		}
		string line;
		while (getline(inFile, line)) cout << line << endl;

		inFile.close();
		return 0;
	}

	if (str[0] == '!') {
		str.erase(str.begin());
		if (str == "all") {
			if(!history.empty())rep(i, 0, (history.size() - 1), 1) cout << i + 1 << ": " << history[i].in << " = " << history[i].result << "\n";
			return 0;
		}
		int xv = calc(str);
		if (xv > history.size() || xv< 1) {
			out("Error\n", RED);
			cout << RED << "Error." << RESET << endl;
			return 0;
		}
		if (xv != 0)cout << history[xv - 1].in << " = " << history[xv - 1].result << "\n";
		else out("Error.\n", RED);
		return 0;
	}

	if (str == "help") {
	    cout << CYAN << "MATHX COMMAND REFERENCE\n"
	         << "=======================\n"
	         << "Calculation:\n"
	         << "  [expression]          Evaluate math expression\n"
	         << "  let [var] [value]     Define variable\n"
	         << "  list vars             List all variables\n"
	         << "  clear vars            Clear all variables\n\n"
	         
	         << "Precision Control:\n"
	         << "  setprecision [1-80]   Set decimal places\n"
	         << "  high precision on     Enable arbitrary-precision\n"
	         << "  high precision off    Disable arbitrary-precision\n\n"
	         
	         << "Display Modes:\n"
	         << "  rpn on                Show RPN translation\n"
	         << "  rpn off               Hide RPN\n"
	         << "  timing on             Show execution time\n"
	         << "  timing off            Hide timing\n\n"
	         
	         << "History System:\n"
	         << "  !n                    Recall nth history item\n"
	         << "  !all                  Show full history\n"
	         << "  clear history         Clear memory history\n"
	         << "  save history [file]   Save to file\n"
	         << "  load history [file]   Load from file\n\n"
	         
	         << "Learning System:\n"
	         << "  study [topic]         Show topic details\n"
	         << "  search [domain]       List domain functions\n\n"
	         
	         << "System Commands:\n"
	         << "  clear                 Clear screen\n"
	         << "  exit                  Quit program\n\n"
	         
	         << "SUPPORTED FUNCTIONS:\n"
	         << "Basic: + - * / ^ %\n"
	         << "Trig: sin cos tan asin acos atan\n"
	         << "Hyperbolic: sinh cosh tanh\n"
	         << "Log/Exp: log ln log10 exp\n"
	         << "Other: sqrt abs round gcd lcm\n\n"
	         
	         << "CONSTANTS:\n"
	         << "pi e gamma phi sqrt2\n"
	         << "=======================\n" << RESET;
	    return 0;
	}

	if (str.rfind("let ", 0) == 0) {
	    istringstream iss(str);
	    string cmd, var_name;
	    iss >> cmd >> var_name;
	    string value_str;
	    getline(iss, value_str);
	    value_str.erase(0, value_str.find_first_not_of(" \t"));
	    
	    if (!all_of(var_name.begin(), var_name.end(), [](char c) { return isalpha(c); })) {
	        out("Error: Invalid variable name. Use letters only.\n", RED);
	        return 0;
	    }
	    if (constants.find(value_str) != constants.end()) variables[var_name] = constants.at(value_str);
	    else if (value_str.find_first_not_of("0123456789.-") == string::npos) variables[var_name] = value_str;
	    
	    else {
	        if (high_precision) variables[var_name] = high_precision_calc(value_str);
	        else {
	            ld result = calc(value_str);
	            variables[var_name] = to_string(result);
			}
	    }
	    
	    out("Variable '" + var_name + "' set to " + variables[var_name] + "\n", BLUE);
	    return 0;
	}

	if (str == "list vars") {
		int x = 1;
		for (auto v : variables)
			cout << x << ". " << v.first << " = " << v.second << endl;
		return 0;
	}

	if (str.rfind("delete", 0) == 0) {
		variables.erase(str.substr(6));
		out("Finished.\n", GREEN);
		return 0;
	}

	if (str == "clear vars") {
		variables.clear();
		out("Finished\n", GREEN);
		return 0;
	}

	if (str == "clear window" || str == "clear") {
		system("cls");
		return 0;
	}

	if (str == "timing on" || str == "timing") {
		timing = 1;
		out("Opened timing.\n", CYAN);
		return 0;
	}

	if (str == "timing off") {
		timing = 0;
		out("Closed timing.\n", "\033[90m");
		return 0;
	}

	if (str == "rpn on" || str == "RPN on" || str == "rpn") {
		RPN = 1;
		out("Opened RPN.\n", CYAN);
		return 0;
	}
	
	if (str == "rpn off" || str == "RPN off") {
		RPN = 0;
		out("Closed RPN.\n", "\033[90m");
		return 0;
	}
	
	if (str == "high precision on" || str == "high precision" || str == "high"){
		high_precision = 1;
		out("Opened high precision.\nNote: Only addition (+), subtraction (-), multiplication (*), and division (/) can be used in high precision mode, and only integers can be input and output in this mode.\n", CYAN);
		return 0; 
	} 
	
	if (str == "high precision off" || str == "high off"){
		high_precision = 0;
		out("Closed high precision.\n", "\033[90m");
		return 0; 
	}
	
	if (str.rfind("language", 0) == 0) {
	    str.erase(0, 8);
	    str.erase(0, str.find_first_not_of(" \t"));
	    
	    if (str.empty()) {
	        out("Usage: language [en|zh|ru]\n", MAGENTA);
	        out("Example: language zh (切换为中文知识库)\n", MAGENTA);
	        return 0;
	    }
	
	    transform(str.begin(), str.end(), str.begin(), ::tolower);
	
	    if (str == "en") {
	        language = 1;
	        out("Knowledge base switched to English.\n", GREEN);
	    } else if (str == "zh") {
	        language = 2;
	        out("知识库已切换为中文。\n", GREEN);
	    } else if (str == "ru") {
	        language = 3;
	        out("База знаний переключена на русский.\n", GREEN);
		}
	    else out("Unsupported language. Available options: en, zh, ru\n", RED);
	    inledge();
	    return 0;
	}
	
	if (str.rfind("study", 0) == 0) {
	    string func_name = str.substr(5);
	    if(func_name.empty()) {
	        int x=1;
	        for(auto v:ledge) out(to_string(x++)+": "+v.first+"\n",CYAN);
	        return 0;
	    }
	    func_name.erase(func_name.begin());
	    string best_match;
	    int max_score = 0;
	    
	    for(auto &v: ledge) {
	        for(auto &u: v.second) {
	            const string &target = u.first;
	            int score = 0;
	            bool full_match = true;
	            if(!func_name.empty() && !target.empty() && 
	               tolower(func_name[0]) == tolower(target[0])) 
	                score += 2;
	            for(int i=0; i<max(func_name.size(), target.size()); i++) {
	                if(i < func_name.size() && i < target.size() && 
	                   tolower(func_name[i]) == tolower(target[i])) {
	                    score += 1;
	                } else if(i < min(func_name.size(), target.size())) {
	                    full_match = false;
	                }
	            }
	            
	            if(full_match&&u.first.size()==func_name.size()) {
	                for(auto t:u.second) out(t+"\n", CYAN);
	                return 0;
	            }
	            if(score > max_score) {
	                max_score = score;
	                best_match = target;
	            }
	        }
	    }
	    
	    if(max_score > func_name.size()/2) out("Did you mean \"" + best_match + "\"?\n", YELLOW);
	    else out("No matching concept found.\n", RED);
	    return 0; 
	}
	if(str.rfind("search ", 0) == 0){
		string func_name = str.substr(7);
		if(ledge.find(func_name)!=ledge.end()){
			int x=1;
			for(auto v:ledge[func_name]) out(to_string(x++)+": "+v.first+"\n",CYAN);
		}else out("No this one.\n", YELLOW);
		return 0;
	}
	return 1;
}

void Run() {
	cout << ">> " << flush;
	string str;
	getline(cin, str);

	if(!chao(str))return;
	
	bool f = 1;
	for (char c : str) {
		if (isalpha(c) || isdigit(c)) {
			f = 0;
			break;
		}
	}
	if (f)return;
	
	auto start = chrono::high_resolution_clock::now();	

    // 原来的计算逻辑
    if(!high_precision){
        ld result = calc(str);
        cout << "Result = " << fixed << setprecision(preci) << result << "\n";
        history.push_back({str, to_string(result)});
    } else {
        string result = high_precision_calc(str);
        cout << "Result = " << result << "\n";
        history.push_back({str, result});
    }
	

	auto end = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<std::chrono::microseconds>(end - start);

	if (timing)cout << "Use:" << duration.count() / 1000 << " ms." << "\n";
}



void Ready(){
	inledge();
	int x=1;
	rep(i,0,49,1){
		if(func[i]=="min")x=2;
		if(func[i]!="")yuan[i]=x;
	}
	cout << "========================================\n"
	     << "|           MathX 1.8 Beta             |\n"
	     << "|      Enter [help] for commands       |\n"
	     << "|       Example: 2 * sin(pi/4)         |\n"
	     << "========================================\n";
}

int main() {
	ios::sync_with_stdio(0);
	cin.tie(0);
	cout.tie(0);
	Ready();
	while (1) Run();
}
