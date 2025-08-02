// Developers: ShiYuze 
#pragma GCC optimize("O3,unroll-loops")
#include<bits/stdc++.h>
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
int lagg = 1;
int baoliu=4;
bool timing = 1;
bool RPN = 0;
bool high_precision = 0;
bool science = 0;

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
	"hypot",
	"sum" 
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
        if (c == ',') result += ' ';
        else result += c;
        

        if ((isalnum(c) || c == '.') && i < len - 1) {
            char next = expr[i + 1];
            if (string("+-*/^()").find(next) != string::npos) result += ' ';
        }

        if (string("+*/^()").find(c) != string::npos && i < len - 1) {
            if (expr[i + 1] != ' ') result += ' ';
        }
        
        if (c == '-' && i < len - 1) {
            bool isNegativeSign = (i == 0) || (string("+-*/^( ").find(expr[i - 1]) != string::npos);
            int j = i + 1;
            string lin;
            while (j < len && isalpha(expr[j])) lin += expr[j++];
            
            if (prio(lin) == 5) isNegativeSign = true;
            if (!isNegativeSign && expr[i + 1] != ' ') result += ' ';
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
        
        if(i+2<len&&result[i]==')'&&result[i+2]=='('){
        	final+=result[i];
        	final+=" * ";
        	i++;
        	continue;
		}
        
        final += result[i];
        i++;
    }

    return final;
}

void out(string str, string color) {
	cout << color << str << RESET;
}

ld calc(string str) {
	str = NS(str);
//	cout<<str<<endl; 
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
		        if(lagg==1) out("The usage of this function is to add a non-negative integer after the function name.\n", RED);
		        else if(lagg==2) out("用法：请在函数名后输入一个非负整数。\n", RED);
		        else if(lagg==3) out("Использование: укажите целое неотрицательное число после имени функции.\n", RED);
		        return 0;
		    }
		    ld sum = stold(token);
		    if (int(sum) != sum) {
		        if(lagg==1) out("The usage of this function is to add a non-negative integer after the function name.\n", RED);
		        else if(lagg==2) out("错误：必须输入整数。\n", RED);
		        else if(lagg==3) out("Ошибка: требуется целое число.\n", RED);
		        return 0;
		    }
		    if (sum > 80) {
		        if(lagg==1) out("Error: precision exceeds long double limit.\n", RED);
		        else if(lagg==2) out("错误：精度值超过long double限制。\n", RED);
		        else if(lagg==3) out("Ошибка: точность превышает пределы long double.\n", RED);
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
			if(lagg==1)out("Expression error.\n", RED);
			else if(lagg==2)out("伙计，表达式都写错了。\n",RED);
			else if(lagg==3)out("Ошибка выражения.\n", RED);
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
				if(lagg==1)out("Error: Insufficient " + token + " operators.\n", RED);
				else if(lagg==2)out("没有足够的"+token+"操作数，即有左没右或者反之。\n", RED);
				else if(lagg==3)out("Ошибка: не хватает операторов " + token + " для выполнения операции.\n", RED);
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
	    while (pos > 0 && result[pos] == 0) --pos;
	    
	
	    string res;
	    if (negative && !(pos == 0 && result[0] == 0)) res += '-';
	    
	    for (int i = pos; i >= 0; --i) res += to_string(result[i]);
	    
	
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
			if(lagg==1)out("Expression error.\n", RED);
			else if(lagg==2)out("伙计，表达式都写错了。\n",RED);
			else if(lagg==3)out("Ошибка выражения.\n", RED);
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
	
	stack<string> stk;
	
	while (!q.empty()) {
		string token = q.front();
		q.pop();
//		cout<<token<<endl;
		if (isdigit(token[0]) || (token[0] == '-' && token.length() > 1 && isdigit(token[1]))) stk.push(token);
		else if (prio(token) < 4) {
//			cout<<"Here"<<endl;
			if (stk.size() < 2) {
				out("Error: Insufficient " + token + " operators.\n", RED);
				return 0;
			}
			string r = stk.top();
			stk.pop();
			string l = stk.top();
			stk.pop();
			string result;
			if (token == "+")result = high_precision_func::add(l,r);
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
				if(lagg==1)out("Expression error.\n", RED);
				else if(lagg==2)out("伙计，表达式都写错了。\n",RED);
				else if(lagg==3)out("Ошибка выражения.\n", RED);
				return 0;
			}
			stk.push(result);
		}
//		cout<<stk.top()<<" "<<token<<endl;
	}
	return stk.top();
}

void inledge(){
	ledge.clear();
	string path; 
	if(lagg==1)path="knowledge_en.txt";
	else if(lagg==2)path="knowledge_zh.txt";
	else if(lagg==3)path="knowledge_ru.txt";
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
		if(lagg==1)out("Finished.\n", GREEN);
		else if(lagg==2)out("完毕。唉，但愿我不要再说类似“已完成”这样的话了。\n", GREEN);
		else if(lagg==3)out("Завершено.\n", GREEN); 
		history.clear();
		return 0;
	}
	
	if(str=="science off"){
		if(lagg==1)out("The Scientific notation mode has been turned off.\n", "\033[90m");
		else if(lagg==2)out("已关闭科学计数法模式。\n", "\033[90m");
		else if(lagg==3)out("Научный метод подсчёта закрыт.\n","\033[90m");
		science=0;
		return 0;
	}
	
	if(str.rfind("science",0)==0){
		str.erase(0,8);
		if(!str.empty()){
			bool f=1;
			for(char c :str)if(!isdigit(c)){f=0;break;}
			if(f){
				baoliu = stoi(str);
				if(lagg==1)out("The number of reserved digits is adjusted to "+to_string(baoliu)+"\n",CYAN);
				else if(lagg==2)out("保留位数调整为 "+to_string(baoliu)+"\n",CYAN);
				else if(lagg==3)out("Оговорка была изменена на"+to_string(baoliu)+"\n",CYAN);
			}
		} 
		if(lagg==1)out("The Scientific Notation mode has been enabled. Note: It can only be used in high-precision mode and will not display the full digits.\n",BLUE);
		else if(lagg==2)out("已开启科学计数法模式，注意：只能在高精度模式下使用，且会不显示完整数位。\n", BLUE);
		else if(lagg==3)out("Включен режим научного счёта, заметьте: его можно использовать только в режиме высокой точности, и он не будет отображать полный разряд.\n",BLUE);
		science=1;
		return 0;
	} 

	if (str == "clear local history") {
		ofstream outFile("History.txt", ios::trunc);
		if (!outFile) {
			if(lagg==1)out("Error: Unable to clear history file.\n", RED);
			else if(lagg==1)out("无法清除历史文件。\n",RED);
			else if(lagg==3)out("Не могу стереть исторические документы.\n",RED);
			return 0;
		}
		outFile.close();
		if(lagg==1)out("Local history file (History.txt) has been cleared.\n", BLUE);
		else if(lagg==2)out("记录已清除~\n", CYAN);
		else if(lagg==3)out("Локальный файл истории (History.txt) был очищен.\n", BLUE);
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
			if(lagg==1)out("Error: Unable to open file '" + path + "' for writing.\n", RED);
			else if(lagg==2)out("你小子又把文件放哪了？？？\n", RED);
			else if(lagg==3)out("Ошибка: Не удалось открыть файл '" + path + "' для записи.\n", RED);
			return 0;
		}

		time_t now = time(0);
		char* dt = ctime(&now);
		if(lagg==1)outFile << "\n=== History saved on: " << dt;
		else if(lagg==2)outFile << "\n===保存记录的时间：" << dt;
		else if(lagg==3)outFile << "\n===Время записи:" << dt;

		for (size_t i = 0; i < history.size(); ++i) {
			outFile << i + 1 << ". " << history[i].in << " = "
			        << fixed << setprecision(preci) << history[i].result << endl;
		}

		outFile.close();
		if(lagg==1)out("History saved to '" + path + "' (" + to_string(history.size()) + " records)\n", GREEN);
		else if(lagg==2)out("历史记录存于"+path+"，共"+to_string(history.size())+"条。\n",GREEN);
		else if(lagg==3)out("История сохранена в '" + path + "' (" + to_string(history.size()) + " записей)\n", GREEN);
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
			if(lagg==1)out("No history file found.\n", RED);
			else if(lagg==2)out("你文件放哪了？？？？？\n",RED);
			else if(lagg==3)out("Исторических документов нет.\n",RED);
			return 0;
		}
		string line;
		while (getline(inFile, line))cout << line << endl;

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
			if(lagg==1)out("Out of Range.\n", RED);
			else if(lagg==2)out("超出范围。\n",RED);
			else if(lagg==3)out("Вне диапазона.\n",RED);
			return 0;
		}
		if (xv != 0)cout << history[xv - 1].in << " = " << history[xv - 1].result << "\n";
		else {
			if(lagg==1)out("The expression is incorrect.\n", RED);
			else if(lagg==2)out("表达式错误。\n",RED);
			else if(lagg==3)out("Ошибка выражения.\n",RED);
		}
		return 0;
	}

	if (str == "help") {
		if(lagg==1){
			cout << CYAN << "MATHX COMMAND REFERENCE\n"
	         << "========================================\n"
	         << "BASIC CALCULATION\n"
	         << "  [expression]          Calculate any math expression\n"
	         << "  let x 5               Define variable 'x' with value 5\n"
	         << "  list vars             Show all variables\n"
	         << "  clear vars            Delete all variables\n\n"
	         
	         << "ADVANCED FEATURES\n"
	         << "  high precision on     Enable integer-only arbitrary precision\n"
	         << "  setprecision 15       Set decimal places (1-80)\n\n"
	         
	         << "DISPLAY SETTINGS\n"
	         << "  rpn on/off            Toggle Reverse Polish Notation mode\n"
	         << "  timing on/off         Show/hide calculation time\n"
	         << "  language en/zh/ru     Switch interface language\n\n"
	         
	         << "HISTORY SYSTEM\n"
	         << "  !3                    Recall 3rd calculation result\n"
	         << "  !all                  Show complete history\n"
	         << "  save history math.txt Export history to file\n\n"
	         
	         << "LEARNING MODE\n"
	         << "  study derivative      Learn about calculus derivatives\n"
	         << "  search trigonometry   List trig functions\n\n"
	         
	         << "SYSTEM\n"
	         << "  clear                 Clean the screen\n"
	         << "  exit                  Quit MathX\n"
	         << "========================================\n" << RESET;
		} else if(lagg==2){
			cout << CYAN << "MATHX 命令大全\n"
	         << "========================================\n"
	         << "基础计算\n"
	         << "  [表达式]          计算数学表达式\n"
	         << "  let x 5           定义变量x并赋值为5\n"
	         << "  list vars         显示所有变量\n"
	         << "  clear vars        清除全部变量\n\n"
	         
	         << "高级功能\n"
	         << "  high precision on 开启高精度整数模式\n"
	         << "  setprecision 15   设置小数位数(1-80)\n\n"
	         
	         << "显示设置\n"
	         << "  rpn on/off         切换逆波兰表达式模式\n"
	         << "  timing on/off      显示/隐藏计算时间\n"
	         << "  language en/zh/ru 切换界面语言\n\n"
	         
	         << "历史记录\n"
	         << "  !3                 调取第3条计算结果\n"
	         << "  !all               显示完整历史记录\n"
	         << "  save history math.txt 导出历史到文件\n\n"
	         
	         << "学习模式\n"
	         << "  study derivative   学习导数知识\n"
	         << "  search trigonometry 列出三角函数\n\n"
	         
	         << "系统命令\n"
	         << "  clear              清屏\n"
	         << "  exit               退出程序\n"
	         << "========================================\n" << RESET;
		} else if(lagg==3){
			cout << CYAN << "СПРАВКА ПО MATHX\n"
	         << "========================================\n"
	         << "ОСНОВНЫЕ ВЫЧИСЛЕНИЯ\n"
	         << "  [выражение]       Вычислить математическое выражение\n"
	         << "  let x 5           Определить переменную x = 5\n"
	         << "  list vars         Показать все переменные\n"
	         << "  clear vars        Удалить все переменные\n\n"
	         
	         << "ПРОДВИНУТЫЕ ФУНКЦИИ\n"
	         << "  high precision on Включить режим произвольной точности\n"
	         << "  setprecision 15    Установить знаки после запятой (1-80)\n\n"
	         
	         << "НАСТРОЙКИ ОТОБРАЖЕНИЯ\n"
	         << "  rpn on/off        Переключить ОПН (обратную польскую запись)\n"
	         << "  timing on/off     Показать/скрыть время вычислений\n"
	         << "  language en/zh/ru Сменить язык интерфейса\n\n"
	         
	         << "ИСТОРИЯ ВЫЧИСЛЕНИЙ\n"
	         << "  !3                Показать 3-й результат\n"
	         << "  !all              Показать всю историю\n"
	         << "  save history math.txt Экспорт истории в файл\n\n"
	         
	         << "ОБУЧАЮЩИЙ РЕЖИМ\n"
	         << "  study derivative  Изучить производные\n"
	         << "  search trigonometry Список тригонометрических функций\n\n"
	         
	         << "СИСТЕМА\n"
	         << "  clear             Очистить экран\n"
	         << "  exit              Выйти из программы\n"
	         << "========================================\n" << RESET;
		}
	    
	    return 0;
	}

	if (str.rfind("let ", 0) == 0) {
	    istringstream iss(str);
	    string cmd, var_name;
	    iss >> cmd >> var_name;
	    string value_str;
	    getline(iss, value_str);
	    value_str.erase(0, value_str.find_first_not_of(" \t"));
	
	    if(isdigit(var_name[0])){
	        if(lagg==1)out("The variable name cannot be used in Numbers at first.\n",YELLOW);
	        else if(lagg==2)out("变量名开头不能用数字。。\n",YELLOW);
	        else if(lagg==3)out("Переменные не могут начинаться с чисел.\n",YELLOW); 
	        return 0;
	    }
	
	    ld result = calc(value_str);
	    ostringstream oss;
	    oss << fixed << setprecision(80) << result;
	    variables[var_name] = oss.str();
	
	    if(lagg==1)out("Variable '" + var_name + "' set to " + variables[var_name] + "\n", BLUE);
	    else if(lagg==2)out("变量"+var_name+"已经设为"+variables[var_name]+"\n", BLUE);
	    else if(lagg==3)out("Variable '" + var_name + "' set to " + variables[var_name] + "\n", BLUE);
	    return 0;
	}

	if (str == "list vars") {
		int x = 1;
		for (auto v : variables)
			cout << x++ << ". " << v.first << " = " << v.second << endl;
		return 0;
	}

	if (str.rfind("delete", 0) == 0) {
	    str.erase(0, 6);
	    str.erase(0, str.find_first_not_of(" \t"));
	    str.erase(str.find_last_not_of(" \t") + 1);
	    
	    if (str.empty()) {
	        if (lagg == 1) out("Please specify a variable name to delete.\n", RED);
	        else if (lagg == 2) out("请指定要删除的变量名。\n", RED);
	        else if (lagg == 3) out("Укажите имя переменной для удаления.\n", RED);
	        return 0;
	    }
	    
	    if (variables.find(str) == variables.end()) {
	        if (lagg == 1) out("Variable '" + str + "' not found.\n", YELLOW);
	        else if (lagg == 2) out("变量 '" + str + "' 不存在。\n", YELLOW);
	        else if (lagg == 3) out("Переменная '" + str + "' не найдена.\n", YELLOW);
	        return 0;
	    }
	    
	    variables.erase(str);
	    if (lagg == 1) out("Deleted '" + str + "'\n", GREEN);
	    else if (lagg == 2) out("已删除变量 '" + str + "'\n", GREEN);
	    else if (lagg == 3) out("Переменная '" + str + "' удалена.\n", GREEN);
	    return 0;
	}

	if (str == "clear vars") {
		variables.clear();
		if(lagg==1)out("Finished\n", GREEN);
		else if(lagg==2)out("已完成。\n", GREEN);
		else if(lagg==3)out("Готово.\n",GREEN);
		return 0;
	}

	if (str == "clear window" || str == "clear") {
		system("cls");
		return 0;
	}

	if (str == "timing on" || str == "timing") {
		timing = 1;
		if(lagg==1)out("Enable timing.\n", CYAN);
		else if(lagg==2)out("已开启计时模式。\n", CYAN);
		else if(lagg==3)out("Таймер активирован.\n",CYAN); 
		return 0;
	}

	if (str == "timing off") {
		timing = 0;
		if(lagg==1)out("timing disabled.\n", "\033[90m");
		else if(lagg==2)out("已关闭计时模式。\n", "\033[90m");
		else if(lagg==3)out("Режим таймера выключен.\n", "\033[90m");
		return 0;
	}

	if (str == "rpn on" || str == "RPN on" || str == "rpn") {
		RPN = 1;
		if(lagg==1)out("Enable RPN.\n", CYAN);
		else if(lagg==2)out("已开启RPN模式。\n",CYAN);
		else if(lagg==3)out("RPN включён.\n",CYAN); 
		return 0;
	}
	
	if (str == "rpn off" || str == "RPN off") {
		RPN = 0;
		if(lagg==1)out("RPN disabled.\n", "\033[90m");
		else if(lagg==2)out("已关闭RPN模式。\n","\033[90m");
		else if(lagg==3)out("RPN отключён.\n","\033[90m");
		return 0;
	}
	
	if (str == "high precision on" || str == "high precision" || str == "high"){
		high_precision = 1;
		if(lagg==1)out("Opened high precision.\nNote: Only addition (+), subtraction (-), multiplication (*), and division (/), power(^)can be used in high precision mode, and only integers can be input and output in this mode.\n", CYAN);
		else if(lagg==2)out("已开启高精度模式。\n注意：高精度模式下仅支持加法（+）、减法（-）、乘法（*）和除法（/）、次方（^）运算，且该模式下仅允许输入和输出整数。\n",CYAN);
		else if(lagg==3)out("Активирован режим высокой точности.\nПримечание: В этом режиме доступны только операции сложения (+), вычитания (-), умножения (*) и деления (/),Степен (^), а ввод и вывод возможны только в целых числах.\n",CYAN);
		return 0; 
	} 
	
	if (str == "high precision off" || str == "high off"){
		high_precision = 0;
		if(lagg==1)out("Closed high precision.\n", "\033[90m");
		else if(lagg==2)out("已关闭高精度模式。\n", "\033[90m");
		else if(lagg==3)out("Режим высокой точности отключен.\n","\033[90m");
		return 0; 
	}
	
	if (str.rfind("language", 0) == 0) {
	    str.erase(0, 8);
	    str.erase(0, str.find_first_not_of(" \t"));
	    
	    if (str.empty()) {
	    	if(lagg==1){
	    		out("Usage: language [en|zh|ru]\n", MAGENTA);
	        	out("Example: language zh (Switch to Chinese)\n", MAGENTA);
			}else if(lagg==2){
				out("用法：language [en|zh|ru]\n",MAGENTA);
				out("比如：language ru (切换到俄文)\n",MAGENTA);
			}else if(lagg==3){
				out("Синтаксис: language [en|zh|ru]\n",MAGENTA);
				out("Например: language en (переключить на английский)\n",MAGENTA);
			} 
	        return 0;
	    }
	
	    transform(str.begin(), str.end(), str.begin(), ::tolower);
	
	    if (str == "en") {
	        lagg = 1;
	        out("It has been switched to English.\n", GREEN);
	    } else if (str == "zh") {
	        lagg = 2;
	        out("已切换为中文。\n", GREEN);
	    } else if (str == "ru") {
	        lagg = 3;
	        out("Переключился на русский.\n", GREEN);
		}
	    else {
	    	if(lagg==1)out("Unsupported language. Available options: en, zh, ru\n", RED);
	    	else if(lagg==2) out("这是不支持的语言。可选选项：en, zh, ru",RED); 
	    	else if(lagg==3) out("Неподдерживаемый язык. Доступные варианты: en, zh, ru.\n",RED);
		}
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
	            if(!func_name.empty() && !target.empty() && func_name[0] == target[0])
	                score += 2;
	            for(int i=0; i<max(func_name.size(), target.size()); i++) {
	                if(i < func_name.size() && i < target.size() && 
	                   func_name[i] == target[i]) {
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
	    
	    if(max_score > func_name.size()/2) {
	    	if(lagg==1)out("Did you mean \"" + best_match + "\"?\n", YELLOW);
	    	else if(lagg==2)out("您想说的是\""+best_match+"\"？\n", YELLOW);
	    	else if(lagg==3)out("Вы хотели сказать\""+best_match+"\"?\n", YELLOW);
		}
	    else {
	    	if(lagg==1)out("No matching concept found.\n", RED);
	    	else if(lagg=2)out("未找到匹配项。\n",RED);
	    	else if(lagg=2)out("Совпадений не найдено.\n",RED);
		}
	    return 0; 
	}
	if(str.rfind("search ", 0) == 0){
		string func_name = str.substr(7);
		if(ledge.find(func_name)!=ledge.end()){
			int x=1;
			for(auto v:ledge[func_name]) out(to_string(x++)+": "+v.first+"\n",CYAN);
		}else {
			if(lagg==1)out("No matching concept found.\n", RED);
	    	else if(lagg=2)out("未找到匹配项。\n",RED);
	    	else if(lagg=2)out("Совпадений не найдено.\n",RED);
		}
		return 0;
	}
//	还在测试： 
//	if(str.rfind("sum",0)==0){
//		str.erase(0,2);
//		if(str.empty()){
//			
//		}
//	} 
	
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

    if(!high_precision){
        ld result = calc(str);
        auto end = chrono::high_resolution_clock::now();
		auto duration = chrono::duration_cast<std::chrono::microseconds>(end - start);
        if(lagg==1)cout << "Result = ";
        else if(lagg==2)cout<<"结果 = ";
        else if(lagg==3)cout<<"Результат = "; 
        cout << fixed << setprecision(preci) << result << "\n";
        history.push_back({str, to_string(result)});
        if (timing)cout << "Use:" << duration.count() / 1000 << " ms." << "\n";
    } else {
    	string result = high_precision_calc(str);
    	auto end = chrono::high_resolution_clock::now();
		auto duration = chrono::duration_cast<std::chrono::microseconds>(end - start);
		if(lagg==1)cout << "Result = ";
	    else if(lagg==2)cout<<"结果 = ";
	    else if(lagg==3)cout<<"Результат = "; 
    	if(science){
			if(result.size()<baoliu){
				cout<<result<<endl;
				history.push_back({str, result});
			} else{
				string jia;
				jia+=result[0]+".";
				rep(i,1,(baoliu-1),1)jia+=result[i];
				jia+="e"+to_string(int(result.size()-1));
				history.push_back({str,jia});
				cout<<jia<<"\n";
			}
			
		}else{
			cout << result << "\n";
        	history.push_back({str, result});
		}
		
        if (timing)cout << "Use:" << duration.count() / 1000 << " ms." << "\n";
    }
}

void Ready(){
	inledge();
	int x=1;
	rep(i,0,49,1){
		if(func[i]=="min")x=2;
		if(func[i]!="")yuan[i]=x;
	}
	cout << "========================================\n"
	     << "|            MathX 2.0 Beta            |\n"
	     << "|      Enter [help] for commands       |\n"
	     << "|          Example: 2sin(pi/4)         |\n"
	     << "========================================\n";
}

int main() {
	ios::sync_with_stdio(0);
	cin.tie(0);
	cout.tie(0);
	Ready();
	while (1) Run();
}
