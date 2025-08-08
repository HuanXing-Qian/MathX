#pragma once
#include <map>
#include "Initial.h"

namespace high_precision_func {
	inline long long mul_mod(long long a, long long b, long long mod) {
	    long long res = 0;
	    while (b > 0) {
	        if (b & 1) res = (res + a) % mod;
	        a = a * 2 % mod;
	        b >>= 1;
	    }
	    return res;
	}
	
	inline long long pow_mod(long long a, long long e, long long mod) {
	    long long res = 1;
	    while (e > 0) {
	        if (e & 1) res = mul_mod(res, a, mod);
	        a = mul_mod(a, a, mod);
	        e >>= 1;
	    }
	    return res;
	}
	
	inline bool is_prime(long long n) {
	    if (n < 2) return false;
	    if (n == 2 || n == 3) return true;
	    if (n % 2 == 0) return false;
	
	    long long d = n - 1;
	    int s = 0;
	    while (d % 2 == 0) {
	        d /= 2;
	        s++;
	    }
	
	    const vector<long long> bases = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
	    
	    for (long long a : bases) {
	        if (a >= n) continue;
	        long long x = pow_mod(a, d, n);
	        if (x == 1 || x == n - 1) continue;
	        
	        bool composite = true;
	        for (int i = 0; i < s - 1; ++i) {
	            x = mul_mod(x, x, n);
	            if (x == n - 1) {
	                composite = false;
	                break;
	            }
	        }
	        if (composite) return false;
	    }
	    return true;
	}
	
	inline long long gcd(long long a, long long b) {
	    while (b != 0) {
	        long long t = b;
	        b = a % b;
	        a = t;
	    }
	    return a;
	}
	
	inline long long pollards_rho(long long n) {
	    if (n % 2 == 0) return 2;
	    if (n % 3 == 0) return 3;
	    if (n % 5 == 0) return 5;
	
	    random_device rd;
	    mt19937 gen(rd());
	    uniform_int_distribution<long long> dis(1, n - 1);
	
	    while (true) {
	        long long x = dis(gen);
	        long long y = x;
	        long long c = dis(gen);
	        long long d = 1;
	        
	        if (c == 0 || c == 2) c = 1;
	
	        auto f = [&](const long long de) { return (mul_mod(de, de, n) + c) % n; };
	
	        int power = 1, lam = 1;
	        while (d == 1) {
	            if (power == lam) {
	                x = y;
	                power *= 2;
	                lam = 0;
	            }
	            y = f(y);
	            lam++;
	            d = gcd(abs(x - y), n);
	        }
	        
	        if (d != n) return d;
	    }
	}
	
	inline map<long long, int> factorize(long long n) {
	    map<long long, int> factors;
	    
	    if (n == 1) return factors;
	    if (is_prime(n)) {
	        factors[n]++;
	        return factors;
	    }
	    
	    long long d = pollards_rho(n);
	    auto factors1 = factorize(d);
	    auto factors2 = factorize(n / d);
	    
	    for (auto& p : factors1) factors[p.first] += p.second;
	    for (auto& p : factors2) factors[p.first] += p.second;
	    
	    return factors;
	}
	
	inline string format_factorization(const map<long long, int>& factors) {
	    string result;
	    bool first = true;
	    
	    for (const auto& p : factors) {
	        if (!first) result += " * ";
	        first = false;
	        
	        if (p.second == 1) {
	            result += to_string(p.first);
	        } else {
	            result += to_string(p.first) + "^" + to_string(p.second);
	        }
	    }
	    
	    if (result.empty()) result = "1"; // Special case for n = 1
	    
	    return result;
	}

	inline bool is_prime(u64 n) {
		if (n < 2 || n % 6 % 4 != 1 && n % 6 % 4 != 5)
			return n == 2 || n == 3;
		u64 d = n - 1, s = 0;
		while (d % 2 == 0) d /= 2, s++;
		for (u64 a : {
		         2, 325, 9375, 28178, 450775, 9780504, 1795265022
		     }) {
			if (a >= n) break;
			u64 x = 1, p = d;
			for (; p > 0; p >>= 1, a = a * a % n)
				if (p & 1) x = x * a % n;
			if (x == 1 || x == n - 1) continue;
			for (u64 i = 0; i < s - 1 && x != n - 1; i++)
				x = x * x % n;
			if (x != n - 1) return false;
		}
		return true;
	}
	string add(string a, string b);
	string sub(string a, string b);

	inline string trim_zeros(const string &s) {
		auto pos = s.find_first_not_of('0');
		return (pos == string::npos) ? "0" : s.substr(pos);
	}

	inline string add(string a, string b) {
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

	inline string sub(string a, string b) {
		bool a_neg = (a[0] == '-');
		bool b_neg = (b[0] == '-');
		if (a_neg && b_neg) return sub(b.substr(1), a.substr(1));
		if (a_neg) return "-" + add(a.substr(1), b);
		if (b_neg) return add(a, b.substr(1));

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

	inline string div(string a, string b) {
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

	const double mm_pi = acos(-1.0);

	inline void parallel_fft(vector<complex<double>>& a, bool invert) {
		int n = a.size();
		if (n == 1) return;
		vector<complex<double>> a0(n / 2), a1(n / 2);
		for (int i = 0; i < n; i += 2) {
			a0[i / 2] = a[i];
			a1[i / 2] = a[i + 1];
		}

		parallel_fft(a0, invert);
		parallel_fft(a1, invert);

		double ang = 2 * mm_pi / n * (invert ? -1 : 1);
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

	inline string mul(string num1, string num2) {
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
			const int value = static_cast<int>(round(fa[i].real())) + carry;
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

	inline string pow(string base, string exponent) {
		string res = "1";
		while (exponent != "0") {
			if ((exponent.back() - '0') % 2 == 1) res = mul(res, base);
			base = mul(base, base);
			exponent = div(exponent, "2");
		}
		return res;
	}
}

inline string high_precision_calc(string str) {
	str = ns(str);
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
		if (isdigit(token[0]) || (token[0] == '-' && isdigit(token[1])))q.push(token);
		else if (token == "^") {
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
			if (lagg == 1)out("Expression error.\n", RED);
			else if (lagg == 2)out("伙计，表达式都写错了。\n", RED);
			else if (lagg == 3)out("Ошибка выражения.\n", RED);
			return 0;
		}
		i++;
	}

	while (!st.empty()) {
		q.push(st.top());
		st.pop();
	}

	if (rpn) {
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
			if (token == "+")result = high_precision_func::add(l, r);
			else if (token == "-") result = high_precision_func::sub(l, r);
			else if (token == "*") result = high_precision_func::mul(l, r);
			else if (token == "/") {
				if (r == "0") {
					out("inf\n", YELLOW);
					return 0;
				}
				result = high_precision_func::div(l, r);
			} else if (token == "^") result = high_precision_func::pow(l, r);
			else {
				if (lagg == 1)out("Expression error.\n", RED);
				else if (lagg == 2)out("伙计，表达式都写错了。\n", RED);
				else if (lagg == 3)out("Ошибка выражения.\n", RED);
				return 0;
			}
			stk.push(result);
		}
//		cout<<stk.top()<<" "<<token<<endl;
	}
	return stk.top();
}

