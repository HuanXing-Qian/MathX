#pragma once
#include "Initial.h"
#include "high_precision.h"
#include <cmath>
#include<stack>
#include<queue>

inline ld calc(string str) {
	str = ns(str);
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
				if (lagg == 1) out("The usage of this function is to add a non-negative integer after the function name.\n", RED);
				else if (lagg == 2) out("用法：请在函数名后输入一个非负整数。\n", RED);
				else if (lagg == 3) out("Использование: укажите целое неотрицательное число после имени функции.\n", RED);
				return 0;
			}
			ld sum = stold(token);
			if (ll(sum) != sum) {
				if (lagg == 1) out("The usage of this function is to add a non-negative integer after the function name.\n", RED);
				else if (lagg == 2) out("错误：必须输入整数。\n", RED);
				else if (lagg == 3) out("Ошибка: требуется целое число.\n", RED);
				return 0;
			}
			if (sum > 80) {
				if (lagg == 1) out("Error: precision exceeds long double limit.\n", RED);
				else if (lagg == 2) out("错误：精度值超过long double限制。\n", RED);
				else if (lagg == 3) out("Ошибка: точность превышает пределы long double.\n", RED);
				return 0;
			}
			precision = sum;
			return 0;
		}

		if (variables.find(token) != variables.end()) q.push(variables[token]);
		else if (isdigit(token[0]) || (token[0] == '-' && isdigit(token[1])))q.push(token);
		else if (prio(token) == 4 || prio(token) == 6) st.push(token);
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

	stack<ld> stk;
	while (!q.empty()) {
		string token = q.front();
		q.pop();
		if (isdigit(token[0]) || (token[0] == '-' && isdigit(token[1]))) stk.push(stold(token));
		else if (prio(token) < 4) {
			if (stk.size() < 2) {
				if (lagg == 1)out("Error: Insufficient " + token + " operators.\n", RED);
				else if (lagg == 2)out("没有足够的" + token + "操作数，即有左没右或者反之。\n", RED);
				else if (lagg == 3)out("Ошибка: не хватает операторов " + token + " для выполнения операции.\n", RED);
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
				if (!r) { // 

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
			ll index = distance(begin(func), find(begin(func), end(func), token));
			if (yuan[index] == 1) {
				ld x = stk.top();
				stk.pop();
				if (token == "sin") res = sin(remainder(x, 2*acos(-1)));
				else if (token == "cos") res = cos(remainder(x, 2*acos(-1)));
				else if (token == "tan") res = tan(remainder(x, 2*acos(-1)));
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
				else if (token == "atr") res = x / 180 * stold(constants.at("pi"));
				else if (token == "rta") res = x / stold(constants.at("pi")) * 180;
				else if (token == "isprime")
				{
					if (find(primes.begin(),primes.end(),ll(x))!=primes.end())res=1;
					else res = high_precision_func::is_prime(u64(x));
				} 
				else if (token == "prime") {
					if (x <= primes.size()) res = primes[x - 1];
					else {
						for (ll j = primes.back()==2?primes.back()+1:primes.back()+2;; j+=2) {
							if (high_precision_func::is_prime(j)) primes.push_back(j);
							if (x==primes.size()) {
								res = j;
								break;
							}
						}
					}
				}

			} else if (yuan[index] == 2) {
				ld x = stk.top();
				stk.pop();
				ld y = 0;
				if (!stk.empty())y = stk.top();
				if (!stk.empty()) {
					if (token == "min")res = min(x, y);
					else if (token == "max")res = max(x, y);
					else if (token == "arctan2")res = atan2(y, x);
					else if (token == "gcd") {
						if (ll(x) != x ||ll(y) != y) {
							out("Error: gcd requires integer inputs\n", RED);
							return 0;
						}
						res = high_precision_func::gcd(ll(x), ll(y));
					} else if (token == "lcm") {
						if (ll(x) != x || ll(y) != y) {
							out("Error: lcm requires integer inputs\n", RED);
							return 0;
						}
						res = ll(x*y) / high_precision_func::gcd(ll(x), ll(y));
					} else if (token == "hypot")res = sqrt(x * x + y * y);
					else if (token == "C")
					{
						res = tgamma(y+1) / (tgamma(x+1) * tgamma(y - x+1));
						Achievement[11].first.second++;
					}
					else if (token == "rand") {
						ll a = x, b = y;
						if(a>b)swap(a,b);
						random_device rd;
						mt19937 gen(rd());
						uniform_int_distribution<ll> dist(a,b);
						res = ld(dist(gen));
						Achievement[6].first.second++;
					}
					stk.pop();
				} else {
					out("Insufficient parameters.\n", RED);
					return 0;
				}
			}
//			if(!isfinite(res)||isnan(res)) return 0; 
			stk.push(res);
		} else if (prio(token) == 6) {
		    if (stk.empty()) {
		        out("No " + token + " function.\n", RED);
		        return 0;
		    }
		
		    auto it = deffunc.find(token);
		    if (it == deffunc.end()) {
		        out("Function " + token + " not found.\n", RED);
		        return 0;
		    }
		    const string& var_name = it->second.first; 
		    const string& expr = it->second.second; 
		    ld x = stk.top();
		    stk.pop();
		    string old_value;
		    bool existed = (variables.find(var_name) != variables.end());
		    if (existed) old_value = variables[var_name];
		    variables[var_name] = to_string(x);
		    ld res = calc(expr);
		    stk.push(res);
		    if (existed)variables[var_name] = old_value;
		    else variables.erase(var_name);
		}
	}
	return stk.top();
}

inline void big_pow(const string& str)
{
	istringstream iss(str);
	string cmd;
	string qw,qe;
	long double base, exponent;
	iss >> cmd >> qw >> qe;
	base=calc(qw), exponent = calc(qe);
	auto start = chrono::high_resolution_clock::now();
	long double logBase = log10l(base);
	long double logResult = exponent * logBase;
	long double B_float = floorl(logResult);
	ll B = static_cast<ll>(B_float);
	long double A = powl(10.0L, logResult - B_float);
	if (A >= 10.0L) {
		A /= 10.0L;
		B += 1;
	}

	if (A < 1.0L) {
		A *= 10.0L;
		B -= 1;
	}

	ll win = 1, a = static_cast<ll>(base), b = static_cast<ll>(exponent);
	ll mod = 1000000000000000000LL;
	a %= mod;

	auto mod_mul = [mod](ll x, ll y) -> ll {
		x %= mod;
		y %= mod;
		ll res = 0;
		while (y > 0) {
			if (y & 1)
				res = (res + x) % mod;
			x = (x * 2) % mod;
			y >>= 1;
		}
		return res;
	};

	while (b > 0) {
		if (b & 1)
			win = mod_mul(win, a);
		a = mod_mul(a, a);
		b >>= 1;
	}

	auto end = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<std::chrono::microseconds>(end - start);
	cout << fixed << setprecision(baoliu) << A << "e" << B << "\n";
	cout << "最后几位数字：" << win << "\n";
	if (timing) cout << "Use:" << duration.count() / 1000 << " ms." << "\n";
}

inline cld operator+(const cld& cldh, const ld rhs)
{
	return cldh+cld(rhs);	
}

inline cld complex_gamma(cld z) {
	if (z.imag()==0) return tgamma(z.real());
	
	constexpr ld g = 7.0L;
	constexpr ld sqrt_2_pi = 2.506628274631000502415765284811L;
    
	static constexpr ld p[] = {
		0.99999999999980993L,
		676.5203681218851L,
		-1259.1392167224028L,
		771.32342877765313L,
		-176.61502916214059L,
		12.507343278686905L,
		-0.13857109526572012L,
		9.9843695780195716e-6L,
		1.5056327351493116e-7L
	};
    
	if (z.real() <= 0.5L) return PI / (sin(PI * z) * complex_gamma(cld(1.0L) - z));

	z -= 1.0L;
	cld x = p[0];
	for (ll i = 1; i < 9; ++i) x += p[i] / (z + cld(i, 0));
    
	cld t = z + g + 0.5L;
	return sqrt_2_pi * pow(t, z + 0.5L) * exp(-t) * x;
}

inline cld diff_calc(string str,const cld z,const string&var_name)
{
	complex_variables[var_name]=z;
	str=ns(str);
	stack<string> st;
	queue<string> q;
	ll i=0;
	while (i<str.length())
	{
		string token;
		while (i<str.length() && str[i] != ' ') token += str[i++];
		if (token.empty())
		{
			i++;
			continue;
		}
		if (complex_variables.find(token) != complex_variables.end())q.push(token);
		else if (isdigit(token[0])|| (token[0] == '-' && isdigit(token[1]))||token=="i")q.push(token);
		else if (prio(token) == 4||prio(token)==6) st.push(token);
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
		}
		else {
			if (lagg == 1)out("Expression error.\n", RED);
			else if (lagg == 2)out("伙计，表达式都写错了。\n", RED);
			else if (lagg == 3)out("Ошибка выражения.\n", RED);
			return 0;
		}
		i++;
	}
	while (!st.empty())
	{
		q.push(st.top());
		st.pop();
	}
	if (rpn)
	{
		cout<<"RPN: ";
		auto nq = q;
		while (!nq.empty())
		{
			cout<<nq.front()<<" ";
			nq.pop();
		}
		cout<<"\n";
	}
	stack<pair<cld,cld>> stk;
	while (!q.empty())
	{
		string token = q.front();
		q.pop();
		if (complex_variables.find(token) != complex_variables.end())stk.push({complex_variables[token],token!=var_name?0:1});
		else if (isdigit(token[0])|| (token[0] == '-' && isdigit(token[1])) || token=="i")
		{
			cld n;
			if (token=="i")n.imag(1);
			else n.real(stold(token));
			stk.push({n,0});
		}else if (prio(token)<4)
		{
			if (stk.size() < 2) {
				if (lagg == 1)out("Error: Insufficient " + token + " operators.\n", RED);
				else if (lagg == 2)out("没有足够的" + token + "操作数，即有左没右或者反之。\n", RED);
				else if (lagg == 3)out("Ошибка: не хватает операторов " + token + " для выполнения операции.\n", RED);
				return 0;
			}
			auto r=stk.top();
			stk.pop();
			auto l=stk.top();
			stk.pop();
			pair<cld,cld> result;
			if (token == "+") result = {l.first+r.first,l.second+r.second};
			else if (token == "-") result = {l.first-r.first,l.second-r.second};
			else if (token == "*") result = {l.first*r.first,l.second*r.first+l.first*r.second};
			else if (token == "/") {
				if (r.first==cld(0))
				{
					out("inf\n",YELLOW);
					return 0;
				}
				result = {l.first/r.first, (r.first*l.second-l.first*r.second)/r.first/r.first};
			} else if (token == "^") result = {pow(l.first,r.first),pow(l.first,r.first)*(r.second*log(l.first)+r.first*l.second/l.first)};
			stk.push(result);
		}else if (prio(token)==4)
		{
			if (stk.empty()) {
				out("No " + token + " function.\n", RED);
				return 0;
			}
			pair<cld,cld> res;
			ll index = distance(begin(func), find(begin(func), end(func), token));
			if (yuan[index]==1)
			{
				auto x = stk.top();
				stk.pop();
				cld z=x.first;
				cld dz = x.second;
				if (token=="sin")res = {sin(z),cos(z)*dz};
				else if (token=="cos")res = {cos(z),-sin(z)*dz};
				else if (token=="tan")res = {tan(z),(cld(1)+tan(z)*tan(z))*dz};
				else if (token=="arcsin")res = {asin(z),dz/sqrt(cld(1)-z*z)};
				else if (token=="arccos")res = {acos(z), -dz / sqrt(cld(1) - z * z)};
				else if (token=="arctan")res = {atan(z),dz / (cld(1) + z * z)};
				else if (token=="sinh")res = {sinh(z),cosh(z)*dz};
				else if (token=="cosh")res = {cosh(z),sinh(z)*dz};
				else if (token=="tanh")res = {tanh(z),(cld(1) - tanh(z) * tanh(z)) * dz};
				else if (token=="sqrt")res = {sqrt(z),dz/(cld(2)*sqrt(z))};
				else if (token=="exp")res = {exp(z),exp(z)*dz};
				else if (token=="log"||token=="ln")res = {log(z),dz/z};
				else if (token=="log10")res = {log10(z),dz/z/cld(log(10))};
				else if (token == "abs") {
					cld abs_u = abs(z);
					if (abs_u.imag()==0&&abs_u.real()==0) {
						out("Error: Derivative of abs(0) is undefined\n", RED);
						return 0;
					}
					cld deriv = (z * conj(dz)).real() / abs_u;
					res = { abs_u, deriv };
				}
				else if (token=="sinc")res = {z.imag()==0&&z.real()==0 ? cld(1) : sin(z)/z,dz*(z*cos(z)-sin(z)/z/z)};
				else if (token=="Gamma")
				{
					const ld h = 5e-7;
					res = {complex_gamma(z),dz*(complex_gamma(z+h)-complex_gamma(z-h))/cld(2)};
				}
			}
			stk.push(res);
		}
	}
	return stk.top().second;
}

inline cld complex_calc(string str)
{
	str=ns(str);
	stack<string> st;
	queue<string> q;
	ll i=0;
	while (i<str.length())
	{
		string token;
		while (i<str.length() && str[i] != ' ') token += str[i++];
		if (token.empty())
		{
			i++;
			continue;
		}
		if (complex_variables.find(token) != complex_variables.end())q.push(token);
		else if (isdigit(token[0])|| (token[0] == '-' && isdigit(token[1]))||token=="i")q.push(token);
		else if (prio(token) == 4||prio(token)==6) st.push(token);
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
		}
		else {
			if (lagg == 1)out("Expression error.\n", RED);
			else if (lagg == 2)out("伙计，表达式都写错了。\n", RED);
			else if (lagg == 3)out("Ошибка выражения.\n", RED);
			return 0;
		}
		i++;
	}
	while (!st.empty())
	{
		q.push(st.top());
		st.pop();
	}
	if (rpn)
	{
		cout<<"RPN: ";
		auto nq = q;
		while (!nq.empty())
		{
			cout<<nq.front()<<" ";
			nq.pop();
		}
		cout<<"\n";
	}
	stack<cld> stk;
	while (!q.empty())
	{
		string token = q.front();
		q.pop();
		if (complex_variables.find(token) != complex_variables.end())stk.push(complex_variables[token]);
		else if (isdigit(token[0])|| (token[0] == '-' && isdigit(token[1])) || token=="i")
		{
			cld n;
			if (token=="i")n.imag(1);
			else n.real(stold(token));
			stk.push(n);
		}
		else if (prio(token)<4)
		{
			if (stk.size() < 2) {
				if (lagg == 1)out("Error: Insufficient " + token + " operators.\n", RED);
				else if (lagg == 2)out("没有足够的" + token + "操作数，即有左没右或者反之。\n", RED);
				else if (lagg == 3)out("Ошибка: не хватает операторов " + token + " для выполнения операции.\n", RED);
				return 0;
			}
			cld r=stk.top();
			stk.pop();
			cld l=stk.top();
			stk.pop();
			cld result;
			if (token == "+") result = l + r;
			else if (token == "-") result = l - r;
			else if (token == "*") result = l * r;
			else if (token == "/") {
				if (r.imag() == 0&&r.real()==0) {
					out("inf\n", YELLOW);
					return 0;
				}
				result = l / r;
			} else if (token == "^") result = pow(l, r);
			stk.push(result);
		}
		else if (prio(token)==4)
		{
			if (stk.empty()) {
				out("No " + token + " function.\n", RED);
				return 0;
			}
			cld res;
			ll index = distance(begin(func), find(begin(func), end(func), token));
			if (yuan[index]==1)
			{
				cld z = stk.top();
				stk.pop();
				if (token=="sin")res = sin(z);
				else if (token=="cos")res = cos(z);
				else if (token=="tan")res = tan(z);
				else if (token=="arcsin")res = asin(z);
				else if (token=="arccos")res = acos(z);
				else if (token=="arctan")res = atan(z);
				else if (token=="sinh")res = sinh(z);
				else if (token=="cosh")res = cosh(z);
				else if (token=="tanh")res = tanh(z);
				else if (token=="sqrt")res = sqrt(z);
				else if (token=="exp")res = exp(z);
				else if (token=="log"||token=="ln")res = log(z);
				else if (token=="log10")res = log10(z);
				else if (token=="abs")res = abs(z);
				else if (token=="sinc")res = (z.imag()==0&&z.real()==0) ? cld(1) : sin(z)/z;
				else if (token=="Gamma")res = complex_gamma(z);
			}
			stk.push(res);
		}else if (prio(token) == 6) {
			if (stk.empty()) {
				out("No " + token + " function.\n", RED);
				return 0;
			}
		
			auto it = deffunc.find(token);
			if (it == deffunc.end()) {
				out("Function " + token + " not found.\n", RED);
				return 0;
			}
			const string& var_name = it->second.first; 
			const string& expr = it->second.second; 
			cld x = stk.top();
			stk.pop();
			cld old_value;
			bool existed = complex_variables.find(var_name) != complex_variables.end();
			if (existed) old_value = complex_variables[var_name];
			complex_variables[var_name] = x;
			cld res = complex_calc(expr);
			stk.push(res);
			if (existed)complex_variables[var_name] = old_value;
			else complex_variables.erase(var_name);
		}
	}
	return stk.top();
}

inline string diff_calc_fu(string str,string var_name)
{
	str=ns(str);
	stack<string> st;
	queue<string> q;
	ll i=0;
	while (i<str.length())
	{
		string token;
		while (i<str.length() && str[i] != ' ') token += str[i++];
		if (token.empty())
		{
			i++;
			continue;
		}
		if (complex_variables.find(token) != complex_variables.end()||token==var_name)q.push(token);
		else if (isdigit(token[0])|| (token[0] == '-' && isdigit(token[1]))||token=="i")q.push(token);
		else if (prio(token) == 4||prio(token)==6) st.push(token);
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
		}
		else {
			if (lagg == 1)out("Expression error.\n", RED);
			else if (lagg == 2)out("伙计，表达式都写错了。\n", RED);
			else if (lagg == 3)out("Ошибка выражения.\n", RED);
			return 0;
		}
		i++;
	}
	while (!st.empty())
	{
		q.push(st.top());
		st.pop();
	}
	if (rpn)
	{
		cout<<"RPN: ";
		auto nq = q;
		while (!nq.empty())
		{
			cout<<nq.front()<<" ";
			nq.pop();
		}
		cout<<"\n";
	}
	stack<pair<string,string>> stk;
	while (!q.empty())
	{
		string token = q.front();
		q.pop();
		if (complex_variables.find(token) != complex_variables.end()&&token!=var_name)stk.push({to_string(complex_variables[token].real())+"+"+to_string(complex_variables[token].imag())+"i","0"});
		else if (token==var_name)stk.push({var_name,"1"});
		else if (isdigit(token[0])|| (token[0] == '-' && isdigit(token[1])) || token=="i")
		{
			cld n;
			if (token=="i")n.imag(1);
			else n.real(stold(token));
			if (n.imag()==0)stk.push({to_string(n.real()),"0"});
			else stk.push({to_string(n.real())+"+"+to_string(n.imag())+"i","0"});
		}
		else if (prio(token)<4)
		{
			if (stk.size() < 2) {
				if (lagg == 1)out("Error: Insufficient " + token + " operators.\n", RED);
				else if (lagg == 2)out("没有足够的" + token + "操作数，即有左没右或者反之。\n", RED);
				else if (lagg == 3)out("Ошибка: не хватает операторов " + token + " для выполнения операции.\n", RED);
				return 0;
			}
			auto r=stk.top();
			stk.pop();
			auto l=stk.top();
			stk.pop();
			pair<string,string> result;
			if (token == "+") result = {l.first+"+"+r.first,l.second+"+"+r.second};
			else if (token == "-") result = {l.first+"-"+r.first,l.second+"-"+r.second};
			else if (token == "*")
			{
				string hou;

				// 处理 l.second
				if (l.second != "1") {
					hou = l.second;
					if (r.first != "1" && r.first != "-1") {
						hou += "*" + r.first;  // 只有当 r.first 不是 1 或 -1 时才添加
					}
				} else if (r.first != "1" && r.first != "-1") {
					hou = r.first;
				}

				hou += "+";

				// 处理 l.first
				if (l.first != "1") {
					if (r.second != "1" && r.second != "-1") {
						hou += l.first + "*" + r.second;  // 只有当 r.second 不是 1 或 -1 时才添加
					} else {
						hou += l.first;
					}
				} else if (r.second != "1" && r.second != "-1") {
					hou += r.second;
				}

				// 构建结果
				string first_part = (l.first != "1" ? l.first : "") + (r.first != "1" ? "*" + r.first : "");
				string second_part = hou;

				result = {first_part, second_part};
			}
			else if (token == "/") {
			    if (r.first == "0")
			    {
			        out("inf\n", YELLOW);
			        return 0;
			    }
			    string hou = "((" + r.first + ")*(" + l.second + ") - (" + l.first + ")*(" + r.second + ")) / (" + r.first + ")^2";
			    
			    // 简化 hou
			    if (r.first == "1") {
			        hou = "(" + l.second + " - (" + l.first + ")*(" + r.second + "))";
			    }
			    if (l.second == "0") {
			        hou = "-" + l.first + "*" + r.second + " / " + r.first + "^2";
			    }
			    if (r.second == "0") {
			        hou = l.second + " / " + r.first;
			    }
			    if (l.second == "0" && r.second == "0") {
			        hou = "0";
			    }

			    result = { l.first + "/" + r.first, hou };
			}
			else if (token == "^")
			{
				string first_;
				string second_;
				for (auto v:l.first)if (!isdigit(v)&&v!='.')
				{
					first_ = "("+l.first+")";
					break;
				}
				if (first_.empty())first_ = l.first;
				bool f=1;
				for (auto v:r.first)if (!isdigit(v)&&v!='.')
				{
					first_ += "^("+r.first+")";
					f=0;
					break;
				}
				if (f)
				{
					if (r.first=="1")first_=l.first;
					else if (r.first=="0")first_="1";
					else first_ += r.first;
				}
				bool type1=0,type2=0;
				
				for (auto v:l.first)if (!isdigit(v)&&v!='.'&&v!='i'&&v!='+'&&v!='('&&v!=')')
				{
					type1=1;
					break;
				}
				for (auto v:r.first)if (!isdigit(v)&&v!='.'&&v!='i'&&v!='+'&&v!='('&&v!=')')
				{
					type2=1;
					break;
				}
				// cout << type1 << " " << type2 << "\n";
				if (type1&&!type2) // x^a r.first是常数
				{
					cld cc = complex_calc(r.first);
					if (l.first.size()!=1)
					{
						if (cc.imag()==0)second_ = r.first+"("+l.first+")";
						else second_ = "("+to_string(cc.real())+"+"+to_string(cc.imag())+"i)("+l.first+")";
					}
					else
					{
						if (cc.imag()==0)second_ = r.first+l.first;
						else second_ = "("+to_string(cc.real())+"+"+to_string(cc.imag())+"i)"+l.first;
					}
					if (cc!=cld(2))
					{
						if (cc.imag()==0) second_ += "^"+to_string(cc.real()-1); // 实数
						else second_ += "^("+to_string(cc.real()-1)+"+"+to_string(cc.imag())+"i)";
					}
				}
				else if (!type1&&type2) // a^x, l.first^r.first l.first^r.first*ln(l.first)*r.second
				{
					cld cc = complex_calc(l.first);
					if (cc.imag()==0)
					{
						second_ = to_string(cc.real())+"^";
						if (r.first.size()==1)second_+=r.first;
						else second_+="("+r.first+")";
						second_+="*ln("+to_string(cc.real())+")";
						if (r.second!="1")second_+="("+r.second+")";
					}else
					{
						second_ = "("+to_string(cc.real())+"+"+to_string(cc.imag())+"i)^";
						if (r.first.size()==1)second_+=r.first;
						else second_+="("+r.first+")";
						second_+="*ln("+to_string(cc.real())+"+"+to_string(cc.imag())+"i)";
						if (r.second!="1")second_+="("+r.second+")";
					}
					if (cc==cld(1))second_="0";
				}
				else // f(x)^g(x)
				{
					if (l.first.size()==1)second_+=l.first+"^";
					else second_+="("+l.first+")^";
					if (r.first.size()==1)second_+=r.first+"*(";
					else second_+="("+r.first+")*(";
					if (l.first!="1"){
						if (r.second.find_first_not_of("0123456789") != string::npos)second_+="("+r.second+")";
						else if (r.second!="1")second_+=r.second;
						if (r.second!="0")second_+="ln("+l.first+")+";
					}
					if (r.first==l.first)second_+="1";
					else
					{
						if (r.first.size()!=1)second_+="("+r.first+")";
						else second_+=r.first;
						if (l.second.size()!=1)second_+="("+l.second+")";
						else if (l.second!="1")second_+="*"+l.second;
						second_+="/";
						if (l.first.size()!=1)second_+="("+l.first+")";
						else second_+=l.first;
					}
					
					second_+=")";
				}
				result = {first_,second_};
			}
			stk.push(result);
		}
		else if (prio(token)==4)
		{
			if (stk.empty()) {
				out("No " + token + " function.\n", RED);
				return 0;
			}
			pair<string,string> res;
			auto x = stk.top();
			stk.pop();
			string z=x.first;
			string dz = x.second;
			if (token=="sin")
				{
					string hou;
					if (dz!="0")
					{
						hou="cos("+z+")";
						if (dz!="1")hou+="*("+dz+")";
					}
					else hou="0";
					
					res = {"sin("+z+")",hou};
				}
			else if (token=="cos")
				{
					string hou;
					if (dz!="0")
					{
						hou="(-sin("+z+"))";
						if (dz!="1")hou+="*("+dz+")";
					}
					else hou="0";
					res = {"cos("+z+")",hou};
				}
			else if (token=="tan") {
				    string hou;
				    if (dz != "0") {
				    	for (auto v:dz)if (!isdigit(v))
				    	{
				    		hou = "("+dz+")";
				    		break;
				    	}
				    	if (hou.empty())hou=dz;
				    	hou+="/(cos("+z+")^2)";
				    } else hou = "0";
				    res = {"tan(" + z + ")", hou};
				} else if (token=="arcsin") {
					string hou;
					if (dz != "0") {
						for (auto v:dz)if (!isdigit(v))
						{
							hou = "("+dz+")";
							break;
						}
						if (hou.empty())hou=dz;
						hou+="/(2*sqrt("+z+"))";
					} else hou = "0";
					res = {"arcsin(" + z + ")", hou};
				} else if (token=="arccos") {
				    string hou;
				    if (dz != "0") {
				        hou = "-" + dz + "/sqrt(1-(" + z + ")^2)";
				        if (dz == "-1" || dz == "1") hou = dz + "/sqrt(1-(" + z + ")^2)";
				    } else hou = "0";
				    res = {"arccos(" + z + ")", hou};
				} else if (token=="arctan") {
				    string hou;
				    if (dz != "0") {
				        hou = "(" + dz + ")/(1+(" + z + ")^2)";
				        if (dz == "-1" || dz == "1") hou = dz + "/(1+(" + z + ")^2)";
				    } else hou = "0";
				    res = {"arctan(" + z + ")", hou};
				} else if (token=="sinh") {
				    string hou;
				    if (dz != "0") {
				        hou = "cosh(" + z + ")*(" + dz + ")";
				        if (dz == "1") hou = "cosh(" + z + ")";
				        else if (dz == "-1") hou = "-cosh(" + z + ")";
				    } else hou = "0";
				    res = {"sinh(" + z + ")", hou};
				} else if (token=="cosh") {
				    string hou;
				    if (dz != "0") {
				        hou = "sinh(" + z + ")*(" + dz + ")";
				        if (dz == "1") hou = "sinh(" + z + ")";
				        else if (dz == "-1") hou = "-sinh(" + z + ")";
				    } else hou = "0";
				    res = {"cosh(" + z + ")", hou};
				} else if (token=="tanh") {
				    string hou;
				    if (dz != "0") {
				        hou = "(1-tanh(" + z + ")^2)*(" + dz + ")";
				        if (dz == "1") hou = "1-tanh(" + z + ")^2";
				        else if (dz == "-1") hou = "tanh(" + z + ")^2-1";
				    } else hou = "0";
				    res = {"tanh(" + z + ")", hou};
				} else if (token=="sqrt") {
				    string hou;
				    if (dz != "0") {
				        hou = "(" + dz + ")/(2*(" + z + "))";
				        if (dz == "1") hou = "1/(2*sqrt(" + z + "))";
				        else if (dz == "-1") hou = "-1/(2*sqrt(" + z + "))";
				    } else hou = "0";
				    res = {"sqrt(" + z + ")", hou};
				} else if (token=="exp") {
				    string hou;
				    if (dz != "0") {
				        hou = "exp(" + z + ")*(" + dz + ")";
				        if (dz == "1") hou = "exp(" + z + ")";
				        else if (dz == "-1") hou = "-exp(" + z + ")";
				    } else hou = "0";
				    res = {"exp(" + z + ")", hou};
				} else if (token=="log" || token=="ln") {
				    string hou;
				    if (dz != "0") {
				        hou = "(" + dz + ")/(" + z + ")";
				        if (dz == "1") hou = "1/(" + z + ")";
				        else if (dz == "-1") hou = "-1/(" + z + ")";
				    } else hou = "0";
				    res = {"ln(" + z + ")", hou};
				} else if (token=="log10") {
				    string hou;
				    if (dz != "0") {
				        hou = "(" + dz + ")/((" + z + ")*ln(10))";
				        if (dz == "1") hou = "1/((" + z + ")*ln(10))";
				        else if (dz == "-1") hou = "-1/((" + z + ")*ln(10))";
				    } else hou = "0";
				    res = {"log10(" + z + ")", hou};
				}else if (token == "abs") {
					// 值计算: |u|
					string val_expr = "abs(" + z + ")";
    
					// 导数计算: (u * u') / |u|
					string deriv_expr = 
						"((" + z + ") * (" + dz + ")) / abs(" + z + ")";
    
					// 检查 u=0 的情况（可选，符号计算可能无法静态检测）
					res = {val_expr, deriv_expr};
				}else if (token == "sinc") {
					string val_expr = "sinc("+z+")";
					string deriv_expr = 
						"(((" + z + ")*cos(" + z + ") - sin(" + z + ")) / (" + 
						z + ")^2) * (" + dz + ")";
    
					res = {val_expr, deriv_expr};
				}
			
			stk.push(res);
		}
		else if (prio(token) == 6) {
			if (stk.empty()) {
				out("No " + token + " function.\n", RED);
				return 0;
			}
		
			auto it = deffunc.find(token);
			if (it == deffunc.end()) {
				out("Function " + token + " not found.\n", RED);
				return "0";
			}
			const string& basic_string = it->second.first; 
			const string& expr = it->second.second; 
			stk.emplace(expr,diff_calc_fu(expr,basic_string));
		}
	}
	return stk.top().second;
}

inline cld complex_diff(const string& func, const cld z, const string& var_name) {
	auto f = [&](cld y) {
		complex_variables[var_name] = y;
		return complex_calc(func);
	};

	const ld h = 5e-7;
	return (f(z+h)-f(z-h))/(2*h);
}

inline cld complex_integral(const string& func, const string& zx, const string& can_var, const ld minn, const ld maxn, const string& var_name) {
    auto f = [&](cld y) { 
        complex_variables[var_name] = y;
        return complex_calc(func);
    };
    auto z = [&](ld t) -> cld {
        // Save current parameter value if it exists
        cld old_param_val;
        bool param_existed = complex_variables.find(can_var) != complex_variables.end();
        if (param_existed) old_param_val = complex_variables[can_var];
        
        // Set the parameter value
        complex_variables[can_var] = t;
        
        // Evaluate the path
        cld path_val;
        try {
            path_val = complex_calc(zx);
        } catch (...) {
            if (param_existed) complex_variables[can_var] = old_param_val;
            else complex_variables.erase(can_var);
            throw;
        }
        
        // Restore original parameter value
        if (param_existed) complex_variables[can_var] = old_param_val;
        else complex_variables.erase(can_var);
        
        return path_val;
    };

    // Numerical derivative for path calculation
    auto dz_dt = [&](ld t) -> cld {
    	complex_variables[can_var]=t;
        return complex_calc(diff_calc_fu(zx,can_var));
    };

    constexpr ll max_iter = 20;
    cld integral = 0;
    cld prev_integral = 0;
    ll n = 1;

    for (ll iter = 0; iter < max_iter; ++iter) {
        constexpr ld tol = 1e-8;
        cld sum = f(z(minn)) * dz_dt(minn) + f(z(maxn)) * dz_dt(maxn);
        ld h = (maxn - minn) / n;

        for (ll i = 1; i < n; ++i) {
            ld ti = minn + i * h;
            ll weight = i % 2 == 1 ? 4 : 2;
            sum += cld(weight) * f(z(ti)) * dz_dt(ti);
        }

        integral = sum * h / 3.0L;
        if (iter > 0 && abs(integral - prev_integral) < tol) {
            break;
        }
        prev_integral = integral;
        n *= 2;
    }

    return integral;
}