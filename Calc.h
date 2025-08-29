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
			if (int(sum) != sum) {
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
//			cout << "token: " << token << endl;
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
			int index = distance(begin(func), find(begin(func), end(func), token));
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
				else if (token == "isprime") res = high_precision_func::is_prime(u64(x));
				else if (token == "prime") {
					int cnt = 0;
					for (u64 i = 2;; i++)if (high_precision_func::is_prime(i))if (++cnt == x) {
								res = i;
								break;
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
					else if (token == "C") res = tgamma(y+1) / (tgamma(x+1) * tgamma(y - x+1));
					else if (token == "rand") {
						ll a = x, b = y;
						if(a>b)swap(a,b);
						random_device rd;
						mt19937 gen(rd());
						uniform_int_distribution<ll> dist(a,b);
						res = ld(dist(gen));
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

inline ld diff(const string& func, const ld x, const string& var_name){
	auto f = [&](ld y) { 
		variables[var_name] = to_string(y);
		return calc(func); 
	};
	ld h = 5e-7;
	return (f(x+h)-f(x-h))/(2*h);
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

inline string integral_calc(string& str, const string& var_name)
{
	str=ns(str);
	stack<string> st;
	queue<string> q;
	ll i=0;
	while (i<str.length())
	{
		string token;
		while (i<str.length() && str[i] != ' ') token += str[i++];
		
		if (token.empty()) {
			i++;
			continue;
		}
		if (token == var_name) st.push(token);
		else if (variables.find(token) != variables.end()) q.push(variables[token]);
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
			return "";
		}
		i++;
	}
	while (!st.empty())
	{
		q.push(st.top());
		st.pop();
	}
	auto nq=q;
	while (!nq.empty())
	{
		cout<<nq.front()<<" ";
		nq.pop();
	}
	cout<<"\n";
	return "";
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
	for (int i = 1; i < 9; ++i) x += p[i] / (z + cld(i, 0));
    
	cld t = z + g + 0.5L;
	return sqrt_2_pi * pow(t, z + 0.5L) * exp(-t) * x;
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
			int index = distance(begin(func), find(begin(func), end(func), token));
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
        return complex_diff(zx,t,can_var);
    };

    constexpr int max_iter = 20;
    cld integral = 0;
    cld prev_integral = 0;
    int n = 1;

    for (int iter = 0; iter < max_iter; ++iter) {
        constexpr ld tol = 1e-8;
        cld sum = f(z(minn)) * dz_dt(minn) + f(z(maxn)) * dz_dt(maxn);
        ld h = (maxn - minn) / n;

        for (int i = 1; i < n; ++i) {
            ld ti = minn + i * h;
            int weight = i % 2 == 1 ? 4 : 2;
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