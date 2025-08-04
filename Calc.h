#pragma once

#include "Initial.h"
#include "high_precision.h"

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
			preci = sum;
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

ld diff(string func, ld x, string var_name){
	auto f = [&](ld x) { 
		variables[var_name] = to_string(x);
		return calc(func); 
	};
	ld h = 1e-6;
	return (f(x+h)-f(x-h))/(2*h);
}