// Developers: Shi Yuze
#pragma GCC optimize("O3,unroll-loops")
#define _CRT_SECURE_NO_WARNINGS
#include "Calc.h"
vector<node>history;

static bool chao(string str) {
	if (str=="easter egg")
	{
		random_device rd;
		mt19937 gen(rd());
		uniform_int_distribution<ll> dist(0,little_knowledge.size()-1);
		int index = dist(gen);
		cout<<little_knowledge.at(index)<<"\n";
		int In = knowledge_vis.size();
		knowledge_vis.insert(index);
		if (In<knowledge_vis.size())Achievement[1].first.second++;
		return 0;
	}
	
	if (str=="list functions")
	{
		int x=1;
		for (auto v:deffunc)cout << x++ << ":" << v.first << "(" << v.second.first << ") = " << v.second.second << endl;
		return 0;
	}
	
	if (str=="time")
	{
		time_t rawtime;
		tm * timeinfo;

		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		if (lagg==1)printf ( "The current date/time is: %s", asctime (timeinfo) );
		else if (lagg==2)printf ( "当前的时间是: %s", asctime (timeinfo) );
		else if (lagg==3)printf ( "Текущее время: %s", asctime (timeinfo) );
		return false;
	}
	
	if (str=="complex"||str=="complex on")
	{
		complex_mode = true;
		if (lagg==1)out("Complex mode on.\n",GREEN);
		else if (lagg==2)out("复数模式已开启。\n",GREEN);
		else if (lagg==3)out("Режим множественного числа включен.\n",GREEN);
		return false;
	}

	if (str=="complex off")
	{
		complex_mode = 0;
		if (lagg==1)out("Complex mode off.\n","\033[90m");
		else if (lagg==2)out("复数模式已关闭。\n","\033[90m");
		else if (lagg==3)out("Система множественных чисел отключена.\n","\033[90m");
		return false;
	}
	
	if (str.rfind("study", 0) == 0) {
		string func_name = str.substr(5);
		if (func_name.empty()) {
			int x = 1;
			for (auto v : ledge) out(to_string(x++) + ": " + v.first + "\n", CYAN);
			return 0;
		}
		func_name.erase(func_name.begin());
		string best_match;
		int max_score = 0;

		for (auto &v : ledge) {
			for (auto &u : v.second) {
				const string &target = u.first;
				int score = 0;
				bool full_match = true;
				if (!func_name.empty() && !target.empty() && func_name[0] == target[0])
					score += 2;
				for (int i = 0; i < max(func_name.size(), target.size()); i++) {
					if (i < func_name.size() && i < target.size() &&
						func_name[i] == target[i]) {
						score += 1;
						} else if (i < min(func_name.size(), target.size())) {
							full_match = false;
						}
				}

				if (full_match && u.first.size() == func_name.size()) {
					for (auto t : u.second) out(t + "\n", CYAN);
					return 0;
				}
				if (score > max_score) {
					max_score = score;
					best_match = target;
				}
			}
		}

		if (max_score > func_name.size() / 2) {
			if (lagg == 1)out("Did you mean \"" + best_match + "\"?\n", YELLOW);
			else if (lagg == 2)out("您想说的是\"" + best_match + "\"？\n", YELLOW);
			else if (lagg == 3)out("Вы хотели сказать\"" + best_match + "\"?\n", YELLOW);
		} else {
			if (lagg == 1)out("No matching concept found.\n", RED);
			else if (lagg == 2)out("未找到匹配项。\n", RED);
			else if (lagg == 2)out("Совпадений не найдено.\n", RED);
		}
		return 0;
	}
	if (str.rfind("search ", 0) == 0) {
		string func_name = str.substr(7);
		if (ledge.find(func_name) != ledge.end()) {
			int x = 1;
			for (auto v : ledge[func_name]) out(to_string(x++) + ": " + v.first + "\n", CYAN);
		} else {
			if (lagg == 1)out("No matching concept found.\n", RED);
			else if (lagg == 2)out("未找到匹配项。\n", RED);
			else if (lagg == 2)out("Совпадений не найдено.\n", RED);
		}
		return 0;
	}
	if (str.empty()) return 0;
	if (str == "exit"||str=="quit") exit(0);
	if (str == "clear history") {
		if (lagg == 1)out("Finished.\n", GREEN);
		else if (lagg == 2)out("完毕。唉，但愿我不要再说类似“已完成”这样的话了。\n", GREEN);
		else if (lagg == 3)out("Завершено.\n", GREEN);
		history.clear();
		return 0;
	}

	if (str.rfind("factor") == 0) {
		str.erase(0, 7);
		if (!str.empty()) {
			ll n = ll(calc(str));
			if (n <= 1) {
				cout << "请输入大于1的正整数" << "\n";
				return 0;
			}
			auto factors = high_precision_func::factorize(n);
			string factorization = high_precision_func::format_factorization(factors);
			cout<<"质因数分解："<< factorization << "\n";
			Achievement[5].first.second++;
		}
		return 0;
	}
	
	if (str.rfind("def")==0){
	    str.erase(0,4);
	    if(!str.empty()){
	    	istringstream iss(str);
	    	string var_name, func_name, func;
	    	iss >> func_name >> var_name;
	    	if (func_name.empty() || var_name.empty()) {
	    		out("Invalid function definition syntax.\n", RED);
	    		return 0;
	    	}
	    	getline(iss, func);
	    	func.erase(0, func.find_first_not_of(" \t"));
	    	func.erase(func.find_last_not_of(" \t") + 1);

	    	string old_value;
	    	bool var_existed = variables.find(var_name) != variables.end();
	    	if (var_existed) old_value = variables[var_name];

	    	variables[var_name] = "0"; 
	        bool f=true;
	    	str = ns(func);
	    	stack<string> st;
	    	queue<string> q;
	    	ll i = 0;
	    	while (i < str.size())
	    	{
	    		string token;
	    		while (i < str.size() && str[i] != ' ') {
	    			token += str[i];
	    			i++;
	    		}
		
	    		if (token.empty()) {
	    			i++;
	    			continue;
	    		}
				
	    		if (variables.find(token) != variables.end() || token == var_name) q.push(variables[token]);
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
	    		} else f=0;
	    	}
	    	if (!f)
	    	{
	    		if (var_existed) variables[var_name] = old_value;
	    		else variables.erase(var_name);
	    		return 0;
	    	}
	    	deffunc[func_name] = {var_name, func};
	    	if (var_existed) variables[var_name] = old_value;
	    	else variables.erase(var_name);

	    	cout << func_name << "(" << var_name << ") = " << func << endl;
	    }
	    return 0;
	}

	/*if (str.rfind("int")==0)
	{
		str.erase(0,4);
		if (!str.empty())
		{
			istringstream iss(str);
			string var_name, func_name, func;
			iss >> func_name >> var_name;
			if (func_name.empty() || var_name.empty())
			{
				out("你输入错了。\n", RED);
				return 0;
			}
			if (var_name.size()>1)
			{
				out("变量名必须是一个字母！\n",RED);
				return 0;
			}
			getline(iss, func);
			func.erase(0, func.find_first_not_of(" \t"));
			func.erase(func.find_last_not_of(" \t") + 1);
			deffunc[func_name] = {var_name, integral_calc(func,var_name)};
			cout << "Integral: " << deffunc[func_name].second << "\n";
			return 0;
		}
	}*/
	
	if (str.rfind("int")==0){
		int n = 5000;
		str.erase(0,4);
		if(!str.empty()){
			istringstream iss(str);
			string var_name, func;
			ld minn, maxn;
			string a,b;
			iss >> var_name >> a >> b;
			minn = calc(a);
			maxn = calc(b);
			getline(iss, func);
			string bl;
			if(!isalpha(var_name[0])){
				if (lagg == 1)out("The variable name cannot be used in Numbers at first.\n", YELLOW);
				else if (lagg == 2)out("变量名开头不能用数字。。\n", YELLOW);
				else if (lagg == 3)out("Переменные не могут начинаться с чисел.\n", YELLOW);
				return 0;
			} 
			bool bao=0;
			if(variables.find(var_name)!=variables.end())bl = variables[var_name], bao=1;
			auto start = chrono::high_resolution_clock::now();
			auto f = [&](ld x) { 
			    variables[var_name] = to_string(x);
			    return calc(func); 
			};
			int huo=1;
			if(maxn<minn)swap(maxn,minn),huo=-1;
			ld h = (maxn - minn)/n;
			ld sum = 0;
			sum += f(maxn)+f(minn);
			rep(i,1,n-1,1){
				ld x = minn+i*h;
				ld fx = f(x);
				if (!isfinite(fx) || isnan(x)) fx = 0;
				if (f(x)>0&&f(x+h)<0||f(x)<0&&f(x+h)>0)fx=0;
			    sum += i&1 ? 4*fx : 2*fx;
			}
			cout << "Result: " << huo*sum*h/3 << "\n";
			variables["ans"] = to_string(sum*h/3);
			complex_variables["ans"]=cld(sum*h/3);
			auto end = chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<std::chrono::microseconds>(end - start);
			if (timing) cout << "Use:" << duration.count() / 1000 << " ms." << "\n";
			if(bao)variables[var_name]=bl;
			else variables.erase(var_name);
			Achievement[3].first.second++;
		}
		return 0;
	}

	if (str.rfind("complex int")==0)
	{
		str.erase(0,12);
		if (!str.empty())
		{
			istringstream iss(str);
			string var_name, func;
			ld minn, maxn;
			string a,b;
			iss >> var_name >> a >> b;
			if(isdigit(var_name[0])){
				if (lagg == 1)out("The variable name cannot be used in Numbers at first.\n", YELLOW);
				else if (lagg == 2)out("变量名开头不能用数字。。\n", YELLOW);
				else if (lagg == 3)out("Переменные не могут начинаться с чисел.\n", YELLOW);
				return 0;
			}
			minn = calc(a);
			maxn = calc(b);
			getline(iss, func);
			cld bl;
			bool bao=false;
			if (complex_variables.find(var_name)!=complex_variables.end())bl=complex_variables[var_name],bao=1;
			string z,var;
			
			if (lagg==1)cout << "Complex function integration requires specifying the integration path: ";
			else if (lagg==2)cout << "复变函数积分需要指定积分路径(变量名 参数化)：";
			else if (lagg==3)cout << "Для интегрирования комплексной функции необходимо указать путь интегрирования: ";
			cin >> var; 
			getline(cin, z); // 读取剩余路径表达式
			z.erase(0, z.find_first_not_of(" \t")); // 去除前导空格
			if (z.empty())return 0;
			cld bl1;
			bool bao1=false;
			if (complex_variables.find(var)!=complex_variables.end())bl=complex_variables[var],bao1=true;
			
			auto start = chrono::high_resolution_clock::now();
			cld ans = complex_integral(func,z,var,minn,maxn,var_name);
			auto end = chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<std::chrono::microseconds>(end - start);
			if (lagg == 1)cout << "Result = ";
			else if (lagg == 2)cout << "结果 = ";
			else if (lagg == 3)cout << "Результат = ";
			cout << fixed << setprecision(precision) << ans.real() << (ans.imag()<0?'-':'+') << abs(ans.imag()) << "i\n";
			complex_variables["ans"]=ans;
			if (timing) cout << "Use:" << duration.count() / 1000 << " ms." << "\n";
			if (bao)complex_variables[var_name]=bl;
			else complex_variables.erase(var_name);
			if (bao1)complex_variables[var]=bl1;
			else complex_variables.erase(var);
			Achievement[3].first.second++;
		}
		return 0;
	}
	
	/*if (str.rfind("diff")==0){
		str.erase(0,5);
		if(!str.empty()){
			istringstream iss(str);
			string var_name, func;
			string input; 
			ld x;
			iss>>var_name>>input;
			x=calc(input);
			getline(iss,func);
			string bl;
			if(!isalpha(var_name[0])){
				if (lagg == 1)out("The variable name cannot be used in Numbers at first.\n", YELLOW);
				else if (lagg == 2)out("变量名开头不能用数字。。\n", YELLOW);
				else if (lagg == 3)out("Переменные не могут начинаться с чисел.\n", YELLOW);
				return 0;
			} 
			bool bao=0;
			if(variables.find(var_name)!=variables.end())bl = variables[var_name], bao=1;
			auto start = chrono::high_resolution_clock::now();
			ld sum=diff(func, x, var_name);
			if (lagg == 1)cout << "Result = ";
			else if (lagg == 2)cout << "结果 = ";
			else if (lagg == 3)cout << "Результат = ";
			cout << fixed << setprecision(precision) << sum << "\n";
			auto end = chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<std::chrono::microseconds>(end - start);
			if (timing) cout << "Use:" << duration.count() / 1000 << " ms." << "\n";
			variables["ans"] = to_string(sum);
			complex_variables["ans"]=cld(sum);
			if(bao)variables[var_name] = bl;
			else variables.erase(var_name);
			Achievement[4].first.second++;
			return 0;
		}
	}*/

	if (str.rfind("diff")==0)
	{
		str.erase(0,5);
		if(!str.empty())
		{
			istringstream iss(str);
			string var_name, func;
			string input;
			cld z;
			iss>>var_name>>input;
			z=complex_calc(input);
			getline(iss,func);
			cld bl;
			if(isdigit(var_name[0]))
			{
				if (lagg == 1)out("The variable name cannot be used in Numbers at first.\n", YELLOW);
				else if (lagg == 2)out("变量名开头不能用数字。。\n", YELLOW);
				else if (lagg == 3)out("Переменные не могут начинаться с чисел.\n", YELLOW);
				return 0;
			}
			bool bao=0;
			if (complex_variables.find(var_name)!=complex_variables.end())bl=complex_variables[var_name], bao=1;
			auto start = chrono::high_resolution_clock::now();
			complex_variables[var_name]=z;
			string result = diff_calc_fu(func,var_name);
			cout << "求导结果（可能未完全化简）：" << result << "\n";
			cld ans = complex_calc(result);
			auto end = chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<std::chrono::microseconds>(end - start);
			if (lagg == 1)cout << "Result = ";
			else if (lagg == 2)cout << "结果 = ";
			else if (lagg == 3)cout << "Результат = ";
			cout << fixed << setprecision(precision) << ans.real() << (ans.imag()<0?'-':'+') << abs(ans.imag()) << "i\n";
			complex_variables["ans"]=ans;
			if (timing) cout << "Use:" << duration.count() / 1000 << " ms." << "\n";
			if (bao)complex_variables[var_name]=bl;
			else complex_variables.erase(var_name);
			Achievement[4].first.second++;
		}
		return 0;
	}
	
	if (str.rfind("sum")==0){
		str.erase(0,4);
		istringstream iss(str);
		string var_name, func;
		string a,b;
		iss >> var_name >> a >> b;
		ll minn, maxn;
		minn = calc(a);
		maxn = calc(b);
		getline(iss, func);
		string bl;
		if(!isalpha(var_name[0])){
			if (lagg == 1)out("The variable name cannot be used in Numbers at first.\n", YELLOW);
			else if (lagg == 2)out("变量名开头不能用数字。。\n", YELLOW);
			else if (lagg == 3)out("Переменные не могут начинаться с чисел.\n", YELLOW);
			return 0;
		} 
		bool bao=0;
		if(variables.find(var_name)!=variables.end())bl = variables[var_name], bao=1;
		auto start = chrono::high_resolution_clock::now();
		auto f = [&](ld x) { 
		    variables[var_name] = to_string(x);
		    return calc(func); 
		};
		ld ans=0;
		rep(i,minn,maxn,1)ans+=f(i);
		if (lagg == 1)cout << "Result = ";
		else if (lagg == 2)cout << "结果 = ";
		else if (lagg == 3)cout << "Результат = ";
		cout << fixed << setprecision(precision) << ans << "\n";
		variables["ans"]=to_string(ans);
		complex_variables["ans"]=cld(ans);
		auto end = chrono::high_resolution_clock::now();
		auto duration = chrono::duration_cast<std::chrono::microseconds>(end - start);
		if (timing) cout << "Use:" << duration.count() / 1000 << " ms." << "\n";
		if(bao)variables[var_name]=bl;
		else variables.erase(var_name);
		return 0;
	}
	
	if (str.rfind("data")==0){
		str.erase(0,5);
		istringstream iss(str);
		string var_name, func;
		string a,b,li;
		iss>>var_name>>a>>b>>li;
		ld minn, maxn, l;
		minn=calc(a);
		maxn=calc(b);
		if(minn>maxn)swap(minn,maxn);
		l=calc(li);
		getline(iss, func);
		string bl;
		if(!isalpha(var_name[0])){
			if (lagg == 1)out("The variable name cannot be used in Numbers at first.\n", YELLOW);
			else if (lagg == 2)out("变量名开头不能用数字。。\n", YELLOW);
			else if (lagg == 3)out("Переменные не могут начинаться с чисел.\n", YELLOW);
			return 0;
		} 
		bool bao=0;
		if(variables.find(var_name)!=variables.end())bl = variables[var_name], bao=1;
		auto f = [&](ld x) { 
		    variables[var_name] = to_string(x);
		    return calc(func); 
		};
		
		vector<string> csv;
		ld ans=0;
		cout<<var_name<<"  result\n";
		csv.push_back(var_name+"  result\n");
		ld time=0;
		rep(i,minn,maxn,l){
			auto start = chrono::high_resolution_clock::now();
			ans = f(i);
			auto end = chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<std::chrono::microseconds>(end - start);
			time += duration.count() / 1000;
			cout << i << "  " << ans << "\n";
			csv.push_back(to_string(i)+"  "+to_string(ans)+"\n");
		}
		
		if (timing) cout << "Use:" << time << " ms." << "\n";
		if(bao)variables[var_name]=bl;
		else variables.erase(var_name);
		Achievement[13].first.second++;
		return 0;
	}

	if (str.rfind("solve") == 0) {
	    str.erase(0, 6);
	    if (!str.empty()) {
		    constexpr ld eps = 1e-8;
	        constexpr int max_iter = 500;
	        
	        istringstream iss(str);
	        string var_name, x0_str, func;
	        ld x0, x1;
	        
	        iss >> var_name >> x0_str;
	        getline(iss, func);
	        func.erase(0, func.find_first_not_of(" \t"));
	    	
	        x0 = calc(x0_str);
	    	
	        string old_value;
	        bool var_existed = (variables.find(var_name) != variables.end());
	        if (var_existed) old_value = variables[var_name];
	        
	        auto start = chrono::high_resolution_clock::now();
	        int iter = 0;
	        bool converged = false;
	        
	        auto f = [&](ld x) {
	            variables[var_name] = to_string(x);
	            return calc(func);
	        };
	        
	        while (iter++ < max_iter) {
	            ld fx = f(x0);
	        	complex_variables[var_name]=x0;
	            ld dfx = complex_calc(diff_calc_fu(func, var_name)).real();
	            
	            // 处理导数为零的情况
	            if (abs(dfx) < 1e-10) {
	                if (lagg == 1) cout << "Warning: Zero derivative at x = " << x0 << endl;
	                else if (lagg == 2) cout << "警告：导数在 x = " << x0 << " 处为零" << '\n';
	                else if (lagg == 3) cout << "Предупреждение: Нулевая производная при x = " << x0 << '\n';
	                break;
	            }
	            
	            x1 = x0 - fx / dfx;
	            
	            // 收敛检查
	            if (abs(x1 - x0) < eps && abs(fx) < eps) {
	                converged = true;
	                break;
	            }
	            
	            x0 = x1;
	        }
	        
	        if (!converged) {
	            if (lagg == 1) cout << "Warning: Did not converge after " << max_iter << " iterations\n";
	            else if (lagg == 2) cout << "警告：经过 " << max_iter << " 次迭代仍未收敛\n";
	            else if (lagg == 3) cout << "Предупреждение: Не сошлось после " << max_iter << " итераций\n";
	        }
	        
	        if (lagg == 1) cout << "Solution: ";
	        else if (lagg == 2) cout << "解：";
	        else if (lagg == 3) cout << "Решение: ";
	        cout << fixed << setprecision(precision) << x0 << "\n";
	    	complex_variables["ans"]=cld(x0);
	    	variables["ans"]=to_string(x0);
	        if (var_existed) variables[var_name] = old_value;
	        else variables.erase(var_name);
	        
	        auto end = chrono::high_resolution_clock::now();
	        auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
	        if (timing) {
	            if (lagg == 1) cout << "Time: ";
	            else if (lagg == 2) cout << "耗时：";
	            else if (lagg == 3) cout << "Время: ";
	            cout << duration.count() / 1000.0 << " ms\n";
	        }
	    }
		Achievement[12].first.second++;
	    return 0;
	}
	
	if (str == "science off") {
		if (lagg == 1)out("The Scientific notation mode has been turned off.\n", "\033[90m");
		else if (lagg == 2)out("已关闭科学计数法模式。\n", "\033[90m");
		else if (lagg == 3)out("Научный метод подсчёта закрыт.\n", "\033[90m");
		science = 0;
		return 0;
	}
	
	if (str.find("pow ") == 0) {
		big_pow(str);
		Achievement[14].first.second++;
		return 0;
	}

	if (str.rfind("science", 0) == 0) {
		str.erase(0, 8);
		if (!str.empty()) {
			bool f = 1;
			for (char c : str)if (!isdigit(c)) {
					f = 0;
					break;
				}
			if (f) {
				baoliu = stoi(str);
				if (lagg == 1)out("The number of reserved digits is adjusted to " + to_string(baoliu) + "\n", CYAN);
				else if (lagg == 2)out("保留位数调整为 " + to_string(baoliu) + "\n", CYAN);
				else if (lagg == 3)out("Оговорка была изменена на" + to_string(baoliu) + "\n", CYAN);
			}
		}
		if (lagg == 1)out("The Scientific Notation mode has been enabled. Note: It can only be used in high-precision mode and will not display the full digits.\n", BLUE);
		else if (lagg == 2)out("已开启科学计数法模式，注意：只能在高精度模式下使用，且会不显示完整数位。\n", BLUE);
		else if (lagg == 3)out("Включен режим научного счёта, заметьте: его можно использовать только в режиме высокой точности, и он не будет отображать полный разряд.\n", BLUE);
		science = 1;
		return 0;
	}

	if (str == "clear local history") {
		ofstream outFile("History.txt", ios::trunc);
		if (!outFile) {
			if (lagg == 1)out("Error: Unable to clear history file.\n", RED);
			else if (lagg == 1)out("无法清除历史文件。\n", RED);
			else if (lagg == 3)out("Не могу стереть исторические документы.\n", RED);
			return 0;
		}
		outFile.close();
		if (lagg == 1)out("Local history file (History.txt) has been cleared.\n", BLUE);
		else if (lagg == 2)out("记录已清除~\n", CYAN);
		else if (lagg == 3)out("Локальный файл истории (History.txt) был очищен.\n", BLUE);
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
			if (lagg == 1)out("Error: Unable to open file '" + path + "' for writing.\n", RED);
			else if (lagg == 2)out("你小子又把文件放哪了？？？\n", RED);
			else if (lagg == 3)out("Ошибка: Не удалось открыть файл '" + path + "' для записи.\n", RED);
			return 0;
		}

		time_t now = time(0);
		char* dt = ctime(&now);
		if (lagg == 1)outFile << "\n=== History saved on: " << dt;
		else if (lagg == 2)outFile << "\n===保存记录的时间：" << dt;
		else if (lagg == 3)outFile << "\n===Время записи:" << dt;

		for (size_t i = 0; i < history.size(); ++i) {
			outFile << i + 1 << ". " << history[i].in << " = "
			        << fixed << setprecision(precision) << history[i].result << endl;
		}

		outFile.close();
		if (lagg == 1)out("History saved to '" + path + "' (" + to_string(history.size()) + " records)\n", GREEN);
		else if (lagg == 2)out("历史记录存于" + path + "，共" + to_string(history.size()) + "条。\n", GREEN);
		else if (lagg == 3)out("История сохранена в '" + path + "' (" + to_string(history.size()) + " записей)\n", GREEN);
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
			if (lagg == 1)out("No history file found.\n", RED);
			else if (lagg == 2)out("你文件放哪了？？？？？\n", RED);
			else if (lagg == 3)out("Исторических документов нет.\n", RED);
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
			if (!history.empty())rep(i, 0, (history.size() - 1), 1) cout << i + 1 << ": " << history[i].in << " = " << history[i].result << "\n";
			Achievement[8].first.second++;
			return 0;
		}
		int xv = calc(str);
		if (xv > history.size() || xv < 1) {
			if (lagg == 1)out("Out of Range.\n", RED);
			else if (lagg == 2)out("超出范围。\n", RED);
			else if (lagg == 3)out("Вне диапазона.\n", RED);
			return 0;
		}
		if (xv != 0)cout << history[xv - 1].in << " = " << history[xv - 1].result << "\n";
		else {
			if (lagg == 1)out("The expression is incorrect.\n", RED);
			else if (lagg == 2)out("表达式错误。\n", RED);
			else if (lagg == 3)out("Ошибка выражения.\n", RED);
		}
		Achievement[8].first.second++;
		return 0;
	}

	if (str == "help") {
		if (lagg == 1) {
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
		} else if (lagg == 2) {
			cout << CYAN << "MathX 命令大全\n"
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
		} else if (lagg == 3) {
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

		if (isdigit(var_name[0])) {
			if (lagg == 1)out("The variable name cannot be used in Numbers at first.\n", YELLOW);
			else if (lagg == 2)out("变量名开头不能用数字。。\n", YELLOW);
			else if (lagg == 3)out("Переменные не могут начинаться с чисел.\n", YELLOW);
			return 0;
		}
		ld result1;
		cld result;
		if (!complex_mode)result1 = calc(value_str);
		else result = complex_calc(value_str);
		ostringstream oss;
		if (!complex_mode)oss << fixed << setprecision(80) << result1;
		else oss << fixed << setprecision(80) << result;
		if (!complex_mode)variables[var_name] = oss.str();
		complex_variables[var_name] = result;

		if (lagg == 1)
		{
			out("Variable '" + var_name + "' set to ", BLUE);
			cout<<BLUE<<complex_variables[var_name]<<"\n"<<RESET;
		}
		else if (lagg == 2)out("变量" + var_name + "已经设为" + variables[var_name] + "\n", BLUE);
		else if (lagg == 3)out("Variable '" + var_name + "' set to " + variables[var_name] + "\n", BLUE);
		Achievement[7].first.second++;
		return 0;
	}

	if (str == "list vars") {
		int x = 1;
		for (auto v : complex_variables)
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

		if (complex_variables.find(str) == complex_variables.end()) {
			if (lagg == 1) out("Variable '" + str + "' not found.\n", YELLOW);
			else if (lagg == 2) out("变量 '" + str + "' 不存在。\n", YELLOW);
			else if (lagg == 3) out("Переменная '" + str + "' не найдена.\n", YELLOW);
			return 0;
		}
		if (complex_variables.find(str) != complex_variables.end())variables.erase(str);
		complex_variables.erase(str);
		if (lagg == 1) out("Deleted '" + str + "'\n", GREEN);
		else if (lagg == 2) out("已删除变量 '" + str + "'\n", GREEN);
		else if (lagg == 3) out("Переменная '" + str + "' удалена.\n", GREEN);
		return 0;
	}

	if (str == "clear vars") {
		variables.clear();
		complex_variables.clear();
		if (lagg == 1)out("Finished\n", GREEN);
		else if (lagg == 2)out("已完成。\n", GREEN);
		else if (lagg == 3)out("Готово.\n", GREEN);
		return 0;
	}

	if (str == "clear window" || str == "clear") {
		clear_screen();
		return 0;
	}

	if (str == "timing on" || str == "timing") {
		timing = 1;
		if (lagg == 1)out("Enable timing.\n", CYAN);
		else if (lagg == 2)out("已开启计时模式。\n", CYAN);
		else if (lagg == 3)out("Таймер активирован.\n", CYAN);
		return 0;
	}

	if (str == "timing off") {
		timing = 0;
		if (lagg == 1)out("timing disabled.\n", "\033[90m");
		else if (lagg == 2)out("已关闭计时模式。\n", "\033[90m");
		else if (lagg == 3)out("Режим таймера выключен.\n", "\033[90m");
		return 0;
	}

	if (str == "rpn on" || str == "RPN on" || str == "rpn") {
		rpn = 1;
		if (lagg == 1)out("Enable RPN.\n", CYAN);
		else if (lagg == 2)out("已开启RPN模式。\n", CYAN);
		else if (lagg == 3)out("RPN включён.\n", CYAN);
		return 0;
	}

	if (str == "rpn off" || str == "RPN off") {
		rpn = 0;
		if (lagg == 1)out("RPN disabled.\n", "\033[90m");
		else if (lagg == 2)out("已关闭RPN模式。\n", "\033[90m");
		else if (lagg == 3)out("RPN отключён.\n", "\033[90m");
		return 0;
	}

	if (str == "high precision on" || str == "high precision" || str == "high") {
		high_precision = 1;
		if (lagg == 1)out("Opened high precision.\nNote: Only addition (+), subtraction (-), multiplication (*), and division (/), power(^)can be used in high precision mode, and only integers can be input and output in this mode.\n", CYAN);
		else if (lagg == 2)out("已开启高精度模式。\n注意：高精度模式下仅支持加法（+）、减法（-）、乘法（*）和除法（/）、次方（^）运算，且该模式下仅允许输入和输出整数。\n", CYAN);
		else if (lagg == 3)out("Активирован режим высокой точности.\nПримечание: В этом режиме доступны только операции сложения (+), вычитания (-), умножения (*) и деления (/),Степен (^), а ввод и вывод возможны только в целых числах.\n", CYAN);
		return 0;
	}

	if (str == "high precision off" || str == "high off") {
		high_precision = 0;
		if (lagg == 1)out("Closed high precision.\n", "\033[90m");
		else if (lagg == 2)out("已关闭高精度模式。\n", "\033[90m");
		else if (lagg == 3)out("Режим высокой точности отключен.\n", "\033[90m");
		return 0;
	}

	if (str.rfind("language", 0) == 0) {
		str.erase(0, 8);
		str.erase(0, str.find_first_not_of(" \t"));

		if (str.empty()) {
			if (lagg == 1) {
				out("Usage: language [en|zh|ru]\n", MAGENTA);
				out("Example: language zh (Switch to Chinese)\n", MAGENTA);
			} else if (lagg == 2) {
				out("用法：language [en|zh|ru]\n", MAGENTA);
				out("比如：language ru (切换到俄文)\n", MAGENTA);
			} else if (lagg == 3) {
				out("Синтаксис: language [en|zh|ru]\n", MAGENTA);
				out("Например: language en (переключить на английский)\n", MAGENTA);
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
		} else {
			if (lagg == 1)out("Unsupported language. Available options: en, zh, ru\n", RED);
			else if (lagg == 2) out("这是不支持的语言。可选选项：en, zh, ru", RED);
			else if (lagg == 3) out("Неподдерживаемый язык. Доступные варианты: en, zh, ru.\n", RED);
		}
		inledge();
		return 0;
	}
	
	return 1;
}

static void run() {
	is_achievement();
	cout << ">> " << flush;
	string str;
	getline(cin, str);

	if (!chao(str))return;

	bool f = true;
	for (const char c : str) {
		if (isalpha(c) || isdigit(c)) {
			f = false;
			break;
		}
	}
	if (f)return;

	const auto start = chrono::high_resolution_clock::now();

	if (!high_precision&&!complex_mode) {
		ld result = calc(str);
		auto end = chrono::high_resolution_clock::now();
		auto duration = chrono::duration_cast<std::chrono::microseconds>(end - start);
		if (lagg == 1)cout << "Result = ";
		else if (lagg == 2)cout << "结果 = ";
		else if (lagg == 3)cout << "Результат = ";
		cout << fixed << setprecision(precision) << result << "\n";
		//        cout << result << endl;
		history.push_back({str, to_string(result)});
		variables["ans"] = to_string(result);
		if (timing)cout << "Use:" << duration.count() / 1000 << " ms." << "\n";
		Achievement[2].first.second++;
	} else if (high_precision) {
		string result = high_precision_calc(str);
		auto end = chrono::high_resolution_clock::now();
		auto duration = chrono::duration_cast<std::chrono::microseconds>(end - start);
		if (lagg == 1)cout << "Result = ";
		else if (lagg == 2)cout << "结果 = ";
		else if (lagg == 3)cout << "Результат = ";
		if (science) {
			if (result.size() < baoliu) {
				cout << result << endl;
				history.push_back({str, result});
			} else {
				string jia = "";
				cout << result[0] << '.';
				jia += result[0];
				jia += '.';
				rep(i, 1, (baoliu - 1), 1)cout << result[i], jia += result[i];
				cout << "e" << result.size() - 1 << "\n";
				jia += "e";
				jia += to_string(result.size() - 1);
				history.push_back({str, result});
			}
		} else {
			cout << result << "\n";
			history.push_back({str, result});
		}
		if (timing)cout << "Use:" << duration.count() / 1000 << " ms." << "\n";
		Achievement[9].first.second++;
	}else if (complex_mode)
	{
		cld result = complex_calc(str);
		auto end = chrono::high_resolution_clock::now();
		auto duration = chrono::duration_cast<std::chrono::microseconds>(end - start);
		if (lagg == 1)cout << "Result = ";
		else if (lagg == 2)cout << "结果 = ";
		else if (lagg == 3)cout << "Результат = ";
		cout << fixed << setprecision(precision) << result.real();
		cout<<(result.imag()<0 ? " - " : " + ")<<abs(result.imag())<<"i\n";
		history.push_back({str, to_string(result.real())+"+"+to_string(result.imag())+"i"});
		complex_variables["ans"] = result;
		if (timing)cout << "Use:" << duration.count() / 1000 << " ms." << "\n";
		Achievement[10].first.second++;
	}
}

int main() {
	ios::sync_with_stdio(false);
	cin.tie(nullptr);
	cout.tie(nullptr);
	Ready();
	while (true) run();
}