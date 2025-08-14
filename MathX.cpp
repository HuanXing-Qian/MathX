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
		else if (lagg==2)printf ( "��ǰ��ʱ����: %s", asctime (timeinfo) );
		else if (lagg==3)printf ( "���֧ܧ��֧� �ӧ�֧ާ�: %s", asctime (timeinfo) );
		return false;
	}
	
	if (str=="complex"||str=="complex on")
	{
		complex_mode = true;
		if (lagg==1)out("Complex mode on.\n",GREEN);
		else if (lagg==2)out("����ģʽ�ѿ�����\n",GREEN);
		else if (lagg==3)out("���֧اڧ� �ާߧ�ا֧��ӧ֧ߧߧ�ԧ� ��ڧ�ݧ� �ӧܧݧ��֧�.\n",GREEN);
		return false;
	}

	if (str=="complex off")
	{
		complex_mode = 0;
		if (lagg==1)out("Complex mode off.\n","\033[90m");
		else if (lagg==2)out("����ģʽ�ѹرա�\n","\033[90m");
		else if (lagg==3)out("���ڧ��֧ާ� �ާߧ�ا֧��ӧ֧ߧߧ�� ��ڧ�֧� ���ܧݧ��֧ߧ�.\n","\033[90m");
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
			else if (lagg == 2)out("����˵����\"" + best_match + "\"��\n", YELLOW);
			else if (lagg == 3)out("���� ����֧ݧ� ��ܧѧ٧ѧ��\"" + best_match + "\"?\n", YELLOW);
		} else {
			if (lagg == 1)out("No matching concept found.\n", RED);
			else if (lagg == 2)out("δ�ҵ�ƥ���\n", RED);
			else if (lagg == 2)out("����ӧ�ѧէ֧ߧڧ� �ߧ� �ߧѧۧէ֧ߧ�.\n", RED);
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
			else if (lagg == 2)out("δ�ҵ�ƥ���\n", RED);
			else if (lagg == 2)out("����ӧ�ѧէ֧ߧڧ� �ߧ� �ߧѧۧէ֧ߧ�.\n", RED);
		}
		return 0;
	}
	if (str.empty()) return 0;
	if (str == "exit"||str=="quit") exit(0);
	if (str == "clear history") {
		if (lagg == 1)out("Finished.\n", GREEN);
		else if (lagg == 2)out("��ϡ�������Ը�Ҳ�Ҫ��˵���ơ�����ɡ������Ļ��ˡ�\n", GREEN);
		else if (lagg == 3)out("���ѧӧ֧��֧ߧ�.\n", GREEN);
		history.clear();
		return 0;
	}

	if (str.rfind("factor") == 0) {
		str.erase(0, 7);
		if (!str.empty()) {
			ll n = ll(calc(str));
			if (n <= 1) {
				cout << "���������1��������" << "\n";
				return 0;
			}
			auto factors = high_precision_func::factorize(n);
			string factorization = high_precision_func::format_factorization(factors);
			cout<<"�������ֽ⣺"<< factorization << "\n";
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
				out("��������ˡ�\n", RED);
				return 0;
			}
			if (var_name.size()>1)
			{
				out("������������һ����ĸ��\n",RED);
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
				else if (lagg == 2)out("��������ͷ���������֡���\n", YELLOW);
				else if (lagg == 3)out("���֧�֧ާ֧ߧߧ�� �ߧ� �ާ�ԧ�� �ߧѧ�ڧߧѧ���� �� ��ڧ�֧�.\n", YELLOW);
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
				else if (lagg == 2)out("��������ͷ���������֡���\n", YELLOW);
				else if (lagg == 3)out("���֧�֧ާ֧ߧߧ�� �ߧ� �ާ�ԧ�� �ߧѧ�ڧߧѧ���� �� ��ڧ�֧�.\n", YELLOW);
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
			else if (lagg==2)cout << "���亯��������Ҫָ������·��(������ ������)��";
			else if (lagg==3)cout << "���ݧ� �ڧߧ�֧ԧ�ڧ��ӧѧߧڧ� �ܧ�ާ�ݧ֧ܧ�ߧ�� ���ߧܧ�ڧ� �ߧ֧�ҧ��էڧާ� ��ܧѧ٧ѧ�� ����� �ڧߧ�֧ԧ�ڧ��ӧѧߧڧ�: ";
			cin >> var; 
			getline(cin, z); // ��ȡʣ��·�����ʽ
			z.erase(0, z.find_first_not_of(" \t")); // ȥ��ǰ���ո�
			if (z.empty())return 0;
			cld bl1;
			bool bao1=false;
			if (complex_variables.find(var)!=complex_variables.end())bl=complex_variables[var],bao1=true;
			
			auto start = chrono::high_resolution_clock::now();
			cld ans = complex_integral(func,z,var,minn,maxn,var_name);
			auto end = chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<std::chrono::microseconds>(end - start);
			if (lagg == 1)cout << "Result = ";
			else if (lagg == 2)cout << "��� = ";
			else if (lagg == 3)cout << "���֧٧�ݧ��ѧ� = ";
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
				else if (lagg == 2)out("��������ͷ���������֡���\n", YELLOW);
				else if (lagg == 3)out("���֧�֧ާ֧ߧߧ�� �ߧ� �ާ�ԧ�� �ߧѧ�ڧߧѧ���� �� ��ڧ�֧�.\n", YELLOW);
				return 0;
			} 
			bool bao=0;
			if(variables.find(var_name)!=variables.end())bl = variables[var_name], bao=1;
			auto start = chrono::high_resolution_clock::now();
			ld sum=diff(func, x, var_name);
			if (lagg == 1)cout << "Result = ";
			else if (lagg == 2)cout << "��� = ";
			else if (lagg == 3)cout << "���֧٧�ݧ��ѧ� = ";
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
				else if (lagg == 2)out("��������ͷ���������֡���\n", YELLOW);
				else if (lagg == 3)out("���֧�֧ާ֧ߧߧ�� �ߧ� �ާ�ԧ�� �ߧѧ�ڧߧѧ���� �� ��ڧ�֧�.\n", YELLOW);
				return 0;
			}
			bool bao=0;
			if (complex_variables.find(var_name)!=complex_variables.end())bl=complex_variables[var_name], bao=1;
			auto start = chrono::high_resolution_clock::now();
			complex_variables[var_name]=z;
			string result = diff_calc_fu(func,var_name);
			cout << "�󵼽��������δ��ȫ���򣩣�" << result << "\n";
			cld ans = complex_calc(result);
			auto end = chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<std::chrono::microseconds>(end - start);
			if (lagg == 1)cout << "Result = ";
			else if (lagg == 2)cout << "��� = ";
			else if (lagg == 3)cout << "���֧٧�ݧ��ѧ� = ";
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
			else if (lagg == 2)out("��������ͷ���������֡���\n", YELLOW);
			else if (lagg == 3)out("���֧�֧ާ֧ߧߧ�� �ߧ� �ާ�ԧ�� �ߧѧ�ڧߧѧ���� �� ��ڧ�֧�.\n", YELLOW);
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
		else if (lagg == 2)cout << "��� = ";
		else if (lagg == 3)cout << "���֧٧�ݧ��ѧ� = ";
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
			else if (lagg == 2)out("��������ͷ���������֡���\n", YELLOW);
			else if (lagg == 3)out("���֧�֧ާ֧ߧߧ�� �ߧ� �ާ�ԧ�� �ߧѧ�ڧߧѧ���� �� ��ڧ�֧�.\n", YELLOW);
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
	            
	            // ������Ϊ������
	            if (abs(dfx) < 1e-10) {
	                if (lagg == 1) cout << "Warning: Zero derivative at x = " << x0 << endl;
	                else if (lagg == 2) cout << "���棺������ x = " << x0 << " ��Ϊ��" << '\n';
	                else if (lagg == 3) cout << "����֧է���֧اէ֧ߧڧ�: ����ݧ֧ӧѧ� ����ڧ٧ӧ�էߧѧ� ���� x = " << x0 << '\n';
	                break;
	            }
	            
	            x1 = x0 - fx / dfx;
	            
	            // �������
	            if (abs(x1 - x0) < eps && abs(fx) < eps) {
	                converged = true;
	                break;
	            }
	            
	            x0 = x1;
	        }
	        
	        if (!converged) {
	            if (lagg == 1) cout << "Warning: Did not converge after " << max_iter << " iterations\n";
	            else if (lagg == 2) cout << "���棺���� " << max_iter << " �ε�����δ����\n";
	            else if (lagg == 3) cout << "����֧է���֧اէ֧ߧڧ�: ���� ����ݧ��� ����ݧ� " << max_iter << " �ڧ�֧�ѧ�ڧ�\n";
	        }
	        
	        if (lagg == 1) cout << "Solution: ";
	        else if (lagg == 2) cout << "�⣺";
	        else if (lagg == 3) cout << "���֧�֧ߧڧ�: ";
	        cout << fixed << setprecision(precision) << x0 << "\n";
	    	complex_variables["ans"]=cld(x0);
	    	variables["ans"]=to_string(x0);
	        if (var_existed) variables[var_name] = old_value;
	        else variables.erase(var_name);
	        
	        auto end = chrono::high_resolution_clock::now();
	        auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
	        if (timing) {
	            if (lagg == 1) cout << "Time: ";
	            else if (lagg == 2) cout << "��ʱ��";
	            else if (lagg == 3) cout << "����֧ާ�: ";
	            cout << duration.count() / 1000.0 << " ms\n";
	        }
	    }
		Achievement[12].first.second++;
	    return 0;
	}
	
	if (str == "science off") {
		if (lagg == 1)out("The Scientific notation mode has been turned off.\n", "\033[90m");
		else if (lagg == 2)out("�ѹرտ�ѧ������ģʽ��\n", "\033[90m");
		else if (lagg == 3)out("���ѧ��ߧ�� �ާ֧��� ���է��ק�� �٧ѧܧ���.\n", "\033[90m");
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
				else if (lagg == 2)out("����λ������Ϊ " + to_string(baoliu) + "\n", CYAN);
				else if (lagg == 3)out("���ԧ�ӧ��ܧ� �ҧ�ݧ� �ڧ٧ާ֧ߧ֧ߧ� �ߧ�" + to_string(baoliu) + "\n", CYAN);
			}
		}
		if (lagg == 1)out("The Scientific Notation mode has been enabled. Note: It can only be used in high-precision mode and will not display the full digits.\n", BLUE);
		else if (lagg == 2)out("�ѿ�����ѧ������ģʽ��ע�⣺ֻ���ڸ߾���ģʽ��ʹ�ã��һ᲻��ʾ������λ��\n", BLUE);
		else if (lagg == 3)out("���ܧݧ��֧� ��֧اڧ� �ߧѧ��ߧ�ԧ� ���ק��, �٧ѧާ֧����: �֧ԧ� �ާ�اߧ� �ڧ���ݧ�٧�ӧѧ�� ���ݧ�ܧ� �� ��֧اڧާ� �ӧ���ܧ�� ����ߧ����, �� ��� �ߧ� �ҧ�է֧� ����ҧ�ѧاѧ�� ���ݧߧ�� ��ѧ٧���.\n", BLUE);
		science = 1;
		return 0;
	}

	if (str == "clear local history") {
		ofstream outFile("History.txt", ios::trunc);
		if (!outFile) {
			if (lagg == 1)out("Error: Unable to clear history file.\n", RED);
			else if (lagg == 1)out("�޷������ʷ�ļ���\n", RED);
			else if (lagg == 3)out("���� �ާ�ԧ� ���֧�֧�� �ڧ����ڧ�֧�ܧڧ� �է�ܧ�ާ֧ߧ��.\n", RED);
			return 0;
		}
		outFile.close();
		if (lagg == 1)out("Local history file (History.txt) has been cleared.\n", BLUE);
		else if (lagg == 2)out("��¼�����~\n", CYAN);
		else if (lagg == 3)out("����ܧѧݧ�ߧ�� ��ѧۧ� �ڧ����ڧ� (History.txt) �ҧ�� ���ڧ�֧�.\n", BLUE);
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
			else if (lagg == 2)out("��С���ְ��ļ������ˣ�����\n", RED);
			else if (lagg == 3)out("����ڧҧܧ�: ���� ��էѧݧ��� ���ܧ���� ��ѧۧ� '" + path + "' �էݧ� �٧ѧ�ڧ��.\n", RED);
			return 0;
		}

		time_t now = time(0);
		char* dt = ctime(&now);
		if (lagg == 1)outFile << "\n=== History saved on: " << dt;
		else if (lagg == 2)outFile << "\n===�����¼��ʱ�䣺" << dt;
		else if (lagg == 3)outFile << "\n===����֧ާ� �٧ѧ�ڧ��:" << dt;

		for (size_t i = 0; i < history.size(); ++i) {
			outFile << i + 1 << ". " << history[i].in << " = "
			        << fixed << setprecision(precision) << history[i].result << endl;
		}

		outFile.close();
		if (lagg == 1)out("History saved to '" + path + "' (" + to_string(history.size()) + " records)\n", GREEN);
		else if (lagg == 2)out("��ʷ��¼����" + path + "����" + to_string(history.size()) + "����\n", GREEN);
		else if (lagg == 3)out("�������ڧ� �����ѧߧ֧ߧ� �� '" + path + "' (" + to_string(history.size()) + " �٧ѧ�ڧ�֧�)\n", GREEN);
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
			else if (lagg == 2)out("���ļ������ˣ���������\n", RED);
			else if (lagg == 3)out("�������ڧ�֧�ܧڧ� �է�ܧ�ާ֧ߧ��� �ߧ֧�.\n", RED);
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
			else if (lagg == 2)out("������Χ��\n", RED);
			else if (lagg == 3)out("���ߧ� �էڧѧ�ѧ٧�ߧ�.\n", RED);
			return 0;
		}
		if (xv != 0)cout << history[xv - 1].in << " = " << history[xv - 1].result << "\n";
		else {
			if (lagg == 1)out("The expression is incorrect.\n", RED);
			else if (lagg == 2)out("���ʽ����\n", RED);
			else if (lagg == 3)out("����ڧҧܧ� �ӧ��ѧا֧ߧڧ�.\n", RED);
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
			cout << CYAN << "MathX �����ȫ\n"
			     << "========================================\n"
			     << "��������\n"
			     << "  [���ʽ]          ������ѧ���ʽ\n"
			     << "  let x 5           �������x����ֵΪ5\n"
			     << "  list vars         ��ʾ���б���\n"
			     << "  clear vars        ���ȫ������\n\n"

			     << "�߼�����\n"
			     << "  high precision on �����߾�������ģʽ\n"
			     << "  setprecision 15   ����С��λ��(1-80)\n\n"

			     << "��ʾ����\n"
			     << "  rpn on/off         �л��沨�����ʽģʽ\n"
			     << "  timing on/off      ��ʾ/���ؼ���ʱ��\n"
			     << "  language en/zh/ru �л���������\n\n"

			     << "��ʷ��¼\n"
			     << "  !3                 ��ȡ��3��������\n"
			     << "  !all               ��ʾ������ʷ��¼\n"
			     << "  save history math.txt ������ʷ���ļ�\n\n"

			     << "ѧϰģʽ\n"
			     << "  study derivative   ѧϰ����֪ʶ\n"
			     << "  search trigonometry �г����Ǻ���\n\n"

			     << "ϵͳ����\n"
			     << "  clear              ����\n"
			     << "  exit               �˳�����\n"
			     << "========================================\n" << RESET;
		} else if (lagg == 3) {
			cout << CYAN << "�������������� ���� MATHX\n"
			     << "========================================\n"
			     << "���������������� ��������������������\n"
			     << "  [�ӧ��ѧا֧ߧڧ�]       �����ڧ�ݧڧ�� �ާѧ�֧ާѧ�ڧ�֧�ܧ�� �ӧ��ѧا֧ߧڧ�\n"
			     << "  let x 5           �����֧է֧ݧڧ�� ��֧�֧ާ֧ߧߧ�� x = 5\n"
			     << "  list vars         ����ܧѧ٧ѧ�� �ӧ�� ��֧�֧ާ֧ߧߧ��\n"
			     << "  clear vars        ���էѧݧڧ�� �ӧ�� ��֧�֧ާ֧ߧߧ��\n\n"

			     << "���������������������� ��������������\n"
			     << "  high precision on ���ܧݧ��ڧ�� ��֧اڧ� ����ڧ٧ӧ�ݧ�ߧ�� ����ߧ����\n"
			     << "  setprecision 15    �����ѧߧ�ӧڧ�� �٧ߧѧܧ� ����ݧ� �٧ѧ����� (1-80)\n\n"

			     << "������������������ ����������������������\n"
			     << "  rpn on/off        ���֧�֧ܧݧ��ڧ�� ������ (��ҧ�ѧ�ߧ�� ���ݧ��ܧ�� �٧ѧ�ڧ��)\n"
			     << "  timing on/off     ����ܧѧ٧ѧ��/��ܧ���� �ӧ�֧ާ� �ӧ��ڧ�ݧ֧ߧڧ�\n"
			     << "  language en/zh/ru ���ާ֧ߧڧ�� ��٧�� �ڧߧ�֧��֧ۧ��\n\n"

			     << "�������������� ��������������������\n"
			     << "  !3                ����ܧѧ٧ѧ�� 3-�� ��֧٧�ݧ��ѧ�\n"
			     << "  !all              ����ܧѧ٧ѧ�� �ӧ�� �ڧ����ڧ�\n"
			     << "  save history math.txt ���ܧ����� �ڧ����ڧ� �� ��ѧۧ�\n\n"

			     << "������������������ ����������\n"
			     << "  study derivative  ���٧��ڧ�� ����ڧ٧ӧ�էߧ��\n"
			     << "  search trigonometry ����ڧ��� ���ڧԧ�ߧ�ާ֧��ڧ�֧�ܧڧ� ���ߧܧ�ڧ�\n\n"

			     << "��������������\n"
			     << "  clear             ����ڧ��ڧ�� ��ܧ�ѧ�\n"
			     << "  exit              ����ۧ�� �ڧ� ����ԧ�ѧާާ�\n"
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
			else if (lagg == 2)out("��������ͷ���������֡���\n", YELLOW);
			else if (lagg == 3)out("���֧�֧ާ֧ߧߧ�� �ߧ� �ާ�ԧ�� �ߧѧ�ڧߧѧ���� �� ��ڧ�֧�.\n", YELLOW);
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
		else if (lagg == 2)out("����" + var_name + "�Ѿ���Ϊ" + variables[var_name] + "\n", BLUE);
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
			else if (lagg == 2) out("��ָ��Ҫɾ���ı�������\n", RED);
			else if (lagg == 3) out("���ܧѧاڧ�� �ڧާ� ��֧�֧ާ֧ߧߧ�� �էݧ� ��էѧݧ֧ߧڧ�.\n", RED);
			return 0;
		}

		if (complex_variables.find(str) == complex_variables.end()) {
			if (lagg == 1) out("Variable '" + str + "' not found.\n", YELLOW);
			else if (lagg == 2) out("���� '" + str + "' �����ڡ�\n", YELLOW);
			else if (lagg == 3) out("���֧�֧ާ֧ߧߧѧ� '" + str + "' �ߧ� �ߧѧۧէ֧ߧ�.\n", YELLOW);
			return 0;
		}
		if (complex_variables.find(str) != complex_variables.end())variables.erase(str);
		complex_variables.erase(str);
		if (lagg == 1) out("Deleted '" + str + "'\n", GREEN);
		else if (lagg == 2) out("��ɾ������ '" + str + "'\n", GREEN);
		else if (lagg == 3) out("���֧�֧ާ֧ߧߧѧ� '" + str + "' ��էѧݧ֧ߧ�.\n", GREEN);
		return 0;
	}

	if (str == "clear vars") {
		variables.clear();
		complex_variables.clear();
		if (lagg == 1)out("Finished\n", GREEN);
		else if (lagg == 2)out("����ɡ�\n", GREEN);
		else if (lagg == 3)out("������ӧ�.\n", GREEN);
		return 0;
	}

	if (str == "clear window" || str == "clear") {
		clear_screen();
		return 0;
	}

	if (str == "timing on" || str == "timing") {
		timing = 1;
		if (lagg == 1)out("Enable timing.\n", CYAN);
		else if (lagg == 2)out("�ѿ�����ʱģʽ��\n", CYAN);
		else if (lagg == 3)out("���ѧۧާ֧� �ѧܧ�ڧӧڧ��ӧѧ�.\n", CYAN);
		return 0;
	}

	if (str == "timing off") {
		timing = 0;
		if (lagg == 1)out("timing disabled.\n", "\033[90m");
		else if (lagg == 2)out("�ѹرռ�ʱģʽ��\n", "\033[90m");
		else if (lagg == 3)out("���֧اڧ� ��ѧۧާ֧�� �ӧ�ܧݧ��֧�.\n", "\033[90m");
		return 0;
	}

	if (str == "rpn on" || str == "RPN on" || str == "rpn") {
		rpn = 1;
		if (lagg == 1)out("Enable RPN.\n", CYAN);
		else if (lagg == 2)out("�ѿ���RPNģʽ��\n", CYAN);
		else if (lagg == 3)out("RPN �ӧܧݧ��ק�.\n", CYAN);
		return 0;
	}

	if (str == "rpn off" || str == "RPN off") {
		rpn = 0;
		if (lagg == 1)out("RPN disabled.\n", "\033[90m");
		else if (lagg == 2)out("�ѹر�RPNģʽ��\n", "\033[90m");
		else if (lagg == 3)out("RPN ���ܧݧ��ק�.\n", "\033[90m");
		return 0;
	}

	if (str == "high precision on" || str == "high precision" || str == "high") {
		high_precision = 1;
		if (lagg == 1)out("Opened high precision.\nNote: Only addition (+), subtraction (-), multiplication (*), and division (/), power(^)can be used in high precision mode, and only integers can be input and output in this mode.\n", CYAN);
		else if (lagg == 2)out("�ѿ����߾���ģʽ��\nע�⣺�߾���ģʽ�½�֧�ּӷ���+����������-�����˷���*���ͳ�����/�����η���^�����㣬�Ҹ�ģʽ�½�������������������\n", CYAN);
		else if (lagg == 3)out("���ܧ�ڧӧڧ��ӧѧ� ��֧اڧ� �ӧ���ܧ�� ����ߧ����.\n����ڧާ֧�ѧߧڧ�: �� ����� ��֧اڧާ� �է�����ߧ� ���ݧ�ܧ� ���֧�ѧ�ڧ� ��ݧ�ا֧ߧڧ� (+), �ӧ��ڧ�ѧߧڧ� (-), ��ާߧ�ا֧ߧڧ� (*) �� �է֧ݧ֧ߧڧ� (/),����֧�֧� (^), �� �ӧӧ�� �� �ӧ�ӧ�� �ӧ�٧ާ�اߧ� ���ݧ�ܧ� �� ��֧ݧ�� ��ڧ�ݧѧ�.\n", CYAN);
		return 0;
	}

	if (str == "high precision off" || str == "high off") {
		high_precision = 0;
		if (lagg == 1)out("Closed high precision.\n", "\033[90m");
		else if (lagg == 2)out("�ѹرո߾���ģʽ��\n", "\033[90m");
		else if (lagg == 3)out("���֧اڧ� �ӧ���ܧ�� ����ߧ���� ���ܧݧ��֧�.\n", "\033[90m");
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
				out("�÷���language [en|zh|ru]\n", MAGENTA);
				out("���磺language ru (�л�������)\n", MAGENTA);
			} else if (lagg == 3) {
				out("���ڧߧ�ѧܧ�ڧ�: language [en|zh|ru]\n", MAGENTA);
				out("���ѧ��ڧާ֧�: language en (��֧�֧ܧݧ��ڧ�� �ߧ� �ѧߧԧݧڧۧ�ܧڧ�)\n", MAGENTA);
			}
			return 0;
		}

		transform(str.begin(), str.end(), str.begin(), ::tolower);

		if (str == "en") {
			lagg = 1;
			out("It has been switched to English.\n", GREEN);
		} else if (str == "zh") {
			lagg = 2;
			out("���л�Ϊ���ġ�\n", GREEN);
		} else if (str == "ru") {
			lagg = 3;
			out("���֧�֧ܧݧ��ڧݧ�� �ߧ� �����ܧڧ�.\n", GREEN);
		} else {
			if (lagg == 1)out("Unsupported language. Available options: en, zh, ru\n", RED);
			else if (lagg == 2) out("���ǲ�֧�ֵ����ԡ���ѡѡ�en, zh, ru", RED);
			else if (lagg == 3) out("���֧��էէ֧�اڧӧѧ֧ާ�� ��٧��. ��������ߧ�� �ӧѧ�ڧѧߧ��: en, zh, ru.\n", RED);
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
		else if (lagg == 2)cout << "��� = ";
		else if (lagg == 3)cout << "���֧٧�ݧ��ѧ� = ";
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
		else if (lagg == 2)cout << "��� = ";
		else if (lagg == 3)cout << "���֧٧�ݧ��ѧ� = ";
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
		else if (lagg == 2)cout << "��� = ";
		else if (lagg == 3)cout << "���֧٧�ݧ��ѧ� = ";
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