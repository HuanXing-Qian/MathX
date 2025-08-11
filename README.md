这次就直接报MathX2.0Beta的功能了，多了超多

首先，是基础计算功能，支持36个函数，分别是：
`sqrt`, `sin`, `cos`, `tan`, `log`/`ln`, `exp`, `exp2`, `log10`, `log2`，`arcsin`, `arccos`, `arctan`, `sinh`, `cosh`, `tanh`, `ceil`（向上取整）, `floor`（向下取整）, `round`（四舍五入）, `abs`（绝对值）, `erf`, `Gamma`, `sinc`(sinx/x), `lngamma`(lnGamma), `atr`（角度转弧度）, `rta`（弧度转角度）, `prime`（求第x个素数）, `isprime`（判断整数x是否是素数）, `min(x,y)`（求min(x,y)）, `max(x,y)`(求max(x,y)), `arctan2(x,y)`(求象限的正切，用于极坐标), `gcd(a,b)`（求a和b的最大公因数）, `lcm(a,b)`（求a和b的最小公倍数），`hypot(x,y)`（求两边为x、y的斜边的直角三角形），`C(a,b)`（求组合数），`rand(a,b)`（求a到b区间内的随机整数）

然后是chao()函数里的特殊命令，分别：
- exit：退出
- clear history 清空保存历史
- factor 后面跟一个2^63-1内的正整数，求质因数分解
- int x a b f(x) 求`\int_a^b f(x) dx`，这里因为一些误差，所以如果是收敛的无穷限积分请调小一点，不然就会出现`int x 0 10000 sin(x^2)`是-293.13的情况
- diff x x0 f(x) 求`f'(x)\mid_x0`
- sum x a b f(x) 求`\sum_{x=a}^b f(x)`
- data x a b l f(x) 求[a,b]区间，步长为l中f(x)分别的值。
- solve x x0 f(x) 求f(x)函数的零点，可以用来解方程（如x^-2x+1=0，但是你要确保右边是0）
- science 开启科学计数法模式，仅限高精度模式，这里可以选择在后面加一个正整数表示精度
- science off 关闭科学计数法模式
- pow a b 求a^b的估算值以及最后18位数字，很快，但是在b比较大时（如9e18）的科学计数法前半部分会出现误差，因为快到上限了
- clear local history 清空`History.txt`的数据
- save history path 把历史保存进path中（这里的path可以直接用*.txt，会保存进当前目录），path如果是空会默认保存进`History.txt`
- load history path 描述和上面差不多，把path的文本输出到控制台，不加入当前的`history`
- !n 查询第n条历史，如果n是`all`那么就会展示所有历史
- help 用处不大
- let name value 创建name变量并赋值value，如果已存在那么name的值会变为value，这里的value可以是表达式。注意，ans是自动创建的，表示上一次的计算结果。
- list vars 展示所有变量，包括ans
- delete name 删除name变量，当然别想着删掉ans，因为会自动创建
- clear vars 清空所有变量，ans还是自动创建
- clear window/clear 清空窗口，仅限windows
- timing on/timing 开启计时模式，启动时就开启
- timing off 关闭计时模式
- rpn on/RPN on/rpn 开启显示`RPN`（逆波兰表达式）模式
- rpn off 关闭RPN模式
- high precision on/high precision/high 开启高精度模式，该模式只能`+,-,*,/,^`以及整数运算
- high precision off/high off 关闭高精度模式
- language [language] 切换界面语言，支持[en|zh|ru]，在某些隐秘的地方只支持zh
- study know 展示know的知识点，提供拼写错误提供建议系统，如果know为空就展示所有领域
- search domain 展示domain中的所有知识的标题
最近2025年8月8日我又添加了非常好的一个功能——复数，输入`complex`或者`complex on`就可以打开，在该模式下你可以使用复数运算，支持了普通计算的大部分功能，如`sin`,`cos`,`tan`等，一些比大小、整数函数就不支持，因为根本就支持不了。
- complex int z a b f(z)，按下回车后继续输入z(t)的参数化的变量名t并输入z(t)的表达式，就可以计算参数化z(t)的复变积分`\int_a^b f(z(t))z'(t) dt`
- complex diff z z0 f(z)，和实数diff差不多，就是可以支持复数了。