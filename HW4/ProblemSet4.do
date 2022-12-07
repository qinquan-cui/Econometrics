 
 
*--------------------------
* profile.do 文档
*--------------------------
*- Last update: 07/12/2022
*- By: @qinquan cui
*- email: qcui@london.edu; quark.tsui@gmail.com

*-说明：
* 此文件设定了每次启动 stata 时需要做的一些基本设定
* 你可以在此文件中添加你希望在stata启动时立刻执行的命令

*-不要自动更新
set update_query  off  // on 
				
set type double           // 设定 generate 命令产生的新变量为双精度类型
set matsize 800           // 设定矩阵的维度为 2000x2000
set scrollbufsize 2000000 // 结果窗口中显示的行数上限
set more off, perma       // 关闭分页提示符

set cformat  %4.3f  //回归结果中系数的显示格式
set pformat  %4.3f  //回归结果中 p 值的显示格式      
set sformat  %4.2f  //回归结果中 se值的显示格式     

/*
set showbaselevels off, permanently
set showemptycells off, permanently
set showomitted off, permanently
*/
set fvlabel on, permanently



*-有关这一部分的完整设定命令，请输入 help set 命令进行查看

sysdir set PLUS "/Users/cuiqinquan/Documents/Stata/ado/plus/"  // 外部命令的存放位置
sysdir set PERSONAL "/Users/cuiqinquan/Documents/Stata/ado/personal/"  // 个人文件夹位置

/*
`c(sysdir_stata)'是一个暂元，里面存放了 Stata 的安装路径：输入sysdir后显示的第一个文件路径。
例如，我的 stata17 存放于D盘根目录下，所以，`c(sysdir_stata)' = D:\stata17
*/


*采用相似的方式，可添加其它允许stata搜索的目录
*adopath + ""



* log文件：自动以当前日期为名存放于 stata15\do 文件夹下
* 若 stata1x\ 下没有 do 文件夹，则本程序会自动建立一个 
cap cd `/Users/cuiqinquan/Documents/Stata/ado/'do
if _rc{
   mkdir `/Users/cuiqinquan/Documents/Stata/ado/'do  //检测后发现无 do 文件夹，则自行建立一个
}

local fn = subinstr("`c(current_time)'",":","-",2)
local fn1 = subinstr("`c(current_date)'"," ","",3)
log    using `/Users/cuiqinquan/Documents/Stata/ado/'do\log-`fn1'-`fn'.log, text replace
cmdlog using `/Users/cuiqinquan/Documents/Stata/ado/'do\cmd-`fn1'-`fn'.log, replace




*-定义资料的全局暂元
global D "/Users/cuiqinquan/Documents/Stata/ado/personal/P218" 
*-默认路径
cd "$D"     
set scheme s2color  //设定图形格式为默认格式
  
  
*			  ==========================================
*						     Question 2   
*             ========================================== 
sysuse ccapm.dta, clear

*-set the data as time series format, i.e., let id be time
generate id = _n
tsset id


*-Q2(a): run gmm, start with beta=1, gamma=1 
gmm ({beta=1}*(cratio)^(-1*{gamma=1})*rrate - 1), instruments(L.cratio L.rrate L.e) twostep


*-Q2(b): run gmm, using option "wmatrix(hac ba 5)"
gmm ({beta=1}*(cratio)^(-1*{gamma=1})*rrate - 1), instruments(L.cratio L.rrate L.e) wmatrix(hac ba 5) twostep


*-Q2(c): test hypothesis that beta=0.98 at 95% level
display (0.998-0.98)/0.004  // t = 4.5 > 1.97 (two-tail t value at 95%)



*			  ==========================================
*						     Question 4   
*             ========================================== 
sysuse MURDER.dta, clear

*-Q4(a): run gression
reg mrdrte exec unem d90 d93


*-Q4(b): run gression
xtset id year, delta(3)  // declare the data are panel
             
list id year unem L.unem if id==1 | id==2, sepby(id)

xtreg mrdrte exec unem d90 d93, fe
mat betafe = get(_b)
mat Vfe = get(VCE)

matrix list betafe 
matrix list Vfe

/* Compared with the OLS estimators, in the FE regression: 
(1) the variance of estimators become smaller, 
(2) the sign of exec changes to be negative, 
(3) coefficient of unem becomes insignificant (changes from P =0.002 to P = 0.457)
(4) d90, d93, and constant become significant at 95%

This is because in the panel data, there exists individaul dependent effect, 
and then in the OLS, there is correlation between x_it and varepsilion_it. 
This leads to biased estimation and larger variance of estimators. 
Hence both the sign, estimation, and the variance of estimators 
(thus t statistics) can change in the FE model
*/


*-Q4(c): run regression using Frisch-Waugh and comparison
xi: reg mrdrte i.id //, noconstant
predict double rmrdrte, residuals
summarize rmrdrte

xi: reg exec i.id //, noconstant
predict double rexec, residuals
summarize rexec

xi: reg unem i.id //, noconstant
predict double runem, residuals
summarize runem

xi: reg d90 i.id //, noconstant
predict double rd90, residuals
summarize rd90

xi: reg d93 i.id //, noconstant
predict double rd93, residuals
summarize rd93


reg rmrdrte rexec runem rd90 rd93, noconstant

mat betaFW = get(_b)
matrix list betafe
matrix list betaFW

mat VFW = get(VCE)
mat Vfe_adj = VFW*((51*3-4)/(51*3-51-4))
matrix list Vfe_adj
matrix list Vfe


/* similarities: the estimators are the same in (b) and (c).
Using the Frisch-Waugh theorem, we can obtain the estimator for the Dummy Variable 
Least Squares (or fixed effects); and by algebra, we can find that the estimator 
is numerically the same as the within estimator. 

dis-similarities: (b) and (c) have different standard errors for the estimator. 
This is because when applying the Frisch-Waugh theorem, the computer miscounts the degree of freedom. This can be corrected by multiplying sqrt((nT-k)/(nT-n-k)).
*/


*-Q4(c): run RE regression 


