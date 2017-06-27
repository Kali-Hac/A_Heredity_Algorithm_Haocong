#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
☆*°☆*°(∩^o^)~━━ 2017/1/18 11:09:39        
      (ˉ▽￣～) ~~ 一捆好葱 (*˙︶˙*)☆*°
      Fuction：遗传算法√ ━━━━━☆*°☆*°
"""
import math as m
import random as r
from prettytable import PrettyTable
from decimal import*
import logging
from logging import config
import time
global s,e,length_2,f_share,f_exp_num,ori_change,dic,dif_num,ger_num,xs,max_goal,max_goals
max_goals=[]
max_goal=0
dif_num=0
xs=0
dic={}
ori_change=[]
f_share=[]
f_exp_num=[]
ger_num=0
s=0
e=0
length_2=0

def bin_to_decgoal(b):
	global ger_num,s,e,length_2
	# print 2**length_2-1
	#通过二进制(bn bn-1 bn-2...b0)的大小=s+((e-s)/2**length_2-1)*(bi*(2**i)的和)[i属于0->n]可以推出二进制表示十进制子代
	return s+float((int(b,2)/float(2**length_2-1))*(e-s))

def f(x):#定义遗传算法适应值计算函数
	return x*m.sin(10*m.pi*x)+x*3*m.cos(m.sqrt(x**10))+2.0+1/x

def select_share(f_result):
	global ger_num,f_share,f_exp_num
	num=ger_num
	sum=0
	f_share=[]
	f_exp_num=[]
	for i in f_result:
		sum=sum+abs(float(i))
	sum1=0.0
	# print
	for i in f_result:
		t=abs(float(i)/float(sum))
		t1=m.floor(round(t*100))
		t1=float(t1)/100
		f_share.append(t1)
		sum1=sum1+t1
	for i in xrange(0,len(f_share)):
		print '第'+str(i)+'号染色体比例:'+str(f_share[i])
	print sum1
	#需要做一些比例保留小数的处理，以确保总比例为1
	# share_min=min(f_share)
	# min_op=f_share.index(share_min)
	share_max=max(f_share)
	max_op=f_share.index(share_max)
	# cj=(sum1-1.0)/2.0
	# f_share[min_op]=f_share[min_op]-cj
	getcontext().prec = 2
	x=f_share[max_op]-sum1+1
	x=float(Decimal(str(x)))
	if x<0.01:
		x=0.0
	f_share[max_op]=x
	sum2=0#由于选择比例不是精确值，所以直接按得到的比例算期望数量可能会少于num,最后一个不用比例再算
	logger.info('选择操作日志  √\n')
	for i in f_share[:-1]:
		f_exp_num.append(float(Decimal(str(num*i))))
		sum2=sum2+num*i
	logger.critical('选择比例为:'+str(f_share))
	f_exp_num.append(float(Decimal(str(num-sum2))))
	logger.critical('数量期望值为:'+str(f_exp_num))

def get_sum_average_max(list):
	l=[]
	ave=sum=0
	# if list[0]>100000:
	# 	for i in xrange(0,len(list)):
	# 		list[i]=bin_to_decgoal(list[i])
	for i in list:
		sum=sum+float(i)
	l.append(sum)
	getcontext().prec=3
	ave=sum/len(list)
	ave=round(ave*(10**2))
	#稍后再来处理小数点问题
	ave=ave/100.0
	l.append(ave)
	l.append(max(list))
	return l

def Turntable_bet():
	global f_share,ori_change,ger_num,xs
	num=ger_num
	add_ger=[]
	for i in xrange(0, num):
		add_ger.append(0)
	dna_list=[]
	for i in ori_change:
		dna_list.append(i)
	ori_change=[]
	for i in xrange(0,num):
		# time.sleep(0.1)
		t=r.random()
		# print t
		s=0
		cnt=-1
		while cnt<num-1 and t>s:
			cnt=cnt+1
			s=s+f_share[cnt]
		add_ger[cnt]=add_ger[cnt]+1
		ori_change.append(dna_list[cnt])
	logger.critical('轮盘赌选择的各染色体数量分别为：'+str(add_ger))
	t=PrettyTable(['编号','初始群体','变量x','目标适应值','选择比','期望数量值','实际选择数目'])
	ori_x=[]
	ori_f=[]
	for i in ori_change:
		d=bin_to_decgoal(i)
		# t=(round(d*(10 ** xs)))
		# print 'YES',t
		# t=t/1000000.0
		ori_x.append(d)
		d = round(f(d) * (10 ** xs))
		d = d / 1000000.0
		ori_f.append(d)
	for i in xrange(0,num):
		t.add_row([str(i),str(ori_change[i]),str(ori_x[i]),str(ori_f[i]),str(f_share[i]),str(f_exp_num[i]),str(add_ger[i])])
	l1=get_sum_average_max(ori_f)
	l2=get_sum_average_max(f_share)
	l3=get_sum_average_max(f_exp_num)
	l4=get_sum_average_max(add_ger)
	t.add_row([' ',' ',' ',' ',' ',' ',' '])
	t.add_row(['总和',' ',' ',l1[0],l2[0],l3[0],l4[0]])
	t.add_row(['平均值', ' ', ' ', l1[1], l2[1], l3[1], l4[1]])
	t.add_row(['最大值', ' ', ' ', l1[2], l2[2], l3[2], l4[2]])
	print t

def exchange_dna(pc):
	global ori_change,ger_num,dic,dif_num
	new_dic={}
	num=ger_num
	dna_list=[]
	for i in ori_change:
		dna_list.append(i)
	table=PrettyTable(['编号','交配池','交叉位置','交叉后新的种群','变量x','适应值大小'])
	ex_pos=[]
	ex_pool=[]
	ex_nums=[]
	ex_succ=[]
	ex_x=[]
	ex_fit=[]
	ex_num=int(round(num*pc))
	cnt=0
	choose=[]
	choose_1=[]
	choose_2=[]
	s=str(num)+'.499999'
	s=float(s)
	# print s
	for i in xrange(0,num):
		choose.append('F')
	while cnt<ex_num:
		t=int(round(r.uniform(0.5,s))-1)
		# print t
		if choose[t]=='F':
			choose[t]='T'
			cnt=cnt+1
			if cnt%2==0:
				choose_2.append(dna_list[t])
				if dna_list[t] in ori_change:
					ori_change.remove(dna_list[t])
			elif cnt!=ex_num:
				choose_1.append(dna_list[t])
				if dna_list[t] in ori_change:
					ori_change.remove(dna_list[t])
	for i in xrange(0,len(choose_2)):
		t=int(round(r.uniform(0.5,length_2-1+0.4999999)))
		logger.info(str(choose_1[i])+'与'+str(choose_2[i])+'在第'+str(t)+'位置后交换染色体序列')#测试
		ex_pos.append(t)
		ex_pos.append(t)
		ex_1 = choose_1[i][t:]
		ex_2 = choose_2[i][t:]
		ex_nums.append(dic[choose_1[i]])
		ex_nums.append(dic[choose_2[i]])
		str_1=choose_1[i][0:t]+'|'+ex_1
		str_2 =choose_2[i][0:t] + '|' + ex_2
		ex_pool.append(str_1)
		ex_pool.append(str_2)
		choose_1[i]=choose_1[i][0:t]
		choose_2[i] = choose_2[i][0:t]
		choose_1[i]=choose_1[i]+ex_2
		choose_2[i]=choose_2[i]+ex_1
		logger.critical('交叉交换成功,变为'+str(choose_1[i])+'与'+str(choose_2[i]))#测试
		ex_succ.append(choose_1[i])
		if not dic.has_key(choose_1[i]):
			dif_num=dif_num+1
			dic[choose_1[i]]=dif_num
			new_dic[choose_1[i]] = dif_num
		ex_succ.append(choose_2[i])
		if not dic.has_key(choose_2[i]):
			dif_num=dif_num+1
			dic[choose_2[i]]=dif_num
			new_dic[choose_1[i]] = dif_num
		tt=bin_to_decgoal(choose_1[i])
		ex_x.append(tt)
		ex_fit.append(f(tt))
		tt = bin_to_decgoal(choose_2[i])
		ex_x.append(tt)
		ex_fit.append(f(tt))
		ori_change.append(choose_1[i])
		ori_change.append(choose_2[i])
	# print ex_succ
	# print ex_nums
	# print ex_pool
	# print ex_pos
	# print ex_fit
	# print ex_x
	logger.info('交叉操作日志  √\n')
	for i in xrange(0,len(ex_nums)):
		table.add_row([ex_nums[i],ex_pool[i],ex_pos[i],ex_succ[i],ex_x[i],ex_fit[i]])
		logger.warn('染色体' + str(ex_nums[i]) + '在第' +str(ex_pos[i]) + '号位置后的基因发生交叉互换')
	if new_dic.has_key(ex_succ[i]):
		logger.critical('是否是整个种群新增的染色体？->YES√ 新的染色体编号为:' + str(new_dic[ex_succ[i]]))
	else:
		logger.info('是否是整个种群新增的染色体？->NO  原已存染色体编号为:' + str(dic[ex_succ[i]]))
	if ex_fit[i] > max_goal:
		logger.critical('出现个体适应值比历代高的个体！√适应值为:f=' + str(ex_fit[i])+'\n')
	else:
		logger.info('是否适应值有发生进化？->NO  适应值为:f=' + str(ex_fit[i]))
	# print [ex_nums[i],ex_pool[i],ex_pos[i],ex_succ[i],ex_x[i],ex_fit[i]]
	print table
	# for (key,value) in new_dic.items():
	# 	logger.critical('新增染色体:'+str(key)+'编号为:'+str(value))

def Mutation(pm):
	global ger_num,ori_change,dic,dif_num,max_goal
	# dic_mu={}
	mu_num=[]
	mu_ori=[]
	mu_new=[]
	mu_x=[]
	mu_fit=[]
	table=PrettyTable(['编号','交叉后的种群','变异的位置','变异后新的种群','变量x','适应值'])
	new_dic={}
	change_dna={}
	cnt=0
	cnt_suc=0
	for one in ori_change:
		flag=False
		ori_str=''
		for i in one:
			t=r.random()
			if t<pm:
				flag=True
				ori_str=one
				l=list(one)
				change_dna[one]=[]
				change_dna[one].append(one.find(i)+1)
				l[one.find(i)]=str(1-int(i))
				one=''.join(l)
		ori_change[cnt]=one
		if flag:
			cnt_suc = cnt_suc + 1
			mu_ori.append(ori_str)
			mu_new.append(one)
			decimal=bin_to_decgoal(one)
			mu_x.append(decimal)
			mu_fit.append(f(decimal))
			ori_change[cnt]=one
			if not dic.has_key(one):
				dif_num=dif_num+1
				dic[one]=dif_num
				new_dic[one]=dif_num
			mu_num.append(dic[one])
		cnt=cnt+1
	# print cnt_suc
	# print mu_fit
	# print mu_ori
	# print mu_x
	# print mu_new
	# print mu_num
	logger.info('变异操作日志  √\n')
	for i in xrange(0,cnt_suc):#######????????????√
		s=''.join(str(change_dna[mu_ori[i]]))
		table.add_row([mu_num[i],mu_ori[i],s,mu_new[i],mu_x[i],mu_fit[i]])
		logger.warn('染色体'+str(mu_ori[i])+'在第'+s+'号位置的基因发生变异√，变为'+str(mu_new[i]))
		if new_dic.has_key(mu_new[i]):
			logger.critical('是否是整个种群新增的染色体？->YES√ 新的染色体编号为:'+str(new_dic[mu_new[i]]))
		else:
			logger.info('是否是整个种群新增的染色体？->NO  原已存染色体编号为:'+str(dic[mu_new[i]]))
		if mu_fit[i]>max_goal:
			logger.critical('出现个体适应值比历代高的个体！√适应值为:f='+str(mu_fit[i])+'\n')
		else:
			logger.info('是否适应值有发生进化？->NO  适应值为:f=' + str(mu_fit[i]))
	print table
	# for (key,value) in new_dic.items():
	# 	logger.critical('新增染色体:'+str(key)+'编号为:'+str(value))


def heredity_algorithm(f,num,interval,precision,pm,pc,n_ger):
	global ger_num,s,e,length_2,f_share,f_exp_num,ori_change,dic,dif_num,xs,max_goal,max_goals
	xs=precision
	ger_num=num
	#初始群体数量为num
	#函数定义域的区间初末位
	s=a=interval[0]
	e=b=interval[1]
	width=e-s
	#根据精确度precision(小数点后pre位)以及区间得到二进制表示范围
	pre=precision
	a=m.floor(a*(10**pre))
	b=m.floor(b*(10**pre))
	a=int(a)
	b=int(b)
	# logger.info('根据精度取得的十进制区间为'+str(a)+'和'+str(b)) #测试
	a_2=bin(a)
	b_2=bin(b)
	l_2=bin(b-a)
	logger.info('用二进制表示的区间长度为'+str(l_2)+'表示范围从'+str(a_2)+'到'+str(b_2))#测试
	length_2=len(str(l_2))-2#注意要减去0b
	# logger.info('根据精度表示一个染色体的二进制位数为'+str(length_2))#测试
	ori_g=[]
	width = float(width) * float(float(2 ** length_2) / float((e - s) * (10 ** pre)))
	# print width
	for i in xrange(0,num):
		#是否需要考虑延时随机？
		#请注意 这里分成的等分是二进制区间里面的宽度(略大于十进制)
		# print width
		temp=r.uniform(0,width)
		# print temp
		temp=m.floor(temp*(10**pre))
		temp=int(temp)
		temp=bin(temp)
		# print temp,temp[2:]
		ori_g.append(temp[2:])
		ori_change.append(temp[2:])
		if not dic.has_key(temp[2:]):
			dif_num=dif_num+1
			dic[temp[2:]]=dif_num
	#打印初始群体
	generation=0
	update_ger=[]
	update_x=[]
	update_dna=[]
	temp=0
	while generation<n_ger:
		f_result=[]
		for i in ori_change:
			# print int(i,2)
			d=bin_to_decgoal(i)
			temp1 = d
			# logger.info('二进制表示为:'+str(i)+' 表示的个体为(10进制)'+str(d))#测试
			d=round(f(d)*(10**pre))
			d=d/1000000.0
			f_result.append(d)
			temp=i
		max_index=f_result.index(max(f_result))
		logger.critical('第'+str(generation+1)+'代最优适应值为:'+str(f_result[max_index])+'\n')
		if f_result[max_index]>max_goal:
			logger.critical('相比于之前某代的最优值:'+str(max_goal)+'本代已进化得到更优值:'+str(f_result[max_index])+'\n')
			# time.sleep(50)
			max_goal=f_result[max_index]
			max_goals.append(f_result[max_index])
			update_ger.append(generation+1)
			update_dna.append(temp)
			update_x.append(bin_to_decgoal(ori_change[f_result.index(max_goal)]))
		# print f_result
		# print ori_change
		select_share(f_result)
		Turntable_bet()
		exchange_dna(pc)
		Mutation(pm)
		generation=generation+1
	table=PrettyTable(['世代数','染色体编码','变量x','适应值'])
	for i in xrange(0,len(max_goals)):
		table.add_row([update_ger[i],update_dna[i],update_x[i],max_goals[i]])
		# print ori_change
	print table
	print '经过'+str(n_ger)+'代遗传后最优适应值为:'+str(max_goal)+'\n'
	logger.critical('经过'+str(n_ger)+'代遗传后最优适应值为:'+str(max_goal)+'\n')
# num=6
# ori=[2.1,2.2,2.3,5.5,6.8,1.8]
# select_share(ori)
# Turntable_bet()
logging.config.fileConfig("d://python//haocong.conf")
logger = logging.getLogger("haocong_1")
heredity_algorithm(f,25,[0.1,2],6,0.01,0.25,50)
###