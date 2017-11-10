#!/usr/bin/env python
import sys,os
import numpy as np
import time
import fnmatch
import cPickle 
from operator import mul
import pylab as py
from line_profiler import LineProfiler

def timer(f,*args):
	t0=time.time()
	f(*args)
	t1=time.time()
	print 't = %5.0e '%(t1-t0),f.func_name
	return t1-t0

def checkdir(path):
	if not os.path.exists(path): 
		os.makedirs(path)

def tex(x):
	return r'$\mathrm{'+x+'}$'

def save(data,name):	
	f=open(name,"w")
	cPickle.dump(data, f)
	f.close()

def load(name):	
	f=open(name,"r")
	data=cPickle.load(f)
	f.close()
	return data

def isnumeric(value):
	try:
		int(value)
		return True
	except:
		return False

	return r'$\mathrm{'+x+'}$'

def get_rgba(r,g,b,a):
	r/=255.0
	g/=255.0
	b/=255.0
	return (r,g,b,a)

def get_rgb(r,g,b):
	r/=255.0
	g/=255.0
	b/=255.0
	return (r,g,b)

def colors(name):

	if name=='maroon':return get_rgb(128,0,0)
	elif name=='dark red':return get_rgb(139,0,0)
	elif name=='brown':return get_rgb(165,42,42)
	elif name=='firebrick':return get_rgb(178,34,34)
	elif name=='crimson':return get_rgb(220,20,60)
	elif name=='red':return get_rgb(255,0,0)
	elif name=='tomato':return get_rgb(255,99,71)
	elif name=='coral':return get_rgb(255,127,80)
	elif name=='indian red':return get_rgb(205,92,92)
	elif name=='light coral':return get_rgb(240,128,128)
	elif name=='dark salmon':return get_rgb(233,150,122)
	elif name=='salmon':return get_rgb(250,128,114)
	elif name=='light salmon':return get_rgb(255,160,122)
	elif name=='orange red':return get_rgb(255,69,0)
	elif name=='dark orange':return get_rgb(255,140,0)
	elif name=='orange':return get_rgb(255,165,0)
	elif name=='gold':return get_rgb(255,215,0)
	elif name=='dark golden rod':return get_rgb(184,134,11)
	elif name=='golden rod':return get_rgb(218,165,32)
	elif name=='pale golden rod':return get_rgb(238,232,170)
	elif name=='dark khaki':return get_rgb(189,183,107)
	elif name=='khaki':return get_rgb(240,230,140)
	elif name=='olive':return get_rgb(128,128,0)
	elif name=='yellow':return get_rgb(255,255,0)
	elif name=='yellow green':return get_rgb(154,205,50)
	elif name=='dark olive green':return get_rgb(85,107,47)
	elif name=='olive drab':return get_rgb(107,142,35)
	elif name=='lawn green':return get_rgb(124,252,0)
	elif name=='chart reuse':return get_rgb(127,255,0)
	elif name=='green yellow':return get_rgb(173,255,47)
	elif name=='dark green':return get_rgb(0,100,0)
	elif name=='green':return get_rgb(0,128,0)
	elif name=='forest green':return get_rgb(34,139,34)
	elif name=='lime':return get_rgb(0,255,0)
	elif name=='lime green':return get_rgb(50,205,50)
	elif name=='light green':return get_rgb(144,238,144)
	elif name=='pale green':return get_rgb(152,251,152)
	elif name=='dark sea green':return get_rgb(143,188,143)
	elif name=='medium spring green':return get_rgb(0,250,154)
	elif name=='spring green':return get_rgb(0,255,127)
	elif name=='sea green':return get_rgb(46,139,87)
	elif name=='medium aqua marine':return get_rgb(102,205,170)
	elif name=='medium sea green':return get_rgb(60,179,113)
	elif name=='light sea green':return get_rgb(32,178,170)
	elif name=='dark slate gray':return get_rgb(47,79,79)
	elif name=='teal':return get_rgb(0,128,128)
	elif name=='dark cyan':return get_rgb(0,139,139)
	elif name=='aqua':return get_rgb(0,255,255)
	elif name=='cyan':return get_rgb(0,255,255)
	elif name=='light cyan':return get_rgb(224,255,255)
	elif name=='dark turquoise':return get_rgb(0,206,209)
	elif name=='turquoise':return get_rgb(64,224,208)
	elif name=='medium turquoise':return get_rgb(72,209,204)
	elif name=='pale turquoise':return get_rgb(175,238,238)
	elif name=='aqua marine':return get_rgb(127,255,212)
	elif name=='powder blue':return get_rgb(176,224,230)
	elif name=='cadet blue':return get_rgb(95,158,160)
	elif name=='steel blue':return get_rgb(70,130,180)
	elif name=='corn flower blue':return get_rgb(100,149,237)
	elif name=='deep sky blue':return get_rgb(0,191,255)
	elif name=='dodger blue':return get_rgb(30,144,255)
	elif name=='light blue':return get_rgb(173,216,230)
	elif name=='sky blue':return get_rgb(135,206,235)
	elif name=='light sky blue':return get_rgb(135,206,250)
	elif name=='midnight blue':return get_rgb(25,25,112)
	elif name=='navy':return get_rgb(0,0,128)
	elif name=='dark blue':return get_rgb(0,0,139)
	elif name=='medium blue':return get_rgb(0,0,205)
	elif name=='blue':return get_rgb(0,0,255)
	elif name=='royal blue':return get_rgb(65,105,225)
	elif name=='blue violet':return get_rgb(138,43,226)
	elif name=='indigo':return get_rgb(75,0,130)
	elif name=='dark slate blue':return get_rgb(72,61,139)
	elif name=='slate blue':return get_rgb(106,90,205)
	elif name=='medium slate blue':return get_rgb(123,104,238)
	elif name=='medium purple':return get_rgb(147,112,219)
	elif name=='dark magenta':return get_rgb(139,0,139)
	elif name=='dark violet':return get_rgb(148,0,211)
	elif name=='dark orchid':return get_rgb(153,50,204)
	elif name=='medium orchid':return get_rgb(186,85,211)
	elif name=='purple':return get_rgb(128,0,128)
	elif name=='thistle':return get_rgb(216,191,216)
	elif name=='plum':return get_rgb(221,160,221)
	elif name=='violet':return get_rgb(238,130,238)
	elif name=='magenta / fuchsia':return get_rgb(255,0,255)
	elif name=='orchid':return get_rgb(218,112,214)
	elif name=='medium violet red':return get_rgb(199,21,133)
	elif name=='pale violet red':return get_rgb(219,112,147)
	elif name=='deep pink':return get_rgb(255,20,147)
	elif name=='hot pink':return get_rgb(255,105,180)
	elif name=='light pink':return get_rgb(255,182,193)
	elif name=='pink':return get_rgb(255,192,203)
	elif name=='antique white':return get_rgb(250,235,215)
	elif name=='beige':return get_rgb(245,245,220)
	elif name=='bisque':return get_rgb(255,228,196)
	elif name=='blanched almond':return get_rgb(255,235,205)
	elif name=='wheat':return get_rgb(245,222,179)
	elif name=='corn silk':return get_rgb(255,248,220)
	elif name=='lemon chiffon':return get_rgb(255,250,205)
	elif name=='light golden rod yellow':return get_rgb(250,250,210)
	elif name=='light yellow':return get_rgb(255,255,224)
	elif name=='saddle brown':return get_rgb(139,69,19)
	elif name=='sienna':return get_rgb(160,82,45)
	elif name=='chocolate':return get_rgb(210,105,30)
	elif name=='peru':return get_rgb(205,133,63)
	elif name=='sandy brown':return get_rgb(244,164,96)
	elif name=='burly wood':return get_rgb(222,184,135)
	elif name=='tan':return get_rgb(210,180,140)
	elif name=='rosy brown':return get_rgb(188,143,143)
	elif name=='moccasin':return get_rgb(255,228,181)
	elif name=='navajo white':return get_rgb(255,222,173)
	elif name=='peach puff':return get_rgb(255,218,185)
	elif name=='misty rose':return get_rgb(255,228,225)
	elif name=='lavender blush':return get_rgb(255,240,245)
	elif name=='linen':return get_rgb(250,240,230)
	elif name=='old lace':return get_rgb(253,245,230)
	elif name=='papaya whip':return get_rgb(255,239,213)
	elif name=='sea shell':return get_rgb(255,245,238)
	elif name=='mint cream':return get_rgb(245,255,250)
	elif name=='slate gray':return get_rgb(112,128,144)
	elif name=='light slate gray':return get_rgb(119,136,153)
	elif name=='light steel blue':return get_rgb(176,196,222)
	elif name=='lavender':return get_rgb(230,230,250)
	elif name=='floral white':return get_rgb(255,250,240)
	elif name=='alice blue':return get_rgb(240,248,255)
	elif name=='ghost white':return get_rgb(248,248,255)
	elif name=='honeydew':return get_rgb(240,255,240)
	elif name=='ivory':return get_rgb(255,255,240)
	elif name=='azure':return get_rgb(240,255,255)
	elif name=='snow':return get_rgb(255,250,250)
	elif name=='black':return get_rgb(0,0,0)
	elif name=='dim gray':return get_rgb(105,105,105)
	elif name=='gray':return get_rgb(128,128,128)
	elif name=='dark gray':return get_rgb(169,169,169)
	elif name=='silver':return get_rgb(192,192,192)
	elif name=='light gray':return get_rgb(211,211,211)
	elif name=='gainsboro':return get_rgb(220,220,220)
	elif name=='white smoke':return get_rgb(245,245,245)
	elif name=='white':return get_rgb(255,255,255)

def add_subplot_axes(ax,rect,axisbg='w'):
	fig = py.gcf()
	box = ax.get_position()
	width = box.width
	height = box.height
	inax_position  = ax.transAxes.transform(rect[0:2])
	transFigure = fig.transFigure.inverted()
	infig_position = transFigure.transform(inax_position)    
	x = infig_position[0]
	y = infig_position[1]
	width *= rect[2]
	height *= rect[3]  # <= Typo was here
	subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
	x_labelsize = subax.get_xticklabels()[0].get_size()
	y_labelsize = subax.get_yticklabels()[0].get_size()
	x_labelsize *= rect[2]**0.5
	y_labelsize *= rect[3]**0.5
	subax.xaxis.set_tick_params(labelsize=x_labelsize)
	subax.yaxis.set_tick_params(labelsize=y_labelsize)
	return subax

class SPLIT_AX(object):

	def __init__(self,ax):
		self.ax=ax
		self.get_LR(ax)
		self.is_ylim_set=False

	def get_LR(self,ax,xlabel='x-label',ylabel='y-label'):

		ax.axis('off')
		axisbg='w'
	
		fig = py.gcf()
		box = ax.get_position()
		transFigure = fig.transFigure.inverted()
		width = box.width/2
		height = box.height
		
		# create axL
		inax_position  = ax.transAxes.transform([0,0])
		infig_position = transFigure.transform(inax_position)    
		x = infig_position[0]
		y = infig_position[1]
		axL = fig.add_axes([x,y,width,height],axisbg=axisbg)
		axL.spines['right'].set_visible(False)
		axL.get_yaxis().tick_left()
	
		# create axR
		inax_position  = ax.transAxes.transform([0.5,0])
		infig_position = transFigure.transform(inax_position)    
		x = infig_position[0]
		y = infig_position[1]
		axR = fig.add_axes([x,y,width,height],axisbg=axisbg)
		axR.get_yaxis().tick_left()
		axR.spines['left'].set_visible(False)
		axR.axes.yaxis.set_ticklabels([])
		axR.axes.get_yaxis().set_ticks([])
	
		self.axL=axL
		self.axR=axR

	def plot(self,X,Y,*args,**kwargs):

		# break the arrays for L&R
		I=-1
		for i in range(len(X)):
			if X[i]>=0.1: 
				I=i
				break
		XL,YL=X[:I+1],Y[:I+1]
		XR,YR=X[I:],Y[I:]

		# plot arrays
		self.axR.plot(XR,YR,*args,**kwargs)
		self.axL.plot(XL,YL,*args,**kwargs) 

		# set y-limits
		y1=np.amin(Y)
		y2=np.amax(Y)

		if self.is_ylim_set==False:
			self.y1_=y1
			self.y2_=y2
			self.is_ylim_set=True
		else:
			self.y1_=np.amin([y1,self.y1_])
			self.y2_=np.amax([y2,self.y2_])

		self.axL.set_ylim(self.y1_,self.y2_)
		self.axR.set_ylim(self.y1_,self.y2_)

		# set x-limits
		self.axL.set_xlim(XL[0],0.1)
		self.axR.set_xlim(0.1,XR[-1])

		self.axR.set_xticks([0.3,0.5,0.7,0.9])
		self.axL.semilogx()

	def set_ylabel(self,text,displace=-0.15,**kwargs):
		self.axL.set_ylabel(text)
		self.axL.yaxis.set_label_coords(displace,0.5)

	def set_xlabel(self,text,displace=-0.1,**kwargs):
		self.axL.set_xlabel(text)
		self.axL.xaxis.set_label_coords(1.0,displace)

	def tick_params(self,*args,**kwargs):
		self.axL.tick_params(*args,**kwargs)
		self.axR.tick_params(*args,**kwargs)

	def set_title(self,*args,**kwargs):
		self.axL.set_title(*args,**kwargs)

	def legend(self,*args,**kwargs):
		if any([k=='loc' for k in kwargs.keys()]):
			if kwargs['loc']==1 or kwargs['loc']==4: self.axR.legend(**kwargs)
			if kwargs['loc']==2 or kwargs['loc']==3: self.axL.legend(**kwargs)
		else:
			self.axR.legend(**kwargs)

	def set_ylim(self,*args):
		self.axL.set_ylim(*args)
		self.axR.set_ylim(*args)

	def axhline(self,**kwargs):
		self.axL.axhline(**kwargs)
		self.axR.axhline(**kwargs)

class ProgressBar:
	"""
	for ipython
	"""
	def __init__(self, iterations,comment):
		self.iterations = iterations
		self.comment=comment
		self.prog_bar = '[]'
		self.fill_char = '*'
		self.width = 40
		self.__update_amount(0)
		self.animate = self.animate_ipython

	def animate_ipython(self, iter):
		print '\r', self,
		sys.stdout.flush()
		self.update_iteration(iter + 1)

	def update_iteration(self, elapsed_iter):
		self.__update_amount((elapsed_iter / float(self.iterations)) * 100.0)
		self.prog_bar += '  %d of %s complete' % (elapsed_iter, self.iterations)

	def __update_amount(self, new_amount):
		percent_done = int(round((new_amount / 100.0) * 100.0))
		all_full = self.width - 2
		num_hashes = int(round((percent_done / 100.0) * all_full))
		self.prog_bar = '[' + self.fill_char * num_hashes + ' ' * (all_full - num_hashes) + ']'+'  '+self.comment
		pct_place = (len(self.prog_bar) // 2) - len(str(percent_done))
		pct_string = '%d%%' % percent_done
		self.prog_bar = self.prog_bar[0:pct_place] + \
		(pct_string + self.prog_bar[pct_place + len(pct_string):])

	def __str__ (self):
		return str(self.prog_bar)

def do_profile(follow=[]):
	def inner(func):
		def profiled_func(*args, **kwargs):
			try:
				profiler = LineProfiler()
				profiler.add_function(func)
				for f in follow:
					profiler.add_function(f)
					profiler.enable_by_count()
					return func(*args, **kwargs)
			finally:
				profiler.print_stats()
		return profiled_func
	return inner

def fill_between(x, y1, y2=0, ax=None, **kwargs):
	"""Plot filled region between `y1` and `y2`.
	This function works exactly the same as matplotlib's fill_between, except
	that it also plots a proxy artist (specifically, a rectangle of 0 size)
	so that it can be added it appears on a legend.
	"""
	ax = ax if ax is not None else py.gca()
	ax.fill_between(x, y1, y2, **kwargs)
	if kwargs['facecolor']=='none': kwargs['facecolor']='w'
	p = py.Rectangle((0, 0), 0, 0, **kwargs)
	ax.add_patch(p)
	return p


