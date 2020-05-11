from manimlib.imports import *
from math import *
from scipy.misc import derivative
class MyAxes(Axes):
	CONFIG = {
		"number_line_config": {
			"include_ticks": False,
		}
	}
def tlf(func,point,**kwargs):
	temp=lambda x:func(point)+derivative(func,point,dx=1e-8)*(x-point)
	return FunctionGraph(temp,**kwargs)

def get_mid(start,end,tip_angle=PI/2):
	trans=np.array([[0,-1,0],
	                [1,0,0],
	                [0,0,1]])
	angle=1/2*tip_angle
	k=1/2*(end-start)
	nk=get_norm(k)/tan(angle)
	mid=start+k+nk*normalize(k).dot(trans)
	return mid

def get_intersect(line1,line2):
	x1,y1=line1.get_start()[:-1]
	x2,y2=line1.get_end()[:-1]
	x3,y3=line2.get_start()[:-1]
	x4,y4=line2.get_end()[:-1]
	a1=y1-y2;a2=y3-y4;a3=x1-x2;a4=x3-x4
	A=[a1,a2,a3,a4]
	z=0
	if a1*a4==a2*a3:return line2.get_start()
	elif a3==0:return np.array([x1,y4+a2/a4*(x1-x4),z])
	elif a4==0:return np.array([x3,y2+a1/a3*(x3-x2),z])
	else:
		a=a1/a3-a2/a4
		b=y1-y3+a2/a4*x3-a1/a3*x1
		x=-b/a
		return np.array([x,a1/a3*(x-x1)+y1,z])

class Angle(Sector):
	CONFIG = {
		"radius": 2.5*DEFAULT_DOT_RADIUS,
        "fill_opacity": 0.5,
		"opacity":1,
	}
	def __init__(self,line1,line2,**kwargs):
		stroke_width=2/3*min(line1.get_stroke_width(),
		                     line2.get_stroke_width())
		color=line1.get_color()
		plot_depth=min(line1.plot_depth,line2.plot_depth)-1
		arc_center=get_intersect(line1,line2)
		angle=line2.get_angle()-line1.get_angle()
		start_angle=line1.get_angle()
		digest_config(self,kwargs,locals())
		Sector.__init__(self,angle=self.angle,start_angle=self.start_angle,
			            arc_center=self.arc_center,outer_radius=self.radius)

class RightAngle(Polygon):
	CONFIG = {
		"length":2.5*DEFAULT_DOT_RADIUS,
		"fill_opacity":0.5,
		"opacity":1,
	}
	def __init__(self,line1,line2,**kwargs):
		stroke_width = 2 / 3 * min(line1.get_stroke_width(),
		                           line2.get_stroke_width())
		color=line1.get_color()
		plot_depth = min(line1.plot_depth, line2.plot_depth) - 1
		arc_center = get_intersect(line1, line2)
		digest_config(self, kwargs, locals())
		self.start = self.arc_center + self.length * normalize(line1.get_vector())
		self.end = self.arc_center + self.length * (normalize(line2.get_vector())
		                                            if line1.get_vector().dot(line2.get_vector())>=0
		                                            else -normalize(line2.get_vector()))
		self.mid = get_mid(self.start,self.end,tip_angle=PI/2)
		Polygon.__init__(self,self.start, self.mid, self.end, self.arc_center)

