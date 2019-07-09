#-- my_code_hw03.py
#-- Assignment 03 GEO1015/2018
#-- [YOUR NAME] Yifang Zhao
#-- [YOUR STUDENT NUMBER] 4798899
#-- [YOUR NAME] Jinglan Li
#-- [YOUR STUDENT NUMBER] 4781937

import sys
import math
import numpy
import rasterio
from rasterio import features



def output_viewshed(d, viewpoints, maxdistance, output_file):
    """
    !!! TO BE COMPLETED !!!
     
    Function that writes the output raster
     
    Input:
        d:            the input datasets (rasterio format)  
        viewpoints:   a list of the viewpoints (x, y, height)
        maxdistance:  max distance one can see
        output_file:  path of the file to write as output
        
    Output:
        none (but output GeoTIFF file written to 'output-file')
    """  
    # basic information about the original grid
    npi = d.read(1)
    x_cellsize = d.transform[0]
    y_cellsize = d.transform[4]
    allvs = []
        
    # 1. compute the viewshed for each viewpoint
    for vp in viewpoints:
        vs = numpy.zeros(d.shape, dtype=numpy.int8)
        vrow, vcol = d.index(vp[0], vp[1])
        vs[vrow , vcol] = 2
        pt_vp = Point(vp[0], vp[1])
        c_vp = Circle(pt_vp, maxdistance)
        edge = []# centers of the pixels intersected with the c_vp

        # 1.1. find the edge pixels and mark the pixels outside as 3
        for i in range(d.height):
            for j in range(d.width):
                # viewpoints (2) not considered
                if vs[i, j] == 0:
                    coor = d.xy(i, j)
                    # pixel vertices
                    pt_ll = Point(coor[0] - x_cellsize / 2, coor[1] + y_cellsize / 2)
                    pt_ur = Point(coor[0] + x_cellsize / 2, coor[1] - y_cellsize / 2)
                    pt_lr = Point(coor[0] + x_cellsize / 2, coor[1] + y_cellsize / 2)
                    pt_ul = Point(coor[0] - x_cellsize / 2, coor[1] - y_cellsize / 2)
                    # pixel rectangle
                    rect = Rectangle(pt_ll, pt_ur)
                    
                    # if the pixel is outside: 3
                    if not rect.intersects(c_vp):
                        vs[i, j] = 3
                    # if the pixel are inside: pass
                    elif pt_ul.intersects(c_vp) and pt_ll.intersects(c_vp) and pt_lr.intersects(c_vp) and pt_ur.intersects(c_vp):
                        if i == 0 or i == d.height - 1 or j == 0 or j == d.width - 1:
                            edge.append((coor[0], coor[1]))
                    # if the pixel intersects the circle
                    else:
                        edge.append((coor[0], coor[1]))

        # 1.2. compute the viewshed using tangents for each LoS
        for p in edge:
            tan_cur = -1000000
            pixels_LoS = {}
            
            # 1.2.1. find the pixels that intersect the line of sight, 'vp-p'. In fact it is a line segment.
            for i in range(d.height):
                for j in range(d.width):
                    # viewpoints (2) and pixels outside (3) not considered
                    if vs[i, j] != 2 and vs[i, j] != 3:
                        coor = d.xy(i, j)
                        #pixel vertices
                        pt_ll = (coor[0] - x_cellsize / 2, coor[1] + y_cellsize / 2)
                        pt_ur = (coor[0] + x_cellsize / 2, coor[1] - y_cellsize / 2)
                        pt_lr = (coor[0] + x_cellsize / 2, coor[1] + y_cellsize / 2)
                        pt_ul = (coor[0] - x_cellsize / 2, coor[1] - y_cellsize / 2)
                        pt_center = (coor[0], coor[1])

                        # pixel not intersecting the LoS
                        if (dy(pt_ll, vp, p)>=0 and dy(pt_lr, vp, p)>=0 and dy(pt_ul, vp, p)>=0 and dy(pt_ur, vp, p)>=0) \
                           or (dy(pt_ll, vp, p)<=0 and dy(pt_lr, vp, p)<=0 and dy(pt_ul, vp, p)<=0 and dy(pt_ur, vp, p)<=0) \
                           or (pt_center[0]<vp[0] and pt_center[0]<p[0]) \
                           or (pt_center[0]>vp[0] and pt_center[0]>p[0]) \
                           or (pt_center[1]<vp[1] and pt_center[1]<p[1]) \
                           or (pt_center[1]>vp[1] and pt_center[1]>p[1]):
                            pass
                        # pixel intersecting the LoS
                        else:
                            v_d = distance2(pt_center, vp, p)
                            dis2vp = math.sqrt(distance1(vp, pt_center)**2 - v_d**2)
                            pixels_LoS.update({dis2vp:(i, j)})
                            
            distances = list(pixels_LoS)
            distances.sort()
            
            # 1.2.2. compute the viewshed
            for dis in distances:
                row = pixels_LoS[dis][0]
                col = pixels_LoS[dis][1]
                tan = (npi[row, col] - npi[vrow, vcol] - vp[2]) / dis
                # visible
                if tan > tan_cur:
                    vs[row, col] = 1
                    tan_cur = tan
                    
        allvs.append(vs)

    # 2. map algebra
    npvs = vs_union(allvs)
        
    #-- write this to disk
    with rasterio.open(output_file, 'w', 
                       driver='GTiff', 
                       height=npi.shape[0],
                       width=npi.shape[1], 
                       count=1, 
                       dtype=rasterio.uint8,
                       crs=d.crs, 
                       transform=d.transform) as dst:
        dst.write(npvs.astype(rasterio.uint8), 1)

    print("Viewshed file written to '%s'" % output_file)

def dy(p, p1, p2):
    """
        Compute the vertical distance between point p and line p1p2 in 2D space.
        Line p1p2: y = kx + b
    """
    k = (p2[1] - p1[1]) / (p2[0] - p1[0])
    b = p2[1] - k * p2[0]
    return p[1] - (k * p[0] + b)

def distance1(p1, p2):
    """
        Compute the distance between two points in 2D space.
    """
    return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)

def distance2(p, p1, p2):
    """
        Compute the distance between point p and line p1p2 in 2D space.
        Line p1p2: y = kx + b
    """
    k = (p2[1] - p1[1]) / (p2[0] - p1[0])
    b = p2[1] - k * p2[0]
    return abs((k * p[0] - p[1] + b) / math.sqrt(k**2 + 1))

def vs_union(vs_list):
    vs = numpy.zeros(vs_list[0].shape, dtype=numpy.int8)
    for i in range(vs_list[0].shape[0]):
        for j in range(vs_list[0].shape[1]):
            vals = []
            for vs in vs_list:
                vals.append(vs[i, j])
            if 2 in vals:
                vs[i, j] = 2
            elif 1 in vals:
                vs[i, j] = 1
            elif 0 in vals:
                vs[i, j] = 0
            else:
                vs[i, j] = 3
    return vs

class Point(object):

    def __init__(self, x, y):
        """Constructor. 
        Takes the x and y coordinates to define the Point instance.
        """
        self.x = float(x)
        self.y = float(y)
        
    def intersects(self, other):
        """Checks whether other shape has any interaction with
        interior or boundary of self shape. Uses type based dispatch.
        
        other - Point, Circle or Rectangle
        
        returns - True / False
        """
        if isinstance(other, Point) and self.x == other.x and self.y == other.y:
            return True
        elif isinstance(other, Circle) and self.distance(other.center) <= other.radius:
            return True
        elif isinstance(other, Rectangle) and other.ll.x <= self.x <= other.ur.x and other.ll.y <= self.y <= other.ur.y:
            return True
        else:
            return False

    def distance(self, other):
        """Returns cartesian distance between self and other Point
        """
        return math.sqrt((self.x - other.x)**2 + (self.y - other.y)**2)
    
class Circle(object):

    def __init__(self, center, radius):
        """Constructor. 
        Takes the center point and radius defining the Circle.
        """
        assert radius > 0
        assert isinstance(center, Point)
        self.center = center
        self.radius = radius
        
    def intersects(self, other):
        """Checks whether other shape has any interaction with
        interior or boundary of self shape. Uses type based dispatch.
        
        other - Point, Circle or Rectangle
        
        Returns - True / False
        """
        if isinstance(other, Point):
            return other.intersects(self)
        elif isinstance(other, Circle) and self.center.distance(other.center) <= self.radius + other.radius:
            return True
        elif isinstance(other, Rectangle):
            rect1 = Rectangle(Point(other.ll.x, other.ll.y - self.radius), Point(other.ur.x, other.ll.y)) #lower rectangle
            rect2 = Rectangle(Point(other.ll.x, other.ur.y), Point(other.ur.x, other.ur.y + self.radius)) #upper rectangle
            rect3 = Rectangle(Point(other.ll.x - self.radius, other.ll.y), Point(other.ll.x, other.ur.y)) #left  rectangle
            rect4 = Rectangle(Point(other.ur.x, other.ll.y), Point(other.ur.x + self.radius, other.ur.y)) #right rectangle
            circle1 = Circle(other.ll, self.radius)                      #lower-left  circle
            circle2 = Circle(Point(other.ur.x, other.ll.y), self.radius) #lower-right circle
            circle3 = Circle(other.ur, self.radius)                      #upper-right circle
            circle4 = Circle(Point(other.ll.x, other.ur.y), self.radius) #upper-left  circle
            if self.center.intersects(other) or \
               self.center.intersects(rect1) or self.center.intersects(rect2) or self.center.intersects(rect3) or self.center.intersects(rect4) or \
               self.center.intersects(circle1) or self.center.intersects(circle2) or self.center.intersects(circle3) or self.center.intersects(circle4):
                return True
            else:
                return False
        else:
            return False

class Rectangle(object):

    def __init__(self, pt_ll, pt_ur):
        """Constructor. 
        Takes the lower left and upper right point defining the Rectangle.
        """
        assert isinstance(pt_ll, Point)
        assert isinstance(pt_ur, Point)
        self.ll = pt_ll
        self.ur = pt_ur

    def intersects(self, other):
        """Checks whether other shape has any interaction with
        interior or boundary of self shape. Uses type based dispatch.
        
        other - Point, Circle or Rectangle
        
        Returns - True / False
        """
        if isinstance(other, Point):
            return other.intersects(self)
        elif isinstance(other, Circle):
            return other.intersects(self)
        elif isinstance(other, Rectangle) and \
             math.fabs(self.center().x - other.center().x) <= 0.5 * (self.width() + other.width()) and \
             math.fabs(self.center().y - other.center().y) <= 0.5 * (self.height() + other.height()):
            return True
        else:
            return False
