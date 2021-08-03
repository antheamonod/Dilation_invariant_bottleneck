#!/usr/bin/env python3
import abc
import functools
import collections
import math
import random

# Point = collections.namedtuple("Point", ("x", "y"))

class NegativeCounterError(KeyError):
    pass

class SaneCounter(collections.Counter):
    def __setitem__(self, key, val):
        super().__setitem__(key, val)
        if val < 0:
            raise NegativeCounterError(
                "Cannot store negative count {} of {}".format(val, key))
        elif val is 0:
            del self[key]

class Point(tuple):
    def __new__(cls, *args):
        return tuple.__new__(cls, args)

def infty_metric(p1, p2):
    return max(abs(v1 - v2) for v1, v2 in zip(p1, p2))

def dist_from_diag(pt):
    # (infinity norm, distance from diagonal x=y)
    return abs(pt[0] - pt[1]) / 2

def closest_diag_point(pt):
    coord = (pt[0] + pt[1]) / 2
    return Point(coord, coord)

def ccw(a, b, c):
    # are the given three points in R^2 counterclockwise?
    # raise NotImplementedError()  # find a nice matrix library.  Also make this more general anyway
    num =  (b[0]*c[1] - b[1]*c[0]
            - a[0]*c[1] + a[1]*c[0]
            + a[0]*b[1] - a[1]*b[0])
    if num > 0: return 1
    if num < 0: return -1
    return 0

# 0.5 * / 1  -1 \
#       \ 1   1 /

def to_northeast(point):
    return (point[0] - point[1]) / 2

def to_northwest(point):
    return (point[0] + point[1]) / 2

def to_diagonal(point):
    # [ 0.5  -0.5 ] [x]
    # [ 0.5   0.5 ] [y]
    return Point(to_northwest(point), to_northeast(point))
def from_diagonal(point):
    return Point(point[0] + point[1],
                 -point[0] + point[1])

def sort_convex_vertices_ccw(*vertices):
    if len(vertices) < 3:
        raise ValueError("at least three vertices required")
    # Just pick a horiz. line through them and sort normally below, backward
    # above on the x coord
    y_divide = (vertices[0][1] + vertices[1][1]) / 2
    if y_divide == vertices[0][1]:
        y_divide = (vertices[1][1] + vertices[2][1]) / 2
        if y_divide == vertices[0][1]:
            raise ValueError("vertices are colinear: {}".format(vertices))
    upper = [vert for vert in vertices if vert[1] >= y_divide]
    lower = [vert for vert in vertices if vert[1] < y_divide]
    upper.sort(key=lambda pt: -pt[0])
    lower.sort(key=lambda pt: pt[0])
    return upper + lower

def convex_shape(*vertices, presorted=False):
    # A factory for convex shapes
    if not presorted:
        vertices = sort_convex_vertices_ccw(*vertices)
    if len(vertices) > 3:
        return BigConvexShape(*vertices)
    elif len(vertices) == 3:
        return Triangle(*vertices)
    else:
        raise ValueError("not enough vertices for a polygon: {}".format(vertices))

class BigConvexShape:
    def __init__(self, *vertices, presorted=True):
        if not presorted:
            vertices = sort_convex_vertices_ccw(*vertices)

        if len(vertices) <= 3:
            raise ValueError("Big shape only has {} vertices, {}".format(len(vertices), vertices))
        # split into two smaller shapes
        division = len(vertices) // 2
        self.side1 = convex_shape(*vertices[:division + 1])
        self.side2 = convex_shape(*vertices[division:], vertices[0])
        self.divider = (vertices[0], vertices[division])
        self.off_diagonal_pt = vertices[division - 1]

    def contains(self, point):
        orientation = ccw(*self.divider, point)
        if orientation is -1:
            return self.side1.contains(point)
        elif orientation is 1:
            return self.side2.contains(point)
        else:
            # point may lie on the diagonal
            assert orientation is 0
            return - (ccw(self.off_diagonal_pt, self.divider[0], point) *
                      ccw(self.off_diagonal_pt, self.divider[1], point))

            # for index in len(self.divider):
                
            # for coord in range(len(self.divider[0])):
            #     low_x, high_x = self.divider[0][coord], self.divider[1][coord]
            #     if low_x > high_x:
            #         low_x, high_x = high_x, low_x
            #     if not (low_x <= point[coord] <= high_x):
            #         return False
            # return True

    def intersects(self, foo):
        raise NotImplementedError()

class Triangle:
    def __init__(self, *vertices, presorted=True):
        if not presorted:
            vertices = sort_convex_vertices_ccw(vertices)
        self.vertices = vertices
        if len(vertices) != 3:
            raise ValueError("wrong number of vertices for triangle: {}"
                             .format(vertices))

    def contains(self, point):
        return min(ccw(point, self.vertices[i], self.vertices[i+1])
                   for i in range(-3, 0))


# class AbstractKDTree:
#     __metaclass__ = abc.ABCMeta

#     @abc.abstractmethod
#     def intersect(


# class 

# class KDTree(object):

#     def __init__(self, *points):

def quick_select(xs, k, key=lambda x: x):
    # k should be in range(0, len(xs)), i.e. so smallest is k=0
    if len(xs) <= k:
        raise ValueError("cannot select element {} from list of length {}: {}"
                         .format(k, len(xs), xs))
    if len(xs) == k - 1:
        return max(xs, key=key)
    elif k == 0:
        return min(xs, key=key)
    pivot = random.choice(xs)
    lower = [x for x in xs if key(x) < key(pivot)]
    if len(lower) <= k:
        upper = [x for x in xs if key(x) > key(pivot)]
        if len(xs) - len(upper) > k:
            return pivot
        return quick_select(upper, k + len(upper) - len(xs), key=key)
    else:
        return quick_select(lower, k, key=key)

class SimpleKDTree():
    """Simple in the sense that it chokes on repeated points."""

    def __init__(self, *points, _split_dim=0):
        if points:
            self._empty = False
            self.split_dim = _split_dim
            self.point = quick_select(points, len(points) // 2,
                                      key=lambda pt: pt[self.split_dim])
            self.count = len(points)
            self.deleted = False
            sides = [[], []]
            for pt in points:
                if pt is not self.point:
                    sides[self._upward(pt)].append(pt)
            self.children = [None if not side else
                             SimpleKDTree(*side, _split_dim=(self.split_dim + 1) % len(self.point))
                             for side in sides]
        else:
            self._empty = True
            self.count = 0

    def _upward(self, point):
        # which side is the point on?  Up (1) or down (0)?
        return point[self.split_dim] > self.point[self.split_dim]

    def __len__(self):
        return self.count

    # def _leftward(self, pt):
    #     # should the point go to the left subtree?  Else the right.
    #     return pt[self.split_dim] <= self.point[self.split_dim]

    def _delete(self, point, yell_if_absent):
        # return True if a point is successfully deleted
        def yell():
            if yell_if_absent:
                raise KeyError("Point {} not found for deletion".format(point))
            else:
                return False

        if self._empty:
            return yell()
        elif point == self.point:
            if self.deleted:
                return yell()
            else:
                self.deleted = True
                self.count -= 1
                return True
        elif self.count is 0:
            # This should only happen to the root
            return yell()
        else:
            side = self._upward(point)
            child = self.children[side]
            if child is None:
                yell()
                return False
            result = child._delete(point, yell_if_absent)
            if result:
                self.count -= 1
                if child.count is 0:
                    self.children[side] = None
            # no way we need to yell here
            return result

    def _collect_for_iter(self, items):
        if not self.deleted:
            items.append(self.point)
        for child in self.children:
            if child is not None:
                child._collect_for_iter(items)

    def __iter__(self):
        if self._empty:
            return
        items = []
        self._collect_for_iter(items)
        return iter(items)

    def __contains__(self, point):
        return self.neighbor(point, 0, closed=True)

    def delete(self, point):
        self._delete(point, yell_if_absent=True)

    def neighbors(self, point, radius, closed=True, want=-1, found=None):
        if found is None:
            found = []
        if len(found) == want or self._empty:
            return found
        if (not self.deleted and closed_less_than(
                infty_metric(point, self.point), radius, closed)):
            found.append(self.point)
        side = self._upward(point)
        if self.children[side] is not None:
            self.children[side].neighbors(point, radius, closed, want, found)
        dist_from_split = abs(point[self.split_dim] - self.point[self.split_dim])
        if (self.children[not side] is not None
            # and closed_less_than(dist_from_split, radius, closed)
            # Check out this amazing premature optimization!
            and (dist_from_split < radius or
                 closed and dist_from_split == radius and
                 side is not self._upward(self.point))):
            self.children[not side].neighbor(point, radius, closed, want, found)
        return found
                 

    def neighbor(self, point, radius, closed=True):
        if self._empty:
            return None
        # return a neighbor of the point within the radius, if possible
        dist = infty_metric(point, self.point)
        if not self.deleted:
            if closed_less_than(dist, radius, closed):
                return self.point
        side = self._upward(point)
        result = None
        if self.children[side] is not None:
            result = self.children[side].neighbor(point, radius, closed=closed)
        if result is not None:
            return result
        dist_from_split = abs(point[self.split_dim] - self.point[self.split_dim])
        if (self.children[not side] is not None
            and (dist_from_split < radius or
                 closed and dist_from_split == radius
                 and side is not self._upward(self.point))):
            return self.children[not side].neighbor(point, radius, closed)


class EfratNeighborStructure(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def delete(self, point, mult=1):
        raise NotImplementedError()

    @abc.abstractmethod
    def neighbor(self, point):
        raise NotImplementedError()

    @abc.abstractmethod
    def count(self, point):
        raise NotImplementedError()


class MultiEfratKDTree(EfratNeighborStructure):
    kd_tree_cls = SimpleKDTree
    def __init__(self, *ctr_args, closed=False, radius=math.inf):#, **ctr_kwargs):
        self.radius = radius
        self.closed = closed  # open or closed balls
        self.counter = self._ctr_from_args(*ctr_args)#  SaneCounter(ctr_args)#, **ctr_kwargs)
        self.tree = self.kd_tree_cls(*self.counter.keys())

    def neighbor(self, node):
        return self.tree.neighbor(node, self.radius, closed=self.closed)

    def _ctr_from_args(self, *args):
        if len(args) is 1 and isinstance(args[0], dict):
            return SaneCounter(*args)
        else:
            return SaneCounter(args)

    def delete(self, point, mult=1):
        if mult < -1:
            raise ValueError("Invalid multiplicity {}: must be -1 or greater")
        if mult is -1:
            del self.counter[point]
            self._delete_from_tree(point)
            return
        try:
            self.counter[point] -= mult
        except NegativeCounterError:
            assert mult > 0
            self.counter[point] += mult  # might as well fix it
            raise KeyError("Cannot remove {} of {}. Only {} remain"
                           .format(mult, point, self.counter[point]))
        if self.counter[point] is 0:
            self._delete_from_tree(point)
        #     del self.counter[point]
        # elif self.counter[point] < 0:
        #     assert mult > 0
        #     self.counter[point] += mult
        #     raise KeyError("Cannot remove {} of {}. Only {} remain"
        #                    .format(mult, point, self.counter[point]))

    def count(self, point):
        return self.counter[point]

    def _delete_from_tree(self, point):
        self.tree.delete(point)

def closed_less_than(a, b, closed=True):
    return a <= b if closed else a < b

class EfratTreeWithDiagonal(MultiEfratKDTree):
    def __init__(self, *ctr_args, diag_key=None, other_diag=None, **kwargs):
        counter = self._ctr_from_args(*ctr_args)
        diag_count = counter.pop(diag_key) if diag_key in counter else 0
        super().__init__(*counter.elements(), **kwargs)
        self.near_diagonal = set(pt for pt in self.counter
                                 if closed_less_than(dist_from_diag(pt),
                                                     self.radius,
                                                     closed=self.closed))
        self.diag_key = diag_key
        self.other_diag = other_diag
        if diag_count > 0:
            self.counter[self.diag_key] = diag_count
            self.near_diagonal.add(diag_key)

    def neighbor(self, point):
        if (self.counter[self.diag_key] > 0
            and (point == self.other_diag
                 or closed_less_than(dist_from_diag(point),
                                     self.radius, self.closed))):
            return self.diag_key
        elif point == self.other_diag:
            # return an arbitrary point near the diagonal
            if self.near_diagonal:
                # get an arbitrary point from self.near_diagonal
                val = self.near_diagonal.pop()
                self.near_diagonal.add(val)  # a bit daft
                return val
            else:
                return None
        else:
            return super().neighbor(point)

    def delete(self, point, mult=1):
        super().delete(point, mult=mult)
        if self.counter[point] is 0 and point in self.near_diagonal:
            self.near_diagonal.remove(point)

    def _delete_from_tree(self, point):
        if point != self.diag_key:
            super()._delete_from_tree(point)

    def __repr__(self):
        return repr(self.counter).replace(type(self.counter).__name__, "Efrat")
