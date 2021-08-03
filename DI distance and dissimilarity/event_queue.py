from collections import namedtuple
import functools

import plane_util as pu

Edge = namedtuple("Edge", ("a", "b"))

def birth(e, r):
    (a_x, a_y), (b_x, b_y) = e
    return max(b_x - a_x  - r, b_y - a_y - r)

def death(e, r):
    (a_x, a_y), (b_x, b_y) = e
    return min(b_x - a_x + r, b_y - a_y + r)

class Event:
    pass

def diag_edge_dist(edge):
    return pu.infty_metric(*edge)

class EntryEvent(Event):
    def __init__(self, edge, shift_to_check):
        self.edge = edge
        self.shift_to_check = shift_to_check

    def __str__(self):
        return "EntryEvent({},\n{})".format(self.edge, self.shift_to_check)

class ExitEvent(Event):
    def __init__(self, edge):
        self.edge = edge

    def __str__(self):
        return "ExitEvent({})".format(self.edge)

class QueueStateError(Exception):
    pass

def print_return(func):
    # @functools.wraps(func)
    # def _decorated(*args, **kwargs):
    #     res = func(*args, **kwargs)
    #     print(str(res))
    #     # print(func.__name__, str(args), str(kwargs), "-->", str(res))
    #     return res
    # return _decorated
    return func

class EventQueue:
    def __init__(self, A, B):
        _a_b_edges = [Edge(a, b) for a in A for b in B]
        self._edge_entries = Stack(sorted(_a_b_edges, key=lambda e: birth(e, 0)))
        self._edge_exits = Stack(sorted(_a_b_edges, key=lambda e: death(e, 0)))
        # create a stack `diag_edges` containing all non-skew diagonal edges
        self._diag_edges = Stack(sorted(
            ([Edge(a, pu.closest_diag_point(a)) for a in A] +
             [Edge(pu.closest_diag_point(b), b) for b in B]),
            key=lambda e: -diag_edge_dist(e)))
# pu.infty_metric(e.a, e.b)))

    def __bool__(self):
        return bool(self._edge_entries) and bool(self._edge_exits)

    def push(self, entry_event):
        # used to repeat an edge that caused the radius to decrease
        if not isinstance(entry_event, EntryEvent):
            raise TypeError("expected an instance of {}".format(EntryEvent))
        self._edge_entries.push(entry_event.edge)

    @print_return
    def next_event(self, radius):
        if self._diag_edges and diag_edge_dist(self._diag_edges.top()) >= radius:
            return ExitEvent(self._diag_edges.pop())
        if (death(self._edge_exits.top(), radius)
                <= birth(self._edge_entries.top(), radius)):
            return ExitEvent(self._edge_exits.pop())
        edge = self._edge_entries.pop()
        left_end = birth(edge, radius)
        while (self._edge_entries and
               # used to compare by  (lambda x: birth(x, 0)) -- this does not seem to be the issue
               birth(self._edge_entries.top(), radius) == left_end):
            # could also just throw the edge away, but this way lets us do a nice optimization.
            edge = max(edge, self._edge_entries.pop(), key=lambda e: death(e, radius))
            # self._edge_entries.pop()
        right_end = death(self._edge_exits.top(), radius)
        if self._edge_entries:
            right_end = min(right_end, birth(self._edge_entries.top(), radius))
        return EntryEvent(edge, (left_end + right_end) / 2)

    def next_diagonal_height(self):
        if not self._diag_edges:
            return 0
        else:
            return diag_edge_dist(self._diag_edges.top())

    def next_exit_shift(self, radius):
        return death(self._edge_exits.top(), radius)

class Stack:
    def __init__(self, items):
        self.items = list(items[::-1])

    def pop(self):
        return self.items.pop()

    def push(self, item):
        self.items.append(item)

    def top(self):
        return self.items[-1]

    def __len__(self):
        return len(self.items)

    def bool(self):
        return bool(len(self))

    def sort(self, key=lambda x: x):
        self.items.sort(key=lambda x: -key(x))

    def __repr__(self):
        return repr(self.items)
