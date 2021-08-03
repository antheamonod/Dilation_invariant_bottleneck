from collections import namedtuple

# import dionysus
from bipartite_matching import GeometricBipartiteMatching
import plane_util as pu
import event_queue
from event_queue import Edge, birth, death, Stack

import json

epsilon = 0.00000001

# What if we stop computing the value of the matching and just start ticking 
# the radius down to where the next edge will die?

def default_fudge(r):
    return r * (1 - epsilon)

def upper_bound_on_radius(A, B):
    return max(pu.dist_from_diag(x) for x in A + B)

def shifted_bottleneck_distance(A, B, fudge=default_fudge, analysis=True):
    """Compute the shifted bottleneck distance between two diagrams, A and B (multisets)"""
    A = pu.SaneCounter(A)
    B = pu.SaneCounter(B)
    if not A and not B:
        return 0
    radius = fudge(upper_bound_on_radius(A, B))
    events = event_queue.EventQueue(A, B)
    matching = GeometricBipartiteMatching(A, B)
    print(matching)
    # these counters are for performance monitoring only - they don't affect the logic
    ctr, R_ctr, L_ctr, fail_ctr, win_ctr = 0, 0, 0, 0, 0
    while events and radius > epsilon:
        ctr += 1
        event = events.next_event(radius)
        if isinstance(event, event_queue.ExitEvent):
            R_ctr += 1
            matching.remove_all(event.edge)
        else:
            L_ctr += 1
            if birth(event.edge, radius) >= death(event.edge, radius):
                win_ctr += 1
                continue  # relies on ties being broken with the highest-radius edge
            # assert not matching.diagonal_perfect()
            if matching.diagonal_perfect():
                fail_ctr += 1
                radius = fudge(max(
                    events.next_diagonal_height(),
                    radius - (events.next_exit_shift(radius)
                              - birth(event.edge, radius)) / 2))
                events.push(event)
                continue
            matching.maximize_matching(
                shift=event.shift_to_check,
                radius=radius)
            if matching.diagonal_perfect():
                radius = fudge(matching.value())
                events.push(event)
    if analysis:
        print(len(A) + len(B), "total", ctr, "R", R_ctr, "L", L_ctr, "fail", fail_ctr, "win", win_ctr)
    return radius

def cyq_test(A, B, fudge=default_fudge, analysis=False):
    """I write this part according to Alg 2. in cccg paper////date:2021.4.28"""
    A = pu.SaneCounter(A)
    B = pu.SaneCounter(B)
    if not A and not B:
        return 0
    radius = fudge(upper_bound_on_radius(A, B))
    events = event_queue.EventQueue(A, B)
    matching = GeometricBipartiteMatching(A, B)
    # these counters are for performance monitoring only - they don't affect the logic
    R_ctr, L_ctr= 0, 0
    while events and radius > epsilon:
        event = events.next_event(radius)
        if isinstance(event, event_queue.ExitEvent):  ### R event
            ### remove e.edge from matching
            R_ctr += 1
            matching.remove_all(event.edge)
        else:   
            ### L event
            L_ctr += 1
            matching.maximize_matching(shift = event.shift_to_check, radius = radius)
            if matching.diagonal_perfect():
                events._edge_entries.push(event.edge)
                radius = max(events.next_diagonal_height(),
                             radius - (events.next_exit_shift(radius) 
                             - event_queue.birth(event.edge, radius)) / 2)
 
    if analysis:
        print(len(A) + len(B), "R", R_ctr, "L", L_ctr)
    return radius


def other_shifted_bottleneck_distance(A, B, fudge=default_fudge, analysis=False):
    """Compute the shifted bottleneck distance between two diagrams, A and B (multisets)"""
    A = pu.SaneCounter(A)
    B = pu.SaneCounter(B)
    if not A and not B:
        return 0
    radius = fudge(upper_bound_on_radius(A, B))
    events = event_queue.EventQueue(A, B)
    matching = GeometricBipartiteMatching(A, B)
    # these counters are for performance monitoring only - they don't affect the logic
    ctr, R_ctr, L_ctr, fail_ctr, win_ctr = 0, 0, 0, 0, 0
    while events and radius > epsilon:
        ctr += 1
        event = events.next_event(radius)
        if isinstance(event, event_queue.ExitEvent):
            R_ctr += 1
            matching.remove_all(event.edge)
        else:
            L_ctr += 1
            if birth(event.edge, radius) >= death(event.edge, radius):
                win_ctr += 1
                continue  # relies on ties being broken with the highest-radius edge
            # assert not matching.diagonal_perfect()
            if matching.diagonal_perfect():
                fail_ctr += 1
                radius = fudge(max(
                    events.next_diagonal_height(),
                    radius - (events.next_exit_shift(radius)
                              - birth(event.edge, radius)) / 2))
                events.push(event)
                continue
            matching.maximize_matching(
                shift=event.shift_to_check,
                radius=radius)
            if matching.diagonal_perfect():
                # radius = fudge(matching.value())
                events.push(event)
    if analysis:
        print("other:", len(A) + len(B), "total", ctr, "R", R_ctr, "L", L_ctr, "fail", fail_ctr, "win", win_ctr)
    return radius

def other_crappy_normal_distance(A, B):
    if not A and not B:
        return 0
    radius = max(pu.dist_from_diag(pt) for pt in A + B)
    edges = Stack(([Edge(a, b) for a in A for b in B]
                   + [Edge(a, pu.closest_diag_point(a)) for a in A]
                   + [Edge(pu.closest_diag_point(b), b) for b in B]))
    edges.sort(key=lambda e: -pu.infty_metric(*e))
    matching = GeometricBipartiteMatching(A, B)
    while True:
        while edges and pu.infty_metric(*edges.top()) >= radius:
            matching.remove_all(edges.pop())
        matching.maximize_matching(radius=radius, shift=0, closed=False)
        print("c")
        if not matching.diagonal_perfect():
            return radius
        new_radius = matching.value_at_shift(0)
        assert new_radius < radius
        radius = new_radius

# def botdist(A, B):
#     return dionysus.bottleneck_distance(dionysus.Diagram(A), dionysus.Diagram(B))
# Their underlying metric is Euclidean.  (Ours is L_\infty)

def bin_search(size, direction):
    left, right = 0, size
    while left < right:
        mid = (right + left) // 2
        dir_ = direction(mid)
        if dir_ < 0:
            right = mid
        elif dir_ > 0:
            left = mid + 1
        else:
            return mid
    return left

def simple_botdist(A, B):
    if not A and not B:
        return 0
    edges = sorted([Edge(a, b) for a in A for b in B]
                   + [Edge(a, pu.closest_diag_point(a)) for a in A]
                   + [Edge(pu.closest_diag_point(b), b) for b in B],
                   key=(lambda e: pu.infty_metric(*e)))
    matching = GeometricBipartiteMatching(A, B)
    prev_ix = None
    def edge_len(ix):
        return pu.infty_metric(*edges[ix])

    def direction(ix):
        nonlocal prev_ix, matching
        radius = edge_len(ix)
        # if radius < prev_radius:
        #     matching.maximize
        if prev_ix is not None and prev_ix > ix:
            # we allowed longer edges last time, so remove those
            for edge in edges[ix + 1 : prev_ix + 1]:
                matching.remove_all(edge)
        matching.maximize_matching(radius, 0, closed=True)
        prev_ix = ix
        return -1 if matching.diagonal_perfect() else 1
    return edge_len(bin_search(len(edges), direction))

# def bin_search(xs, target, key=lambda x: x):
#     left, right = 0, len(xs)
#     target_val = key(target)
#     while left < right:
#         mid = (right + left) // 2
#         mid_val = key(mid)
#         if mid_val == target_val:
#             return mid
#         elif mid_val > target_val:
#             left = mid + 1
#         else:
#             right = mid
#     return left

def instance_from_file(file_):

    instance = json.load(file_)
    A = [pu.Point(*pt) for pt in instance['A']]
    B = [pu.Point(*pt) for pt in instance['B']]
    return A, B

if __name__ == "__main__":
    import sys
    A, B = instance_from_file(sys.stdin)
    print(shifted_bottleneck_distance(A, B))
