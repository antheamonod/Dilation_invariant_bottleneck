import abc
import math
import collections
# When nodes have multiplicity, the right notion of "disjoint" for paths is
# the total use of each edge does not exceed the edge's capacity.  Hence
# counters.
import plane_util as pu
from plane_util import SaneCounter as Counter

import logging
logger = logging.getLogger(__name__)
import sys

hdlr = logging.StreamHandler()
hdlr.setLevel(logging.WARNING)
logger.addHandler(hdlr)
logger.setLevel(logging.WARNING)

class Matching:

    def __init__(self):
        self.ctr = Counter()
        self.A_to_B = collections.defaultdict(set)
        self.B_to_A = collections.defaultdict(set)

    def remove_edge(self, a, b):
        count = self.ctr[(a, b)]
        del self.ctr[(a, b)]
        self.A_to_B[a].discard(b)
        self.B_to_A[b].discard(a)
        return count

    def degree(self, pt, in_A=True):
        if in_A:
            return sum(self.ctr[(pt, b)] for b in self.A_to_B[pt])
        else:
            return sum(self.ctr[(a, pt)] for a in self.B_to_A[pt])

    def has_edge(self, a, b):
        # return b in self.A_to_B[a]
        return bool(self.count_edge(a, b))

    def count_edge(self, a, b):
        return self.ctr[(a, b)]

    def __len__(self):
        # Number of edges in the matching
        return sum(self.ctr.values())

    def _augment_edge(self, a, b, mult):
        try:
            self.ctr[(a, b)] += mult
        except pu.NegativeCounterError:
            #logger.error(f"Augmenting edge {(a, b)} by {mult} resulted in negative multiplicity {self.ctr[(a, b)]}!")
            raise
        if self.ctr[(a, b)] == 0:
            self.remove_edge(a, b)
            # self.A_to_B[a].remove(b)
            # self.B_to_A[b].remove(a)
            # del self.ctr[(a, b)]
        else:
            self.A_to_B[a].add(b)
            self.B_to_A[b].add(a)

    def augment_path(self, *nodes, mult=1):
        even = True
        for cur, nxt in zip(nodes, nodes[1:]):
            if even:
                self._augment_edge(cur, nxt, mult)
            else:
                self._augment_edge(nxt, cur, -mult)
            even ^= True

    def __iter__(self):
        yield from self.ctr.elements()

class GeometricBipartiteMatching:
    efrat_cls = pu.EfratTreeWithDiagonal

    A_diag = "A_diag"
    B_diag = "B_diag"

    def __init__(self, A, B):
        # A and B are collections.Counters full of points
        #logger.info("New %s instance with A=%s and B=%s", type(self), A, B)
        self.A = Counter(A)
        self.B = Counter(B)
        self.A[self.A_diag] = sum(self.B.values())
        self.B[self.B_diag] = sum(self.A.values()) - self.A[self.A_diag]
        self.A_exposed = Counter(self.A)
        self.B_exposed = Counter(self.B)
        self.matching = Matching()
        self._prev_path = None

    # def add_edge(source, dest, mult=1):
    #     pass  # edges will be totally implicit anyway

    def remove_all(self, edge):
        #logger.info("remove_all called with edge=%s", edge)
        # Remove all instances of the given edge
        a, b = self._fix_diag_edge(edge)
        assert a in self.A and b in self.B
        count = self.matching.remove_edge(a, b)
        if count > 0:
            self.A_exposed[a] += count
            self.B_exposed[b] += count

    def _fix_diag_edge(self, edge):
        # Translate a diagonal edge to our representation, e.g.
        # turn ((1.5, 1.5), (1, 2)) into (self.A_diag, (1, 2)).
        def float_close(x, y):
            return abs(x - y) < 0.0000000001
        a, b = edge
        if a != self.A_diag and float_close(*a):
            a = self.A_diag
        if b != self.B_diag and float_close(*b):
            b = self.B_diag
        return (a, b)

    def has_edge(self, a, b):
        return self.matching.has_edge(a, b)

    def maximize_matching(self, radius, shift, closed=False):
        #logger.info("maximize_matching called with radius=%s, shift=%s, closed=%s", radius, shift, closed)
        while True:
            #logger.debug("loop asdf")
            image, inverse_image = self._make_translation(shift)
            # reverse_image = self._make_translation(-shift)
            # build a layer subgraph with pu.EfratTreeWithDiagonal
            layer_subgraph = self.build_layer_subgraph(radius, shift, closed=closed)
            #logger.debug("A_exposed: %s", self.A_exposed)
            #logger.debug("B_exposed: %s", self.B_exposed)
            #logger.debug("Matching: %s", self.matching.ctr)
            #logger.debug("Layer subgraph: %s", layer_subgraph)
            if not layer_subgraph[-1]:
                #logger.debug("Fount nothing!")
                # no augmenting paths were found
                return
            else:
                # do a binary search backward in the layer subgraph
                # From here on, I got lazy and handled multiplicity badly
                even_layers = [layer for i, layer in enumerate(layer_subgraph)
                               if i % 2 == 0]
                odd_layers = [layer for i, layer in enumerate(layer_subgraph)
                              if i % 2 == 1]
                def image_of_ctr(ctr):
                    return Counter({image(key): ctr[key] for key in ctr})
                even_efrats = [self.efrat_cls(image_of_ctr(self.A_exposed),
                                              radius=radius,
                                              diag_key=self.A_diag, other_diag=self.B_diag,
                                              closed=closed)]
                even_efrats += [self.efrat_cls(image_of_ctr({x: self.A[x] for x in layer if x in self.A}),
                                               radius=radius, diag_key=self.A_diag,
                                               other_diag=self.B_diag,closed=closed)
                                for layer in even_layers[1:]]
                while True:
                    #logger.debug("this loop awefoih")
                    (aug_path, count) = self._dfs_layers(even_efrats, odd_layers, image, inverse_image)
                    #logger.debug("awefoih %s", aug_path)
                    if not aug_path:
                        break
                    else:
                        assert count > 0
                        self.matching.augment_path(*aug_path, mult=count)
                        # Do some boring maintenance
                        assert aug_path[0] in self.A_exposed
                        assert aug_path[-1] in self.B_exposed
                        self.A_exposed[aug_path[0]] -= count
                        # if self.A_exposed[aug_path[0]] is 0:
                        #     del self.A_exposed[aug_path[0]]
                        self.B_exposed[aug_path[-1]] -= count
                        # if self.B_exposed[aug_path[-1]] is 0:
                        #     del self.B_exposed[aug_path[-1]]

    def _make_translation(self, shift):
        # The two errors we're avoiding here are:
        #  - isinstance(inverse_image(image(pt)), pu.Point) but isinstance(pt, tuple)
        #  - floating point error leading to infinite loops
        inverse_dict = {}
        def translate(a):
            if a is self.A_diag or a is self.B_diag:
                result = a
            else:
                result = pu.Point(*(coord + shift for coord in a))
            inverse_dict[result] = a
            return result
        def undo_translate(pt):
            if pt in inverse_dict:
                return inverse_dict[pt]
            else:
                # maybe it's a tuple
                return inverse_dict[pu.Point(*pt)]
        return translate, undo_translate

    def build_layer_subgraph(self, radius, shift, closed=False):
        image, _ = self._make_translation(shift)

        layers = []
        efrat = self.efrat_cls(*self.B.elements(), diag_key=self.B_diag,
                               other_diag=self.A_diag, radius=radius, closed=closed)
        A_reached = set()  # reached doesn't need to be a Counter, since we don't
        # want cycles.  (We're looking for the *shortest* augmenting paths,
        # which will be simple.)
        # (The Efrat thing keeps track of the reached elements of B for us)
        layers.append(Counter(self.A_exposed))  # TODO -- replace sets with Counters
        A_reached.update(layers[0])
        B_exposed_set = set(self.B_exposed)
        while True:
            #logger.debug("this loop too!")
            # build odd layer
            odd_layer = Counter()
            for point in layers[-1]:
                neighbor = efrat.neighbor(image(point))
                while neighbor is not None:
                    odd_layer[neighbor] += 1
                    efrat.delete(neighbor, mult=1)
                    neighbor = efrat.neighbor(image(point))
            layers.append(odd_layer)
            # A_reached.update(odd_layer)
            if not odd_layer or B_exposed_set.intersection(set(odd_layer)):
                for point in set(odd_layer):
                    if point in B_exposed_set:
                        odd_layer[point] = min(odd_layer[point], self.B_exposed[point])
                        assert odd_layer[point] > 0
                    else:
                        # Not a valid endpoint at all
                        del odd_layer[point]
                return layers
            # now build the even layer
            even_layer = Counter()
            for point in layers[-1]:
                new_stuff = self.matching.B_to_A[point].difference(A_reached)
                for a_point in new_stuff:
                    assert self.matching.has_edge(a_point, point)
                    even_layer[a_point] += self.matching.count_edge(a_point, point)
                A_reached.update(new_stuff)
            layers.append(even_layer)
        # return layers

    def _dfs_layers(self, even_efrats, odd_layers, image, inverse_image):
        # While you are at a node, search from that node.  If you can reach another node,
        # begin searching from that node.

        # If you can't reach anything from your current node, delete it and backtrack.

        # If you reach the end, determine the multiplicity of your path, and reduce the
        # multiplicity of each node by that much.
        # (TODO)
        # Let's un-split these things, like a monster
        layers = []
        assert len(even_efrats) == len(odd_layers)
        for even_layer, odd_layer in zip(even_efrats, odd_layers):
            layers.append(even_layer)
            layers.append(odd_layer)
        #logger.debug("Layers:  %s", layers)
        # for node in layers[-1]:
        layer = len(layers)
        path = [None] * len(layers)
        while layer > 0:
            next_layer = layers[layer - 1]
            assert path[layer - 1] is None
            if not next_layer:
                return (None, 0)
            if layer == len(layers):
                path[layer - 1] = next(iter(next_layer))
                layer -= 1
            elif layer % 2:
                assert path[layer] is not None
                # we're searching from odd to even, so next_layer is an
                # efrat structure
                assert isinstance(next_layer, pu.EfratNeighborStructure)
                query_result = next_layer.neighbor(path[layer])
                if query_result is not None:
                    query_result = inverse_image(query_result)
                path[layer - 1] = query_result
                # print(f"reached {path[layer -1]} from {path[layer]} via efrat (search key: {reverse_image(path[layer])})")
                # print(f"search key for (0, 0): {reverse_image((0, 0))}")
                if path[layer - 1] is None:
                    # delete the current node from the counter
                    del layers[layer][path[layer]]
                    path[layer] = None
                    layer += 1
                else:
                    layer -= 1
            else:
                # searching from even to odd
                assert path[layer] is not None
                assert isinstance(next_layer, Counter)
                next_nodes = self.matching.A_to_B[path[layer]].intersection(next_layer)
                if next_nodes:
                    path[layer-1] = next(iter(next_nodes))
                    assert self.matching.count_edge(path[layer], path[layer-1]) > 0
                    layer -= 1
                else:
                    # delete all occurrences of this useless node
                    layers[layer].delete(image(path[layer]), mult=-1)
                    path[layer] = None
                    layer += 1

        def remove(node, layer, mult):
            if isinstance(layer, pu.EfratNeighborStructure):
                layer.delete(image(node), mult=mult)
            else:
                layer[node] -= mult
                assert layer[node] >= 0
                if layer[node] is 0:
                    del layer[node]
        # TODO multiplicity as on blackboard
        # mult = self.A[path[0]] - self.matching.degree(path[0], in_A=True)
        # mult = self.B[path[-1]] - min(mult, self.matching.degree(path[-1], in_A=False))
        #logger.debug("Finding multiplicity of path %s", path)
        mult = layers[0].count(image(path[0]))
        #logger.debug("multiplicity of %s in %s: %s", path[0], layers[0], layers[0].count(image(path[0])))
        mult = min(mult, layers[-1][path[-1]])
        #logger.debug("multiplicities of endpoints: %s", mult)
        assert mult > 0
        for index in range(1, len(path) - 1, 2):
            # if index % 2 == 1:
                # print(path[index+1], path[index], self.matching.count_edge(path[index+1], path[index]))
            mult = min(mult,
                       self.matching.count_edge(path[index+1],
                                                path[index]))
        assert mult > 0  # again
        # mult = min(get_multiplicity(node, index)
        #            for index, node in enumerate(path))
        # assert mult > 0
        for node, layer in zip(path, layers):
            remove(node, layer, mult)
        #logger.debug("Augmenting: %s %s", path, mult)
        #logger.debug("Modified layer subgraph: %s", layers)
        return (path, mult)

    def diagonal_perfect(self):
        # "diagonal-perfect" is my word for a matching in which the degree of
        # each off-diagonal vertex is one.
        return not (len(self.A_exposed) + len(self.B_exposed))

    def edges(self, repeats=False):
        # Iterate through the edges, with each edge yielded as many times as it is present
        if repeats:
            yield from self.matching
        else:
            yield from self.matching.ctr

    def value(self, force=False):
        # The minimum, over all shifts, of the length of the longest edge in
        # `self.matching`
        if not self.diagonal_perfect() and not force:
            raise ValueError("Matching is not diagonal-perfect. "
                             "To compute value anyway, pass `force=True`")
        else:
            return self.matching_value(self.matching)

    def value_at_shift(self, shift):
        return self.matching_value_fixed_shift(self.matching, shift=shift)

    @classmethod
    def matching_value(cls, matching):
        # Okay.  So for every non-diagonal edge, the plot of value against shift looks like a vee, like a
        # shifted plot of `y=|x|`.  It's easy to find the upper hull of all these vees.
        vee_bottom = None
        longest_diag_edge = 0
        for a, b in matching.ctr:
            if a == cls.A_diag and b == cls.B_diag:
                continue
            elif a == cls.A_diag:
                longest_diag_edge = max(longest_diag_edge, pu.dist_from_diag(b))
            elif b == cls.B_diag:
                longest_diag_edge = max(longest_diag_edge, pu.dist_from_diag(a))
            elif vee_bottom is None:
                vee_bottom = edge_to_vee((a, b))
            else:
                other_vee = edge_to_vee((a, b))
                if (pu.ccw(other_vee,
                           vee_bottom,
                           (other_vee[0] + 1, other_vee[1] + 1)) == 1):
                    vee_bottom = intersect_diagonal_lines(vee_bottom, other_vee)
                if (pu.ccw(other_vee,
                           vee_bottom,
                           (other_vee[0] + 1, other_vee[1] - 1)) == 1):
                    vee_bottom = intersect_diagonal_lines(other_vee, vee_bottom)
        if vee_bottom is None:
            return longest_diag_edge
        else:
            return max(longest_diag_edge, vee_bottom[1])

    @classmethod
    def matching_value_fixed_shift(cls, matching, shift=0):
        def _dist(a, b, shift):
            if a == cls.A_diag and b == cls.B_diag:
                return 0
            if a == cls.A_diag:
                return pu.dist_from_diag(b)
            if b == cls.B_diag:
                return pu.dist_from_diag(a)
            return pu.infty_metric(pu.Point(*(coord + shift for coord in a)), b)
        return max(_dist(a, b, shift) for (a, b) in matching.ctr)

def intersect_diagonal_lines(downslope, upslope):
    # downslope is a point on the line with slope -1
    # upslope is a point on the line with slope 1

    # compute the x intercepts
    down_intercept = downslope[0] + downslope[1]
    up_intercept = upslope[0] - upslope[1]
    # now look halfway between them for the intersection
    diff = (up_intercept - down_intercept) / 2
    return pu.Point(down_intercept + diff, -diff)

def edge_to_vee(edge):
    # Return a point (shift, radius), where
    # `birth(edge, radius) == death(edge, radius) == shift`
    a, b = edge
    # Since we're shifting `a` by some multiple of (1, 1), we find the
    # intersection of the line with that slope through `a` with the
    # perpendicular line through b.
    a_shifted = intersect_diagonal_lines(b, a)
    shift = a_shifted[0] - a[0]
    dist = pu.infty_metric(a_shifted, b)
    return (shift, dist)

# def edge_to_vee(edge):
#     a, b = edge
#     height_diff = pu.dist_from_diag(b) - pu.dist_from_diag(a)
#     closest_y = b[1] - height_diff
#     shift = closest_y - a[1]
#     return (shift, abs(height_diff))
