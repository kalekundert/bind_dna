#!/usr/bin/env python3

from stepwise.reaction.mix import *

from heapq import heapify, heappush, heappop

class PartialMix(frozenset):
    """
    This class is used internally by `iter_complete_mixes()` to identify mixes 
    that may have more reagents added to them.
    """

    def __init__(self, components):

        def merge_partial_mixes(components):
            for component in components:
                if isinstance(component, PartialMix):
                    yield from component
                else:
                    yield component

        super().__init__(merge_partial_mixes(components))

class Levels:
    # - Dict with keys stored in heap, for fast access to top levels.
    # - Easy ability to fork.
    # - Each level a collection of reagents, mixes, and partial mixes
    # - Invariant: reagents that vaary together are grouped into partial mixes

    def __init__(self, components, combos):
        self._levels = levels = {}

        for component in components:
            n = count_combos(component, combos)
            levels.setdeault(n, set()).add(component)

        for n in levels:
            self._merge_identically_varying_components(n)

        self._keys = list(levels.keys()); heapify(self._keys)
        self._combos = combos

    def __len__(self):
        return len(self._levels)

    def copy(self):
        """
        Make a copy of this data structure that can be modified independently 
        from the original.

        This method assumes that all of the components making up each level 
        (e.g. reagents, mixes, partial mixes) are immutable, and thus can be 
        shared between copies.  It also assumes that the given combos are 
        will not be modified.
        """
        copy = self.__class__.__new__(self.__class__)
        copy._levels = copy(self._levels)
        copy._keys = copy(self._keys)
        copy._combos = combos
        return copy

    def peek(self):
        k = self._keys[0]
        return self._levels[k]

    def pop(self):
        k = heappop(self._keys)
        return self._levels.pop(k)

    def update(self, partial_mix):
        """
        Add the given partial mix to the indicated level, and remove all of the 
        components comprising that mix from whichever levels they appear in.
        """

        def drop_components(components):
            top_level = self.peek()

            remaining_level -= top_level - components
            remaining_components = components - top_level

            if remaining_level:
                k = self._keys[0]
                self._levels[k] = remaining_level
                assert not remaining_components
            else:
                self.pop()

            return remaining_components

        remaining_components = partial_mix
        while remaining_components:
            remaining_components = drop_components(remaining_components)

        n = count_combos(partial_mix, self._combos)

        if n in self._levels:
            self._levels[n].add(component)
            self._merge_identically_varying_components(n)
        else:
            self._levels[n] = {partial_mix}
            heappush(self._keys, n)

    def contains_mix(self, mix):
        for level in self._levels.values():
            if any(mix == x for x in iter_all_mixes(level)):
                return True
        return False

    def _merge_identically_varying_components(self, n):
        groups = []
        group_map = {}

        for a, b in combinations(self._levels[n], 2):
            n_ab = count_combos({a, b}, combos)
            if n == n_ab:
                try:
                    group = group_map[a]
                except KeyError:
                    try:
                        group = group_map[b]
                    except KeyError:
                        group = set()
                        groups.append(group)

                group |= {a, b}
                group_map[a] = group_map[b] = group

        for group in groups:
            self._levels[n] -= group
            self._levels[n].add(PartialMix(group))

    __repr__ == repr_from_init(
            attrs={'levels': '_levels'},
    )

# Rules:
# - Always combine reagents with the fewest variants
# - Reagents that vary together will be added together (either as a mix or as 
#   reagents)

# This seems more principled than what I currently have.

# Filters:
# - Check solvent/volume
# - Discard if has any mixes that aren't the most efficient

def iter_complete_mixes(levels, combos, filters=None):
    """
    Yield mixes (possibly nested) that include every component from all of the 
    given levels.

    You can think of this function as generating every possible way (in terms 
    of master mixes) to mix the given components.  In reality, some common 
    sense is used to avoid generating unreasonable results.  Specifically:

    - Reagents with the fewest number of variants are the first to be merged.
    - Reagents that vary together will always be added together.
    - If the same reagents can be mixed in multiple different ways, only keep 
      those that take the fewest pipetting steps.

    Filters can also be provided to perform more targeted pruning.
    """
    filters = filters or []

    def iter_possible_components(x):
        """
        Yield the components associated with the given argument.

        Reagents and mixes are yielded directly.  Partial mixes are yielded as 
        both a mix and a collection of reagents.
        """
        if isintance(x, PartialMix):
            yield [Mix(x)]
            yield x
        else:
            yield [x]

    def make_partial_mixes(pairs):
        for c0, c1 in pairs:
            c0s = iter_possible_components(c0)
            c1s = iter_possible_components(c1)
            for c0, c1 in product(c0s, c1s):
                yield PartialMix(*c0, *c1)

    def pass_all_filters(components):
        return all(f(components) for f in always_iterable(filters))

    candidates = [Levels(components, combos)]
    best_num_pipetting_steps = {}

    while candidates:
        # Do a breadth-first search, to get the most possible chances to prune 
        # sub-optimal solutions before yielding them.  This means treating 
        # `candidates` as a FIFO queue.
        
        levels = candidates.pop(0)
        top_level = levels.peek()

        if len(levels) == 1 and len(top_level) == 1:
            yield from top_level
            continue

        if len(top_level) > 1:
            pairs = itertools.combinations(top_level, 2)
        else:
            levels.pop()
            next_level = levels.peek()
            pairs = itertools.product(top_level, next_level)

        for partial_mix in make_partial_mixes(pairs):
            if not pass_all_filters(partial_mix):
                continue

            # Prune any solutions that arrive at any mix in a sub-optimal way.
            # I kind of want to break this logic into its own function, to 
            # logically separate searching from pruning and to make it easier 
            # to implement additional pruning strategies.  But do so would also 
            # add a lot of complexity, so I'll wait until I actually have the 
            # need to do so.

            suboptimal = set()

            for mix in iter_mixes(partial_mix):
                n_pipet = count_pipetting_steps(mix, combos, memo)
                n_best = best_num_pipetting_steps[mix.all_reagents]

                if n_pipet < n_best:
                    suboptimal.add(mix.all_reagents)
                    best_num_pipetting_steps[mix.all_reagents] = n_pipet
                if n_pipet > n_best:
                    break

            else:
                if suboptimal:
                    candidates = [
                            x for x in candidates
                            if not any(x.contains_mix(y) for y in suboptimal)
                    ]

                candidate = levels.copy().update(partial_mix)
                candidates.append(candidate)

