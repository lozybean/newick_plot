#!/usr/bin/env python
# -*- coding: utf-8 -*- \#
"""
@author = 'liangzb'
@date = '2016-01-08'
a tree data structure, referenced from ete3
"""

from collections import defaultdict, deque
from ete3 import NodeStyle
from ete3.coretype.tree import TreeError


class Tree(object):
    def __init__(self, root=None):
        if root is None:
            self.__root_node = TreeNode('root', level='root')
        else:
            self.__root_node = root
        self.__root_node.tree = self
        self.nodes = defaultdict(TreeNode)
        self.nodes['root'] = self.root

    def __str__(self):
        return str(self.root) + ';'

    def __iter__(self):
        for node in self.root.traverse():
            yield node

    def iter_level(self, level):
        for node in self.root.traverse():
            if node.level == level:
                yield node

    def get_level(self, level):
        return [node for node in self.iter_level(level)]

    def adjust_profile(self):
        for node in self:
            if not node.get_sisters():
                continue
            percent = node.profile / sum(map(lambda s: s.profile, node.get_sisters(include=True)))
            if node.up.profile:
                node.profile = percent * node.up.profile
            else:
                node.profile = percent

    @property
    def root(self):
        return self.__root_node


class TreeNode(object):
    MIN_SIZE = 5

    def _iter_descendants_postorder(self, is_leaf_fn=None):
        to_visit = [self]
        if is_leaf_fn is not None:
            _leaf = is_leaf_fn
        else:
            _leaf = self.__class__.is_leaf

        while to_visit:
            node = to_visit.pop(-1)
            try:
                node = node[1]
            except TypeError:
                # PREORDER ACTIONS
                if not _leaf(node):
                    # ADD CHILDREN
                    to_visit.extend(reversed(node.children + [[1, node]]))
                else:
                    yield node
            else:
                # POSTORDER ACTIONS
                yield node

    def _iter_descendants_levelorder(self, is_leaf_fn=None):
        """
        Iterate over all desdecendant nodes.
        """
        to_visit = deque([self])
        while len(to_visit) > 0:
            node = to_visit.popleft()
            yield node
            if not is_leaf_fn or not is_leaf_fn(node):
                to_visit.extend(node.children)

    def _iter_descendants_preorder(self, is_leaf_fn=None):
        """
        Iterator over all descendant nodes.
        """
        to_visit = deque()
        node = self
        while node is not None:
            yield node
            if not is_leaf_fn or not is_leaf_fn(node):
                to_visit.extendleft(reversed(node.children))
            try:
                node = to_visit.popleft()
            except IndexError:
                node = None

    def __init__(self, name, level=None):
        self.name = name
        self.up = None
        self.children = []
        self.__min_size = self.MIN_SIZE
        self.__style = None
        self.__level = level
        self.__profile = 0
        self.__dist = 1
        self.__size = self.MIN_SIZE
        self.__tree = None

    def __str__(self):
        if self.children:
            result = '('
            result += ','.join(map(str, self.children))
            result += ')%s' % self.name
        else:
            result = self.name

        return result

    # Topology management
    def add_child(self, child=None, name=None, level=None):
        """
        Adds a new child to this node. If child node is not suplied
        as an argument, a new node instance will be created.

        :argument None child: the node instance to be added as a child.
        :argument None name: the name that will be given to the child.
        :argument None level: the taxanomy level that will be given to the child.

        :returns: The child node instance

        """
        if child is None:
            child = self.__class__()
        if name is not None:
            child.name = name
        if level is not None:
            child.level = level
        if child not in self.children:
            self.children.append(child)
        child.up = self
        return child

    def add_sister(self, sister=None, name=None, level=None):
        """
        Adds a sister to this node. If sister node is not supplied
        as an argument, a new TreeNode instance will be created and
        returned.
        :param sister: the node instance to be added as a sister
        :param name: the name that will be added to the sister
        :param level: the level that will be added to the sister
        :return: None
        """
        if self.up is None:
            raise TreeError("A parent node is required to add a sister")
        else:
            return self.up.add_child(child=sister, name=name, level=level)

    def get_children(self):
        """
        Returns an independent list of node's children.
        """
        return [ch for ch in self.children]

    def get_sisters(self, include=False):
        """
        Returns an indepent list of sister nodes.
        :param include: determine weather self is in
        """
        if self.up is not None:
            if include:
                return [ch for ch in self.up.children]
            else:
                return [ch for ch in self.up.children if ch != self]
        else:
            return []

    def get_same_level(self):
        same_level = []
        for node in self.tree:
            if node.level == self.level:
                same_level.append(node)
        return same_level

    def is_leaf(self):
        """
        Return True if current node is a leaf.
        """
        return len(self.children) == 0

    @property
    def min_size(self):
        return self.__min_size

    @min_size.setter
    def min_size(self, value):
        try:
            self.__min_size = float(value)
        except ValueError:
            pass

    @property
    def style(self):
        return self.__style

    @style.setter
    def style(self, value):
        if isinstance(value, NodeStyle):
            self.__style = value

    @property
    def tree(self):
        return self.__tree

    @tree.setter
    def tree(self, t):
        if isinstance(t, Tree):
            self.__tree = t

    @property
    def level(self):
        return self.__level

    @level.setter
    def level(self, value):
        self.__level = value

    @property
    def profile(self):
        return self.__profile

    @profile.setter
    def profile(self, value):
        self.__profile = value
        self.size = self.profile ** 0.5 * 50
        self.__dist = len(self.name) * 7 + self.size + 4

    @property
    def size(self):
        return self.__size

    @size.setter
    def size(self, value):
        try:
            value = float(value)
        except ValueError:
            pass
        if value >= self.min_size:
            self.__size = value

    @property
    def dist(self):
        return self.__dist

    @dist.setter
    def dist(self, value):
        try:
            self.__dist = float(value)
        except ValueError:
            pass

    @property
    def branch_length(self):
        def f(node):
            return node.dist

        level_dists = map(f, self.get_same_level())
        if level_dists:
            max_dist = max(level_dists)
            return max_dist - self.size
        return self.dist

    def traverse(self, strategy="levelorder", is_leaf_fn=None):
        """
        Returns an iterator to traverse the tree structure under this
        node.

        :argument "levelorder" strategy: set the way in which tree
           will be traversed. Possible values are: "preorder" (first
           parent and then children) 'postorder' (first children and
           the parent) and "levelorder" (nodes are visited in order
           from root to leaves)

        :argument None is_leaf_fn: If supplied, ``is_leaf_fn``
           function will be used to interrogate nodes about if they
           are terminal or internal. ``is_leaf_fn`` function should
           receive a node instance as first argument and return True
           or False. Use this argument to traverse a tree by
           dynamically collapsing internal nodes matching
           ``is_leaf_fn``.
        """
        if strategy == "preorder":
            return self._iter_descendants_preorder(is_leaf_fn=is_leaf_fn)
        elif strategy == "levelorder":
            return self._iter_descendants_levelorder(is_leaf_fn=is_leaf_fn)
        elif strategy == "postorder":
            return self._iter_descendants_postorder(is_leaf_fn=is_leaf_fn)
