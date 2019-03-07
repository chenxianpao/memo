class TreeNode(object):
    def __init__(self, value=None, left=None, right=None):
        self.value = value
        self.left = left
        self.right = right


# 深度
def depth(tree):
    if tree is None:
        return 0
    left, right = depth(tree.left), depth(tree.right)
    return max(left, right) + 1


# 前序
def pre_order(tree):
    if tree is None:
        return
    print(tree.value)
    pre_order(tree.left)
    pre_order(tree.right)


# 中序
def mid_order(tree):
    if tree is None:
        return
    mid_order(tree.left)
    print(tree.value)
    mid_order(tree.right)


# 后序
def post_order(tree):
    if tree is None:
        return
    post_order(tree.left)
    post_order(tree.right)
    print(tree.value)


# 层次遍历
def level_order(tree):
    if tree is None:
        return
    q = list()
    q.append(tree)
    while q:
        current = q.pop(0)
        print(current.value)
        if current.left is not None:
            q.append(current.left)
        if current.right is not None:
            q.append(current.right)


# 按层次打印
def level2_order(tree):
    if tree is None:
        return
    q = list()
    q.append(tree)
    results = dict()
    level = 0
    current_level_num = 1
    next_level_num = 0
    d = list()
    while q:
        current = q.pop(0)
        current_level_num -= 1
        d.append(current.value)
        if current.left is not None:
            q.append(current.left)
            next_level_num += 1
        if current.right is not None:
            q.append(current.right)
            next_level_num += 1
        if current_level_num == 0:
            current_level_num = next_level_num
            next_level_num = 0
            results[level] = d
            d = list()
            level += 1
    print(results)


if __name__ == '__main__':
    tree = TreeNode('D', TreeNode('B', TreeNode('A'), TreeNode('C')),
                    TreeNode('E', right=TreeNode('G', TreeNode('F'))))
    """
        D
      B   E
    A  C    G
           F
    """
    print(depth(tree))
    # pre_order(tree)   # DBACEGF
    # mid_order(tree)     # ABCDEFG
    # post_order(tree)   # ACBFGED
    # level_order(tree)   # DBEACGF
    # level2_order(tree)  # {0: ['D'], 1: ['B', 'E'], 2: ['A', 'C', 'G'], 3: ['F']}