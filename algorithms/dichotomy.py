"""
1. 二分查找是有条件的，首先是有序，其次因为二分查找操作的是下标，所以要求是顺序表
2. 最优时间复杂度：O(1)
3. 最坏时间复杂度：O(logn)
"""
index = None


def binary_search_recursion(q_list, data, left, right):
    if left > right:
        return None
    mid = int((left + right) / 2)
    if q_list[mid] > data:
        return binary_search_recursion(q_list, data, left, mid)
    elif q_list[mid] < data:
        return binary_search_recursion(q_list, data, mid + 1, right)
    else:
        return mid


def binary_search_loop(q_list, data):
    left, right = 0, len(q_list) - 1
    while left <= right:
        mid = (left + right) // 2
        if q_list[mid] < data:
            left = mid + 1
        elif q_list[mid] > data:
            right = mid - 1
        else:
            return mid
    return None
# def binary_search(q_list, left, right, num):
#     if left > right:
#         return -1
#     mid = (left + right) // 2
#     if num < q_list[mid]:
#         right = mid - 1
#     elif num > q_list[mid]:
#         left = mid + 1
#     else:
#         return mid
#     return binary_search(q_list, left, right, num)

if __name__ == '__main__':
    import random
    lst = [random.randint(0, 10000) for _ in range(100000)]
    lst.sort()

    def test_recursion():
        binary_search_recursion(lst, 999, 0, len(lst) - 1)


    def test_loop():
        binary_search_loop(lst, 999)


    import timeit

    t1 = timeit.Timer("test_recursion()", setup="from __main__ import test_recursion")
    t2 = timeit.Timer("test_loop()", setup="from __main__ import test_loop")

    print("Recursion:", t1.timeit())
    print("Loop:", t2.timeit())
    #  Recursion: 4.383864044037182
    #  Loop: 2.415481715986971
