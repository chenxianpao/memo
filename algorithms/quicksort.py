def quicksort(q_list, start, end):
    if start < end:
        i, j = start, end
        base = q_list[i]
        while i < j:
            while (i < j) and q_list[j] >= base:
                j -= 1
            q_list[i] = q_list[j]
            while (i < j) and q_list[i] <= base:
                i += 1

            q_list[j] = q_list[i]
        q_list[i] = base
        # print(i)
        # print(j)
        # if i == j:
        #     print(i)
        quicksort(q_list, start, i - 1)
        quicksort(q_list, i + 1, end)
    return q_list


def quick_sort(array):
    if len(array) < 2:
        return array

    stack = list()
    stack.append(len(array) - 1)
    stack.append(0)
    while stack:
        l = stack.pop()
        r = stack.pop()
        index = partition(array, l, r)
        if l < index - 1:
            stack.append(index - 1)
            stack.append(l)
        if r > index + 1:
            stack.append(r)
            stack.append(index + 1)


def partition(array, start, end):
    pivot = array[start]
    while start < end:
        while start < end and array[end] >= pivot:
            end -= 1
        array[start] = array[end]

        while start < end and array[start] <= pivot:
            start += 1
        array[end] = array[start]
    array[start] = pivot
    return start

a = [1, 4, 7, 1, 5, 5, 3, 85, 34, 75, 23, 75, 2, 0]
# print(quicksort(a, 0, len(a) - 1))
quick_sort(a)
print(a)