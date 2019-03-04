import asyncio
# def lazy_range(up_to):
#     """Generator to return the sequence of integers from 0 to up_to, exclusive."""
#     index = 0
#     while index < up_to:
#         yield index
#         index += 1

# iterator = lazy_range(10)
# for i in iterator:
#     print i


def jumping_range(up_to):
    """Generator for the sequence of integers from 0 to up_to, exclusive.

    Sending a value into the generator will shift the sequence by that amount.
    """
    index = 0
    while index < up_to:
        jump = yield index
        if jump is None:
            jump = 1
        index += jump


def lazy_range(up_to):
    """Generator to return the sequence of integers from 0 to up_to, exclusive."""

    # index = 0

    def gratuitous_refactor():
        index = 0
        while index < up_to:
            yield index
            index += 1

    yield from gratuitous_refactor()


def bottom():
    # Returning the yield lets the value that goes up the call stack to come right back
    # down.
    return (yield 42)


def middle():
    return (yield from bottom())


def top():
    return (yield from middle())


# Borrowed from http://curio.readthedocs.org/en/latest/tutorial.html.
@asyncio.coroutine
def countdown(number, n):
    while n > 0:
        print('T-minus', n, '({})'.format(number))
        yield from asyncio.sleep(1)
        n -= 1



if __name__ == '__main__':
    # iterator = lazy_range(10)
    # for x in iterator:
    #     print(x)
    #
    # # Get the generator.
    # gen = top()
    # value = next(gen)
    # print(value)  # Prints '42'.
    # try:
    #     value = gen.send(value * 2)
    # except StopIteration as exc:
    #     value = exc.value
    # print(value)  # Prints '84'.
    # iterator = jumping_range(5)
    # # 0 1 2 3 4
    # print(next(iterator))  # 0
    # print(iterator.send(2))  # 2
    # print(next(iterator))  # 3
    # print(iterator.send(-1))  # 2
    # for x in iterator:
    #     print(x)  # 3, 4
    loop = asyncio.get_event_loop()
    tasks = [
        asyncio.ensure_future(countdown("A", 2)),
        asyncio.ensure_future(countdown("B", 3))]
    loop.run_until_complete(asyncio.wait(tasks))
    loop.close()
