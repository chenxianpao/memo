def consumer(name):
    print("要开始啃骨头了...")
    while True:
        print("\033[31;1m[consumer] %s\033[0m " % name)
        bone = yield
        print("[%s] 正在啃骨头 %s" % (name, bone))


def producer(obj1, obj2):
    obj1.send(None)
    obj2.send(None)
    n = 0
    while n < 5:
        n += 1
        print("\033[32;1m[producer]\033[0m 正在生产骨头 %s" % n)
        obj1.send(n)
        obj2.send(n)


if __name__ == '__main__':
    con1 = consumer("消费者A")
    con2 = consumer("消费者B")
    producer(con1, con2)
