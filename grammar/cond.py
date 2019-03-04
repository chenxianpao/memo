import threading

cond = threading.Condition()


class Test(threading.Thread):
    def __init__(self):
        threading.Thread.__init__(self)

    def run(self):
        while True:
            cond.acquire()
            cond.wait(60)
            cond.release()
            print("Test Thread in")


if __name__ == '__main__':
    cond.acquire()
    cond.notify()
    cond.release()
