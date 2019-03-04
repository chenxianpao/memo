import asyncio
import threading

@asyncio.coroutine
def hello():
    print("Hello World")
    r = yield from asyncio.sleep(1)
    print("Hello Again")


@asyncio.coroutine
def hello_v1():
    print('Hello World!(%s)' % threading.current_thread())
    yield from asyncio.sleep(1)
    print('Hello Again (%s)' % threading.current_thread())


@asyncio.coroutine
def wget(host):
    print('wget %s...' % host)
    connect = asyncio.open_connection(host, 80)
    reader, writer = yield from connect
    header = 'GET / HTTP/1.0\r\nHost: %s\r\n\r\n' % host
    writer.write(header.encode('utf-8'))
    yield from writer.drain()
    while True:
        line = yield from reader.readline()
        if line == b'\r\n':
            break
        print('%s header > %s' % (host, line.decode('utf-8').rstrip()))
    # Ignore the body, close the socket
    writer.close()



loop = asyncio.get_event_loop()
tasks = [hello_v1(), hello_v1()]
tasks_get = [wget(host) for host in ['www.sina.com.cn', 'www.sohu.com', 'www.163.com']]
# loop.run_until_complete(asyncio.wait(tasks))
loop.run_until_complete(asyncio.wait(tasks_get))
loop.close()
