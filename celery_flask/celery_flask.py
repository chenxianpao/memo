import random
import time
from flask import Flask
from flask import jsonify
from celery import Celery
from celery.utils.log import get_task_logger


def make_celery(app):
    celery = Celery(app.name,
                    broker=app.config['CELERY_BROKER_URL'],
                    backend=app.config['CELERY_RESULT_BACKEND'])
    celery.conf.update(app.config)
    TaskBase = celery.Task
    class ContextTask(TaskBase):
        abstract = True
        def __call__(self, *args, **kwargs):
            with app.app_context():
                return TaskBase.__call__(self, *args, **kwargs)
    celery.Task = ContextTask
    return celery

app = Flask(__name__)
app.config.update(
    CELERY_BROKER_URL='amqp://cxp:123@localhost:5672/myvhost',
    CELERY_RESULT_BACKEND='rpc://'
)

celery = make_celery(app)


@celery.task
def add(x, y):
    return x + y


@celery.task(bind=True)
def long_task(self):
    verb = ['starting up', 'booting', 'repairing', 'loading', 'checking']
    adjective = ['master', 'radiant', 'silent', 'harmonic', 'fast']
    noun = ['solar array', 'particle reshaper', 'cosmic ray', 'orbiter', 'bit']
    msg = ''
    total = random.randint(10, 50)
    for i in range(total):
        if not msg or random.random() < 0.25:
            msg = '{0} {1} {2}...'.format(random.choice(verb),
                                          random.choice(adjective),
                                          random.choice(noun))
        self.update_state(state='PROGRESS', meta={'current': i,
                                                  'total': total,
                                                  'status': msg})
        time.sleep(1)
    return {'current': 100, 'total': 100, 'status': 'task completed',
            'result': 42}


@app.route("/")
def add():
    print "begin"
    result = add.delay(10, 20)
    return "aaa"


@app.route("/longtask")
def longtask():
    task = long_task.apply_async()
    return task.id


@app.route('/status/<task_id>')
def taskstatus(task_id):
    task = long_task.AsyncResult(task_id)
    if task.state == 'PENDING':
        print 'job did not start yet'
        ret = {
            'state': task.state,
            'current': 0,
            'total': 1,
            'status': 'Pending'
        }
    elif task.state != 'FAILURE':
        ret = {
            'state': task.state,
            'current': task.info.get('current', 0),
            'total': task.info.get('total', 1),
            'status': task.info.get('status', '')
        }
        if 'result' in task.info:
            ret['result'] = task.info['result']
    else:
        ret = {
            'state': task.state,
            'current': 1,
            'total': 1,
            'status': str(task.info)
        }
    return jsonify(ret)


@app.route('/stop/<task_id>')
def stoptask(task_id):
    # ret = long_task.AsyncResult(task_id, terminate=True)
    celery.control.revoke(task_id, terminate=True)
    return jsonify("stop ok")


if __name__ == '__main__':
    app.run()