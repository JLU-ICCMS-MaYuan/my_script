"""A decorator that limits the running time of the decorated function

Usage:
    from time_limit import timelimit, TimeoutException

    @timelimit(exectime=None)
    def func():
        ...

    result = func(exectime=exectime)

    # ====== Normal Use ===============================================
    if isinstance(result, TimeoutException):
        print('Reached time limit)

    # ====== When kwargs cannot passed to func ========================

    from functools import partial

    part_func = partial(func, exectime=exectime)
    map(part_func, *iterable)

An additional kwargs 'exectime' can be passed in by decorator or add
when calling decorated function, while decorator args is priority.

Set exectime=False means no time limit.
Return `TimeoutException` if reach time limit.
"""

# Source from <https://blog.csdn.net/qq_42709514/article/details/84001494>
# With license [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
# Change:
#   1. Format on 'if ... is not None'
#   2. Add functools.wrap to make decorator can be pickled
#   3. Add a custom TimeoutException
#   4. print error traceback from decorated function self

import ctypes
import inspect
import threading
import time
import traceback
from functools import wraps
from threading import Thread


class TimeoutException(Exception):
    def __init__(self, errmsg=""):
        self.errmsg = errmsg

    def __str__(self):
        return 'Timeout ! %s' % self.errmsg


def _async_raise(tid, exctype):
    """raises the exception, performs cleanup if needed"""
    tid = ctypes.c_long(tid)
    if not inspect.isclass(exctype):
        exctype = type(exctype)
    res = ctypes.pythonapi.PyThreadState_SetAsyncExc(tid, ctypes.py_object(exctype))
    if res == 0:
        raise ValueError("invalid thread id")
    elif res != 1:
        # """if it returns a number greater than one, you're in trouble,
        # and you should call it again with exc=NULL to revert the effect"""
        ctypes.pythonapi.PyThreadState_SetAsyncExc(tid, None)
        raise SystemError("PyThreadState_SetAsyncExc failed")


def stop_thread(thread):
    _async_raise(thread.ident, SystemExit)


def timelimited(exectime=None):
    def decorator(function):
        @wraps(function)
        def decorator2(*args, **kwargs):
            # print(args, kwargs)
            time_out = (
                exectime if exectime is not None else kwargs.pop("exectime", None)
            )
            if time_out is None:
                return Exception("exec time param missing %s, %s" % (args, kwargs))

            class TimeLimited(threading.Thread):
                def __init__(self, _error=None):
                    Thread.__init__(self)
                    self._error = _error
                    self._result = None

                def run(self):
                    try:
                        result = function(*args, **kwargs)
                        if result is None:
                            self._result = True
                        else:
                            self._result = result
                    except Exception as err:
                        traceback.print_exc()
                        self._error = Exception(err)

            t = TimeLimited()
            t.setDaemon(True)
            t.start()
            if time_out is not False:
                t.join(time_out)
                if isinstance(t._error, Exception):
                    return t._error
                else:
                    if t._result is None:
                        stop_thread(t)
                        return TimeoutException()
                    else:
                        return t._result

        return decorator2

    return decorator


@timelimited(exectime=None)  # False: no restrit
def func():
    time.sleep(1)
    return 'Finished'


if __name__ == "__main__":
    print(func(exectime=1.1))
    print(threading.enumerate())
