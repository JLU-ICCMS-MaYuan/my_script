#!/usr/bin/env python


class FileNotFoundError(Exception):
    def __init__(self, errmsg=""):
        self.errmsg = errmsg

    def __str__(self):
        return 'File %s not Found ! ' % self.errmsg


class TimeoutException(Exception):
    def __init__(self, errmsg=""):
        self.errmsg = errmsg

    def __str__(self):
        return 'Timeout ! %s' % self.errmsg