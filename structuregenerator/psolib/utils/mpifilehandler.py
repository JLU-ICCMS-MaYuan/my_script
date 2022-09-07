# A FileHandler to open and write to a file with mpi4py.MPI.File methods
# Modified and fix from gist <https://gist.github.com/JohnCEarls/8172807>

import logging
import logging.handlers
import os
import sys

from mpi4py import MPI
from mpi4py.futures import MPICommExecutor


class MPIFileHandler(logging.FileHandler):
    def __init__(
        self,
        filename,
        mode=MPI.MODE_WRONLY | MPI.MODE_CREATE | MPI.MODE_APPEND,
        encoding=None,
        delay=0,
        comm=MPI.COMM_WORLD,
    ):
        encoding = None
        self.baseFilename = os.path.abspath(filename)
        self.mode = mode
        self.encoding = encoding
        self.comm = comm
        if delay:
            # We don't open the stream, but we still need to call the
            # Handler constructor to set level, formatter, lock etc.
            logging.Handler.__init__(self)
            self.stream = None
        else:
            logging.StreamHandler.__init__(self, self._open())

    def _open(self):
        # stream = MPILogFile.Open( self.comm, self.baseFilename, self.mode )
        stream = MPI.File.Open(self.comm, self.baseFilename, self.mode)
        stream.Set_atomicity(True)
        return stream

    def close(self):
        if self.stream:
            self.stream.Sync()
            self.stream.Close()
            self.stream = None

    def emit(self, record):
        record.rank = MPI.COMM_WORLD.Get_rank()
        record.size = MPI.COMM_WORLD.Get_size()
        if self.stream is None:
            if self.mode != 'w' or not self._closed:
                self.stream = self._open()
        if self.stream is sys.stderr:
            logging.StreamHandler.emit(self, record)
        else:
            try:
                msg = self.format(record) + self.terminator
                msg = msg.encode()
                stream = self.stream
                stream.Write_shared(msg)
                self.flush()
            except RecursionError:
                raise
            except Exception:
                self.handleError(record)


if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    logger = logging.getLogger("node[%i]" % comm.rank)
    logger.setLevel(logging.DEBUG)

    mh = MPIFileHandler("test.log")
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(rank)s/%(size)s - %(levelname)s - %(message)s'
    )
    mh.setFormatter(formatter)

    logger.addHandler(mh)
    # 'application' code
    logger.debug('debug message')
    logger.info('info message')
    logger.warning('warn message')
    logger.error('error message')
    logger.critical('critical message')

    with MPICommExecutor(comm, 0) as executor:
        if executor is not None:
            logger.debug('In executor debug message')
            logger.info('In executor info message')
            logger.warning('In executor warn message')
            logger.error('In executor error message')
            logger.critical('In executor critical message')
