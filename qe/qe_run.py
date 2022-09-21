import os
import logging
from argparse import ArgumentParser
from pathlib import Path
from itertools import chain

from config import config
from qe_inputpara import * 
from qe_writeinput import qe_writeinput
from qe_writesubmit import qe_writesubmit
from qe_submitjob import qe_submitjob

class qe_relax:

    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config = config(args).read_config()

        # prepare the POSCAR POTCAR  
        self.relax_inputpara  = qe_inputpara.init_from_config(self._config)

        # init the INCAR
        self.qe_writeincar  = qe_writeinput.init_from_relaxinput(self.relax_inputpara)

        # init the submit job script
        self.qe_writesubmit = qe_writesubmit.init_from_relaxinput(self.relax_inputpara)

        # submit the job
        self.qe_submitjob   = qe_submitjob.init_from_relaxinput(self.relax_inputpara)


class qe_scf:

    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config = config(args).read_config()

        # prepare the POSCAR POTCAR  
        self.scf_inputpara  = qe_inputpara.init_from_config(self._config)

        # init the INCAR
        self.qe_writeincar  = qe_writeinput.init_from_scfinput(self.scf_inputpara)

        # init the submit job script
        self.qe_writesubmit = qe_writesubmit.init_from_scfinput(self.scf_inputpara)

        # submit the job
        self.qe_submitjob   = qe_submitjob.init_from_scfinput(self.scf_inputpara)


class qe_phono:

    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config = config(args).read_config()

        # prepare the POSCAR POTCAR  
        self.phono_inputpara = qephono_inputpara.init_from_config(self._config)

        # init the INCAR
        self.qe_writeincar  = qe_writeinput.init_from_phonoinput(self.phono_inputpara)

        # init the submit job script
        self.qe_writesubmit = qe_writesubmit.init_from_phonoinput(self.phono_inputpara)

        # submit the job
        self.qe_submitjob   = qe_submitjob.init_from_phonoinput(self.phono_inputpara)

class qe_superconduct:

    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config = config(args).read_config()

        # prepare the POSCAR POTCAR  
        self.sc_inputpara  = qesc_inputpara.init_from_config(self._config)
        # init the INCAR
        self.qe_writeincar  = qe_writeinput.init_from_scinput(self.sc_inputpara)

        # init the submit job script
        self.qe_writesubmit = qe_writesubmit.init_from_scinput(self.sc_inputpara)

        # submit the job
        self.qe_submitjob   = qe_submitjob.init_from_scinput(self.sc_inputpara)