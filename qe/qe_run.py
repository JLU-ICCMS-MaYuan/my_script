import os
import logging
from argparse import ArgumentParser

from qe.config import config
from qe.qe_inputpara import * 
from qe.qe_writeinput import qe_writeinput
from qe.qe_writesubmit import qe_writesubmit
from qe.qe_submitjob import qe_submitjob



logger = logging.getLogger(__name__)

logger.info("This is an info message")
logger.debug("This is a debug message")

def check_pid_jobid(ids: list, submit_job_system):
    if submit_job_system == "bash":
        i = 0
        while True:
            osawk = """sleep 5 | ps -ef | grep -E "pw.x" |  grep -v grep | awk '{print $2}'""" 
            _jobids = os.popen(osawk).read()  # return a string; such as '423423\n324233\n423424\n'
            jobids = _jobids.strip("\n").split("\n")
            for id in ids:
                if id not in jobids:
                    i += 1
                if i == len(ids):
                    return
    elif submit_job_system == "slurm":
        i = 0
        while True:
            osawk = """sleep 5 | squeue | awk '{print $1}'""" 
            _jobids = os.popen(osawk).read()  # return a string; such as '423423\n324233\n423424\n'
            jobids = _jobids.strip("\n").split("\n")
            for id in ids:
                if id not in jobids:
                    i += 1
                if i == len(ids):
                    return
    elif submit_job_system == "pbs":
        i = 0
        while True:
            osawk = """sleep 5 | qstat | awk '{print $1}' | cut -d . -f1""" 
            _jobids = os.popen(osawk).read()  # return a string; such as '423423\n324233\n423424\n'
            jobids = _jobids.strip("\n").split("\n")
            for id in ids:
                if id not in jobids:
                    i += 1
                if i == len(ids):
                    return


class qe_relax:

    def __init__(self, args: ArgumentParser) -> None:
        
        # read input para
        self._config = config(args).read_config()

        # prepare input parameter
        self.relax_inputpara  = qe_inputpara.init_from_config(self._config)

        # init the input
        self.qe_writeinput  = qe_writeinput(self.relax_inputpara)
        inputfilename = self.qe_writeinput.writeinput()

        # init the submit job script
        self.qe_writesubmit = qe_writesubmit(self.relax_inputpara)
        jobname = self.qe_writesubmit.write_submit_scripts(inputfilename)
        
        # submit the job. If we didn't set the parameter of `queue`, it will be set `None` in `qe_inputpara`
        self.qe_submitjob = qe_submitjob(self.relax_inputpara)
        if self.relax_inputpara.queue is not None:
            self.qe_submitjob.submit_mode1(inputfilename, jobname)


class qe_scf:

    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config = config(args).read_config()

        # prepare input parameter
        self.scf_inputpara  = qe_inputpara.init_from_config(self._config)

        # init the input
        self.qe_writeinput  = qe_writeinput(self.scf_inputpara)
        inputfilename = self.qe_writeinput.writeinput()

        # init the submit job script
        self.qe_writesubmit = qe_writesubmit(self.scf_inputpara)
        jobname = self.qe_writesubmit.write_submit_scripts(inputfilename)

        # submit the job
        self.qe_submitjob = qe_submitjob(self.scf_inputpara)
        if self.scf_inputpara.queue is not None:
            self.qe_submitjob.submit_mode1(inputfilename, jobname)


class qe_phono:

    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config = config(args).read_config()

        # prepare input parameter
        self.phono_inputpara = qephono_inputpara.init_from_config(self._config)
        # init the input
        self.qe_writeinput   = qe_writeinput(self.phono_inputpara)
        # init the submit job script
        self.qe_writesubmit = qe_writesubmit(self.phono_inputpara)
        # submit the job
        self.qe_submitjob = qe_submitjob(self.phono_inputpara)
        
        if self.phono_inputpara.mode == "nosplit":
            inputfilename = self.qe_writeinput.writeinput()
            jobnames = self.qe_writesubmit.write_submit_scripts(inputfilename)
            if self.phono_inputpara.queue is not None:
                self.qe_submitjob.submit_mode2(inputfilename, jobnames)
        elif self.phono_inputpara.mode == "split_dyn0" or self.phono_inputpara.mode == "split_assignQ":
            inputfilename = self.qe_writeinput.writeinput()
            jobnames = self.qe_writesubmit.write_submit_scripts(inputfilename)
            if self.phono_inputpara.queue is not None:
                self.qe_submitjob.submit_mode3(inputfilename, jobnames)
        elif self.phono_inputpara.mode == "merge":
            self.phono_inputpara.merge(self.phono_inputpara.work_path)
        elif self.phono_inputpara.mode == "phonobandwidthsdata":
            gauss = self.phono_inputpara.gauss
            qpoints_freqs, q_number, freq_number = self.phono_inputpara.get_phono_freq()
            phononwidth = self.phono_inputpara.get_gam_lines(gauss, q_number, freq_number)
            self.phono_inputpara.merge_qp_freq_width(qpoints_freqs, phononwidth)
            logger.info("You can use `qp_freq_width.csv` to plot phonon-band")
        elif self.phono_inputpara.mode == "phonobanddata":
            inputfilename = self.qe_writeinput.writeinput()
            jobnames = self.qe_writesubmit.write_submit_scripts(inputfilename)
            self.qe_submitjob.submit_mode1(inputfilename, jobnames)
            logger.info("You can use `prefix`.phonon.bands.dat to plot phonon-band")
        elif self.phono_inputpara.mode == "phonodosdata":
            self.phono_inputpara.get_phonodos()
            logger.info("You can use `phdos_proj2eles.csv` to plot phonon-DOS")
        elif self.phono_inputpara.mode == "gibbsvb":
            self.phono_inputpara.get_gibbs_from_phtdos()
            self.phono_inputpara.get_gibbs_from_freq()
        elif self.phono_inputpara.mode == "hspp":
            logger.info("Get hspp from matdyn")
            self.phono_inputpara.read_hspp_in_matdyn()
            logger.info("Get hspp from ASE package")
            self.phono_inputpara.get_hspp(autoselect=self.phono_inputpara.autoselect)
        else:
            inputfilename = self.qe_writeinput.writeinput()
            jobnames = self.qe_writesubmit.write_submit_scripts(inputfilename)
            logger.debug("If mode=matdyn, please remerber process phono spectrum after finish matdyn.x computation")
            logger.debug("The reason is `systemname`.freq, `systemname`.freq.gq and gam.lines will be rewrited when you calcuate phonodos")
            logger.debug("In most cases, the number of q-points is different between matdyn.in(high symmetry path qpoints sample) and phonondos.in(even qpoints sample)")
            self.qe_submitjob.submit_mode1(inputfilename, jobnames)
        

class qe_eletron:

    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config = config(args).read_config()

        # prepare input parameter
        self.eletron_inputpara = qeeletron_inputpara.init_from_config(self._config)
        self.qe_writeinput = qe_writeinput(self.eletron_inputpara)
        self.qe_writesubmit = qe_writesubmit(self.eletron_inputpara)
        self.qe_submitjob  = qe_submitjob(self.eletron_inputpara)
        
        if self.eletron_inputpara.mode == "hspp":
            self.eletron_inputpara.get_hspp(autoselect=self.eletron_inputpara.autoselect)
            eV_To_Ry = 0.0734986443513116
            ef_scffit,  ef_scf  = self.get_fermi_energy()
            nef_scffit, nef_scf = self.get_Nef(ef_scffit, ef_scf)
            logger.info("    fermi_energy = {:<8.4f} eV in scffit.out, N(Ef) = {:<8.4f} states/eV/(Unit Cell) = {:<10.6f} states/spin/Ry/(Unit Cell)".format(ef_scffit, nef_scffit, nef_scffit/2/eV_To_Ry))
            logger.info("    fermi_energy = {:<8.4f} eV in scf.out,    N(Ef) = {:<8.4f} states/eV/(Unit Cell) = {:<10.6f} states/spin/Ry/(Unit Cell)".format(ef_scf, nef_scf, nef_scf/2/eV_To_Ry))
        elif self.eletron_inputpara.mode == "eleband": 
            logger.debug("    !!!!!!!!!! Remember to run pw.x to get eleband.out before you run bands.x") 
            logger.debug("    Run bands.x to get eleband.dat and  eleband.dat.gnu")
            logger.debug("    eleband.dat.gnu can be used in origin to plot-eletronband")
            inputfilename1 = self.qe_writeinput.writeinput(mode="eleband")
            inputfilename2 = self.qe_writeinput.writeinput(mode="elebanddata")
            inputfilename3 = self.qe_writeinput.writeinput(mode="elebandprojdata")
            jobname = self.qe_writesubmit.write_submit_scripts(inpufilename=[inputfilename1, inputfilename2, inputfilename3])
            if self.eletron_inputpara.queue is not None:
                self.qe_submitjob.submit_mode1(inputfilename=inputfilename1, jobname=jobname)
        elif self.eletron_inputpara.mode == "eledos": 
            logger.debug("!!!!!!!!!! Remember to run pw.x to get nscf.out before you run dos.x and projwfc.x") 
            logger.debug("Run dos.x to get tdos and and run projwfc.x to get pdos")
            inputfilename1 = self.qe_writeinput.writeinput(mode="nscf")
            inputfilename2 = self.qe_writeinput.writeinput(mode="eletdos")
            inputfilename3 = self.qe_writeinput.writeinput(mode="elepdos")
            jobname = self.qe_writesubmit.write_submit_scripts(inpufilename=[inputfilename1, inputfilename2, inputfilename3])
            if self.eletron_inputpara.queue is not None:
                self.qe_submitjob.submit_mode1(inputfilename=inputfilename1, jobname=jobname)
            self.get_fermi_energy()
        elif self.eletron_inputpara.mode == "eleproperties":
            inputfilename1 = self.qe_writeinput.writeinput(mode="eleband")
            inputfilename2 = self.qe_writeinput.writeinput(mode="elebanddata")
            inputfilename3 = self.qe_writeinput.writeinput(mode="elebandprojdata")
            inputfilename4 = self.qe_writeinput.writeinput(mode="nscf")
            inputfilename5 = self.qe_writeinput.writeinput(mode="eletdos")
            inputfilename6 = self.qe_writeinput.writeinput(mode="elepdos")
            jobname = self.qe_writesubmit.write_submit_scripts([inputfilename1, inputfilename2, inputfilename3, inputfilename4, inputfilename6, inputfilename6])
            if self.eletron_inputpara.queue is not None:
                self.qe_submitjob.submit_mode1(inputfilename1, jobname)
            self.get_fermi_energy()
        else:
            # write input parameter
            inputfilename = self.qe_writeinput.writeinput()
            # init the submit job script
            jobname = self.qe_writesubmit.write_submit_scripts(inputfilename)
            # submit the job
            self.qe_submitjob = qe_submitjob(self.eletron_inputpara)
            if self.eletron_inputpara.queue is not None:
                self.qe_submitjob.submit_mode1(inputfilename, jobname)
            self.get_fermi_energy()

    def get_fermi_energy(self):
        scffit_out_path = self.eletron_inputpara.work_path.joinpath("scffit.out")
        if scffit_out_path.exists():
            ef_scffit = float(os.popen(f'grep "Fermi energy" {scffit_out_path}').read().split()[4])
        else:
            ef_scffit = None 
        scf_out_path = self.eletron_inputpara.work_path.joinpath("scf.out")
        if scf_out_path.exists():
            ef_scf = float(os.popen(f'grep "Fermi energy" {scf_out_path}').read().split()[4])
        else:
            ef_scf = None
        return ef_scffit, ef_scf

    def get_Nef(self, ef_scffit, ef_scf):
        eletdos_path = self.eletron_inputpara.work_path.joinpath(self.eletron_inputpara.system_name+".tdos")
        if not eletdos_path.exists():
            logger.warning(f"Sorry, {self.eletron_inputpara.system_name}_phono.dos doesn't exist !")
            nef_scffit, nef_scf = 0, 0
            return nef_scffit, nef_scf
        else:
            eletdos = pd.read_table(
                eletdos_path,
                skiprows=1,  # skiprows=1：跳过文件的第一行，即不将其作为数据的一部分进行读取。
                header=None, # header=None：不将文件的第一行作为列名，而将其视为数据。
                sep='\s+'    # sep='\s+'：使用正则表达式 \s+ 作为列之间的分隔符，表示一个或多个空格字符。
                )
            eletdos.columns = ['E(eV)', 'dos(E)', 'Int dos(E)']
            idx_scffit = (eletdos['E(eV)'] - ef_scffit).abs().idxmin()
            nef_scffit = eletdos['dos(E)'][idx_scffit]
            idx_scf    = (eletdos['E(eV)'] - ef_scf).abs().idxmin()
            nef_scf    = eletdos['dos(E)'][idx_scf]
            return nef_scffit, nef_scf


class qe_superconduct:

    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config = config(args).read_config()

        # prepare input parameter
        self.sc_inputpara  = qesc_inputpara.init_from_config(self._config)
        # 这里注意初始化qesc_inputpara时会初始化两个非常重要的参数self.sc_inputpara.gaussid 和 self.sc_inputpara.gauss
        # 这两个参数时qephono_inputpara中继承过来的

        # 当电荷屏蔽常数为0.10 和 0.13 时
        results = []
        for screen_constant in self.sc_inputpara.screen_constant:
            # McAD
            self.qe_writeinput  = qe_writeinput(self.sc_inputpara)
            self.qe_submitjob = qe_submitjob(self.sc_inputpara)

            inputfilename = self.qe_writeinput.write_lambda_in(self.sc_inputpara.work_path, screen_constant)
            if self.sc_inputpara.queue is not None:
                self.qe_submitjob.submit_mode0(inputfilename, dotx_file="lambda.x")

            # 因此初始化qesc_inputpara时还会根据已经获得的self.sc_inputpara.gaussid 和 self.sc_inputpara.gauss， 
            # 去读取相应展宽的 alpha2F.dat 和a2F.dos*
            self.a2Fdos_data = self.sc_inputpara.get_a2Fdos_data(self.sc_inputpara.gaussid)
            self.alpha2Fdat_data = self.sc_inputpara.get_alpha2Fdat_data(self.sc_inputpara.gaussid, self.sc_inputpara.gauss)

            # lambda
            Lambda_byalpha2f = self.sc_inputpara.get_lambda_from_alpha2f_single_broadening(self.alpha2Fdat_data, self.sc_inputpara.gaussid, self.sc_inputpara.gauss)
            Lambda_bya2Fdos  = self.sc_inputpara.get_lambda_from_a2fdos_single_broadening(self.a2Fdos_data, self.sc_inputpara.gaussid)
            
            # wlog
            wlog_bya2Fdos = self.sc_inputpara.get_wlog_from_a2fdos_single_broadening(self.a2Fdos_data, Lambda_bya2Fdos)

            # w2
            w2_bya2Fdos   = self.sc_inputpara.get_w2_from_a2fdos_single_broadening(self.a2Fdos_data, Lambda_bya2Fdos)

            # McAD-Tc
            Lambda_byqe, wlog_byqe, Tc_McM_byqe = self.sc_inputpara.getTc_McM_byqe(self.sc_inputpara.gaussid)

            # McAD-Tc
            Tc_McM_bya2Fdos, Tc_AD_bya2Fdos = self.sc_inputpara.getTc_by_McAD_from_a2fdos_single_broadening(
                Lambda_bya2Fdos,
                wlog_bya2Fdos,
                w2_bya2Fdos,
                screen_constant
                )

            # eliashberg-Tc
            inputfilename = self.qe_writeinput.write_eliashberg_in(self.sc_inputpara.work_path, screen_constant)
            self.qe_writeinput.write_alpha2f_out(self.sc_inputpara.work_path, self.sc_inputpara.gaussid)
            if self.sc_inputpara.queue is not None:
                self.qe_submitjob.submit_mode0(inputfilename, dotx_file="eliashberg.x")
            Tc_eliashberg = self.sc_inputpara.getTc_by_eliashberg()

            
            results.append({
                "screen_constant": screen_constant,
                "gaussid":self.sc_inputpara.gaussid, 
                "gauss":self.sc_inputpara.gauss, 
                
                "Lambda_byqe":Lambda_byqe, 
                # "Lambda_byalpha2f":Lambda_byalpha2f, 
                "Lambda_bya2Fdos":Lambda_bya2Fdos, 

                "wlog_byqe": wlog_byqe,
                "wlog_bya2Fdos":wlog_bya2Fdos,

                "w2_bya2Fdos":w2_bya2Fdos,

                "Tc_McM_byqe":Tc_McM_byqe, 
                "Tc_McM_bya2Fdos":Tc_McM_bya2Fdos,
                "Tc_AD_bya2Fdos":Tc_AD_bya2Fdos,
                "Tc_eliashberg":Tc_eliashberg,
                })
            self.backupfile(screen_constant)

        self.printinfo(results)

    def printinfo(self, results):
        logger.info("\nNote: --------------------")
        for res in results:
            print(f'    screen_constant = {res["screen_constant"]}')
            print(f'    Converged gaussid = {res["gaussid"]}')
            print(f'    Corresponding gauss = {res["gauss"]}')
            print(f'    Corresponding Lambda_byqe = {res["Lambda_byqe"]}')
            print(f'    Corresponding Lambda_bya2Fdos = {res["Lambda_bya2Fdos"]}')
            print(f'    Corresponding wlog_byqe = {res["wlog_byqe"]}')
            print(f'    Corresponding wlog_bya2Fdos = {res["wlog_bya2Fdos"]}')
            print(f'    Corresponding Tc_McM_byqe = {res["Tc_McM_byqe"]}')
            print(f'    Corresponding Tc_McM_bya2Fdos = {res["Tc_McM_bya2Fdos"]}')
            print(f'    Corresponding Tc_AD_bya2Fdos = {res["Tc_AD_bya2Fdos"]}')
            print(f'    Corresponding Tc_eliashberg = {res["Tc_eliashberg"]}')
            print("\n")

    def backupfile(self, mu):
        files = ['lambda.in', 'lambda.out', 'alpha2F.dat', 'INPUT', 'ALPHA2F.OUT', 
                 'ELIASHBERG.OUT', 'ELIASHBERG_IA.OUT', 'ELIASHBERG_GAP_T.OUT',
                 'ELIASHBERG_GAP_RA.OUT', 'ELIASHBERG_Z_RA.OUT']
        for file in files:
            if self.sc_inputpara.work_path.joinpath(file):
                shutil.copy(
                    self.sc_inputpara.work_path.joinpath(file), 
                    self.sc_inputpara.work_path.joinpath(str(mu)+'-'+file)
                    )
                print(f"{file} backuping finish")
            else:
                print(f"{file} doesn't exist!")


class qe_batch:

    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config = config(args).read_config()

        self.batch_inputpara = qebatch_inputpara.init_from_config(self._config)
        self.qe_writeinput  = qe_writeinput(self.batch_inputpara)
        self.qe_writesubmit = qe_writesubmit(self.batch_inputpara)
        self.qe_submitjob   = qe_submitjob(self.batch_inputpara)

        if self.batch_inputpara.mode == "prepareall":
            # 准备relax.in,  scffit.in,  scf.in 的输入文件
            inputfilename1 = self.qe_writeinput.writeinput(mode="relax-vc")
            inputfilename2 = self.qe_writeinput.writeinput(mode="scffit")
            inputfilename3 = self.qe_writeinput.writeinput(mode="scf")
        
            # init the submit job script
            jobname = self.qe_writesubmit.write_submit_scripts(
                [inputfilename1, inputfilename2, inputfilename3]
                )
            # submit the job
            if self.batch_inputpara.queue is not None:
                self.qe_submitjob.submit_mode1(inputfilename1, jobname)

        elif self.batch_inputpara.mode == "preparescf":
            inputfilename2 = self.qe_writeinput.writeinput(mode="scffit")
            inputfilename3 = self.qe_writeinput.writeinput(mode="scf")
        
            # init the submit job script
            jobname = self.qe_writesubmit.write_submit_scripts(
                [inputfilename2, inputfilename3]
                )
            # submit the job
            if self.batch_inputpara.queue is not None:
                self.qe_submitjob.submit_mode1(inputfilename2, jobname)

        elif self.batch_inputpara.mode == "processphono":
            self.batch_inputpara.merge(self.batch_inputpara.work_path)
            inputfilename1 = self.qe_writeinput.writeinput(mode="q2r")
            inputfilename2 = self.qe_writeinput.writeinput(mode="matdyn")
            inputfilename3 = self.qe_writeinput.writeinput(mode="phonobanddata")
            inputfilename4 = self.qe_writeinput.writeinput(mode="phonodos")
            # init the submit job script
            jobname = self.qe_writesubmit.write_submit_scripts(
                [inputfilename1, inputfilename2, inputfilename3,inputfilename4]
                )
            if self.batch_inputpara.queue is not None:
                self.qe_submitjob.submit_mode1(inputfilename1, jobname)


class qe_epw:

    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config = config(args).read_config()

        # prepare input parameter
        self.epw_inputpara = qeepw_inputpara.init_from_config(self._config)

        #  # init the input
        self.epw_writeinput = qe_writeinput(self.epw_inputpara)
        inputfilename = self.epw_writeinput.writeinput()
        logger.info(inputfilename)
        
        # init the submit job script
        self.qe_writesubmit = qe_writesubmit(self.epw_inputpara)
        jobname = self.qe_writesubmit.write_submit_scripts(inputfilename)

        # submit the job
        self.qe_submitjob = qe_submitjob(self.epw_inputpara)
        if self.epw_inputpara.queue is not None:
            self.qe_submitjob.submit_mode1(inputfilename, jobname)
     
            
class qe_sctk:
    
    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config = config(args).read_config()

        # prepare input parameter
        self.sctk_inputpara = qesctk_inputpara.init_from_config(self._config)
        self.qe_writeinput  = qe_writeinput(self.sctk_inputpara)
        self.qe_submitjob   = qe_submitjob(self.sctk_inputpara)
        self.qe_writesubmit = qe_writesubmit(self.sctk_inputpara)
        
        
        
        if self.sctk_inputpara.mode == "nosplit":
            inputfilename = self.qe_writeinput.writeinput()
            jobname = self.qe_writesubmit.write_submit_scripts(inputfilename)
            # submit the job
            if self.sctk_inputpara.queue is not None:
                self.qe_submitjob.submit_mode2(inputfilename, jobname)
                
        elif self.sctk_inputpara.mode == "split_dyn0" or self.sctk_inputpara.mode == "split_assignQ":
            inputfilenames = self.qe_writeinput.writeinput()
            jobnames = self.qe_writesubmit.write_submit_scripts(inputfilenames)
            if self.sctk_inputpara.queue is not None:
                self.qe_submitjob.submit_mode3(inputfilename, jobnames)

        elif self.sctk_inputpara.mode == "sctk_all":
            inputfilename1 = self.qe_writeinput.writeinput(mode="nscf")
            inputfilename2 = self.qe_writeinput.writeinput(mode="twin")
            inputfilename3 = self.qe_writeinput.writeinput(mode="kel")
            inputfilename4 = self.qe_writeinput.writeinput(mode="lambda_mu_k")
            inputfilename5 = self.qe_writeinput.writeinput(mode="scdft_tc")
            inputfilename6 = self.qe_writeinput.writeinput(mode="deltaf")
            inputfilename7 = self.qe_writeinput.writeinput(mode="qpdos")
            # init the submit job script
            jobnames = self.qe_writesubmit.write_submit_scripts(
                [inputfilename1, inputfilename2, inputfilename3,
                 inputfilename4, inputfilename5, inputfilename6,
                 inputfilename7]
                )
            logger.info(      
                [inputfilename1, inputfilename2, inputfilename3,
                 inputfilename4, inputfilename5, inputfilename6,
                 inputfilename7]
                )
            # submit the job
            if self.sctk_inputpara.queue is not None:
                self.qe_submitjob.submit_mode1(inputfilename1, jobnames)

        else:
            #  # init the input
            self.qe_writeinput = qe_writeinput(self.sctk_inputpara)
            inputfilename = self.qe_writeinput.writeinput()
            logger.info(inputfilename)
            
            # init the submit job script
            jobname = self.qe_writesubmit.write_submit_scripts(inputfilename)

            # submit the job
            if self.sctk_inputpara.queue is not None:
                self.qe_submitjob.submit_mode1(inputfilename, jobname)