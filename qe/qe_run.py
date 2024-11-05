import os
import time
from argparse import ArgumentParser

from qe.config import config
from qe.qe_inputpara import * 
from qe.qe_writeinput import qe_writeinput
from qe.qe_writesubmit import qe_writesubmit
from qe.qe_submitjob import qe_submitjob


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
        self.qe_writeinput  = qe_writeinput.init_from_relaxinput(self.relax_inputpara)
        inputfilename = self.qe_writeinput.writeinput()

        # init the submit job script
        self.qe_writesubmit = qe_writesubmit.init_from_relaxinput(self.relax_inputpara)
        jobname = self.qe_writesubmit.write_submit_scripts(inputfilename)
        
        # submit the job. If we didn't set the parameter of `queue`, it will be set `None` in `qe_inputpara`
        self.qe_submitjob = qe_submitjob.init_from_relaxinput(self.relax_inputpara)
        if self.relax_inputpara.queue is not None:
            self.qe_submitjob.submit_mode1(inputfilename, jobname)


class qe_scf:

    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config = config(args).read_config()

        # prepare input parameter
        self.scf_inputpara  = qe_inputpara.init_from_config(self._config)

        # init the input
        self.qe_writeinput  = qe_writeinput.init_from_scfinput(self.scf_inputpara)
        inputfilename = self.qe_writeinput.writeinput()

        # init the submit job script
        self.qe_writesubmit = qe_writesubmit.init_from_scfinput(self.scf_inputpara)
        jobname = self.qe_writesubmit.write_submit_scripts(inputfilename)

        # submit the job
        self.qe_submitjob = qe_submitjob.init_from_scfinput(self.scf_inputpara)
        if self.scf_inputpara.queue is not None:
            self.qe_submitjob.submit_mode1(inputfilename, jobname)


class qe_phono:

    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config = config(args).read_config()

        # prepare input parameter
        self.phono_inputpara = qephono_inputpara.init_from_config(self._config)

        if self.phono_inputpara.mode == "merge":
            self.phono_inputpara.merge(self.phono_inputpara.work_path)
        elif self.phono_inputpara.mode == "phonobanddata":
            gauss = self.phono_inputpara.gauss
            qpoints_freqs, q_number, freq_number = self.phono_inputpara.get_phono_freq()
            phononwidth = self.phono_inputpara.get_gam_lines(gauss, q_number, freq_number)
            self.phono_inputpara.merge_qp_freq_width(qpoints_freqs, phononwidth)
            print("\nNote: --------------------")
            print("    You can use `qp_freq_width.csv` to plot phonon-band")
        elif self.phono_inputpara.mode == "phonodosdata":
            self.phono_inputpara.get_phonodos()
            print("\nNote: --------------------")
            print("    You can use `phdos_proj2eles.csv` to plot phonon-DOS")
        elif self.phono_inputpara.mode == "gibbsvb":
            self.phono_inputpara.get_gibbs_from_phtdos()
            self.phono_inputpara.get_gibbs_from_freq()
        elif self.phono_inputpara.mode == "hspp":
            print("\nNote: --------------------")
            print("    Get hspp from matdyn")
            self.phono_inputpara.read_hspp_in_matdyn()
            print("\nNote: --------------------")
            print("    Get hspp from ASE package")
            self.phono_inputpara.get_hspp()
        else:
            # init the input
            self.qe_writeinput  = qe_writeinput.init_from_phonoinput(self.phono_inputpara)
            inputfilename = self.qe_writeinput.writeinput()

            # init the submit job script
            self.qe_writesubmit = qe_writesubmit.init_from_phonoinput(self.phono_inputpara)
            jobnames = self.qe_writesubmit.write_submit_scripts(inputfilename)

            # submit the job
            self.qe_submitjob = qe_submitjob.init_from_phonoinput(self.phono_inputpara)
            if self.phono_inputpara.queue is not None :
                if self.phono_inputpara.mode == "nosplit":
                    self.qe_submitjob.submit_mode2(inputfilename, jobnames)
                elif self.phono_inputpara.mode == "split_dyn0" or self.phono_inputpara.mode == "split_assignQ":
                    self.qe_submitjob.submit_mode3(inputfilename, jobnames)
                else:
                    print("\nNote: --------------------")
                    print("    If mode=matdyn, please remerber process phono spectrum after finish matdyn.x computation")
                    print("    The reason is `systemname`.freq, `systemname`.freq.gq and gam.lines will be rewrited when you calcuate phonodos")
                    print("    In most cases, the number of q-points is different between matdyn.in(high symmetry path qpoints sample) and phonondos.in(even qpoints sample)")
                    self.qe_submitjob.submit_mode1(inputfilename, jobnames)
        

class qe_eletron:

    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config = config(args).read_config()

        # prepare input parameter
        self.eletron_inputpara = qeeletron_inputpara.init_from_config(self._config)

        if self.eletron_inputpara.mode == "hspp":
            self.eletron_inputpara.get_hspp()
            eV_To_Ry = 0.0734986443513116
            ef_scffit,  ef_scf  = self.get_fermi_energy()
            nef_scffit, nef_scf = self.get_Nef(ef_scffit, ef_scf)
            print("\nNote: --------------------")
            print("    fermi_energy = {:<8.4f} eV in scffit.out, N(Ef) = {:<8.4f} states/eV/(Unit Cell) = {:<10.6f} states/spin/Ry/(Unit Cell)".format(ef_scffit, nef_scffit, nef_scffit/2/eV_To_Ry))
            print("    fermi_energy = {:<8.4f} eV in scf.out,    N(Ef) = {:<8.4f} states/eV/(Unit Cell) = {:<10.6f} states/spin/Ry/(Unit Cell)".format(ef_scf, nef_scf, nef_scf/2/eV_To_Ry))
        elif self.eletron_inputpara.mode == "elebanddata": 
            self.qe_writeinput = qe_writeinput.init_from_eletroninput(self.eletron_inputpara)
            self.qe_writesubmit = qe_writesubmit.init_from_eletroninput(self.eletron_inputpara)
            self.qe_submitjob  = qe_submitjob.init_from_eletroninput(self.eletron_inputpara)
            print("\nNote: --------------------")
            print("    !!!!!!!!!! Remember to run pw.x to get eleband.out before you run bands.x") 
            print("    Run bands.x to get eleband.dat and  eleband.dat.gnu")
            print("    eleband.dat.gnu can be used in origin to plot-eletronband")
            inputfilename = self.qe_writeinput.writeinput(mode="elebanddata")
            if self.eletron_inputpara.queue is not None:
                self.qe_submitjob.submit_mode0(inputfilename, dotx_file="bands.x")
            inputfilename = self.qe_writeinput.writeinput(mode="elebandprojdata")
            if self.eletron_inputpara.queue is not None:
                self.qe_submitjob.submit_mode0(inputfilename, dotx_file="projwfc.x")
            self.get_fermi_energy()
        elif self.eletron_inputpara.mode == "eledosdata": 
            self.qe_writeinput = qe_writeinput.init_from_eletroninput(self.eletron_inputpara)
            self.qe_submitjob  = qe_submitjob.init_from_eletroninput(self.eletron_inputpara)
            self.qe_writesubmit = qe_writesubmit.init_from_eletroninput(self.eletron_inputpara)
            print("\nNote: --------------------")
            print("    !!!!!!!!!! Remember to run pw.x to get nscf.out before you run dos.x and projwfc.x") 
            print("    Run dos.x to get tdos and and run projwfc.x to get pdos")
            inputfilename = self.qe_writeinput.writeinput(mode="eletdos")
            if self.eletron_inputpara.queue is not None:
                self.qe_submitjob.submit_mode0(inputfilename, dotx_file="dos.x")
            inputfilename = self.qe_writeinput.writeinput(mode="elepdos")
            if self.eletron_inputpara.queue is not None:
                self.qe_submitjob.submit_mode0(inputfilename, dotx_file="projwfc.x")
            self.get_fermi_energy()
        elif self.eletron_inputpara.mode == "eleproperties":
            self.qe_writeinput = qe_writeinput.init_from_eletroninput(self.eletron_inputpara)
            self.qe_writesubmit = qe_writesubmit.init_from_eletroninput(self.eletron_inputpara)
            inputfilename1 = self.qe_writeinput.writeinput(mode="eleband")
            inputfilename2 = self.qe_writeinput.writeinput(mode="elebanddata")
            inputfilename3 = self.qe_writeinput.writeinput(mode="nscf")
            inputfilename4 = self.qe_writeinput.writeinput(mode="eletdos")
            inputfilename5 = self.qe_writeinput.writeinput(mode="elepdos")
            self.qe_submitjob  = qe_submitjob.init_from_eletroninput(self.eletron_inputpara)
            jobname = self.qe_writesubmit.write_submit_scripts([inputfilename1, inputfilename2, inputfilename3, inputfilename4, inputfilename5])
            if self.eletron_inputpara.queue is not None:
                self.qe_submitjob.submit_mode1(inputfilename1, jobname)
            self.get_fermi_energy()
        else:
            self.qe_writeinput  = qe_writeinput.init_from_eletroninput(self.eletron_inputpara)
            self.qe_writesubmit = qe_writesubmit.init_from_eletroninput(self.eletron_inputpara)
            # write input parameter
            inputfilename = self.qe_writeinput.writeinput()
            # init the submit job script
            jobname = self.qe_writesubmit.write_submit_scripts(inputfilename)
            # submit the job
            self.qe_submitjob = qe_submitjob.init_from_eletroninput(self.eletron_inputpara)
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
            print("\nNote: --------------------")
            print(f"    Sorry, {self.system_name}_phono.dos doesn't exist !")
            nef_scffi, nef_scf = 0, 0
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
            self.qe_writeinput  = qe_writeinput.init_from_scinput(self.sc_inputpara)
            self.qe_submitjob = qe_submitjob.init_from_scinput(self.sc_inputpara)

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
        print("\nNote: --------------------")
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
                print(f"    {file} backuping finish")
            else:
                print(f"    {file} doesn't exist!")


class qe_prepare:

    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config = config(args).read_config()

        self.prepare_inputpara  = qeprepare_inputpara.init_from_config(self._config)
        # 准备relax.in,  scffit.in,  scf.in 的输入文件

        if self.prepare_inputpara.mode == "prepareall":
            self.qe_writeinput  = qe_writeinput.init_from_relaxinput(self.prepare_inputpara)
            inputfilename1 = self.qe_writeinput.writeinput(mode="relax-vc")

            self.qe_writeinput  = qe_writeinput.init_from_scfinput(self.prepare_inputpara)
            inputfilename2 = self.qe_writeinput.writeinput(mode="scffit")
            
            self.qe_writeinput  = qe_writeinput.init_from_scfinput(self.prepare_inputpara)
            inputfilename3 = self.qe_writeinput.writeinput(mode="scf")
        
            # init the submit job script
            self.qe_writesubmit = qe_writesubmit.init_from_prepareinput(self.prepare_inputpara)
            jobname = self.qe_writesubmit.write_submit_scripts([inputfilename1, inputfilename2, inputfilename3])
            # submit the job
            self.qe_submitjob   = qe_submitjob.init_from_scinput(self.prepare_inputpara)
            if self.prepare_inputpara.queue is not None:
                self.qe_submitjob.submit_mode1(inputfilename1, jobname)

        elif self.prepare_inputpara.mode == "preparescf":

            self.qe_writeinput  = qe_writeinput.init_from_scfinput(self.prepare_inputpara)
            inputfilename2 = self.qe_writeinput.writeinput(mode="scffit")
            
            self.qe_writeinput  = qe_writeinput.init_from_scfinput(self.prepare_inputpara)
            inputfilename3 = self.qe_writeinput.writeinput(mode="scf")
        
            # init the submit job script
            self.qe_writesubmit = qe_writesubmit.init_from_prepareinput(self.prepare_inputpara)
            jobname = self.qe_writesubmit.write_submit_scripts([inputfilename2, inputfilename3])
            # submit the job
            self.qe_submitjob   = qe_submitjob.init_from_scinput(self.prepare_inputpara)
            if self.prepare_inputpara.queue is not None:
                self.qe_submitjob.submit_mode1(inputfilename2, jobname)


