alias spgl='phonopy --symmetry --tolerance=0.1 |head -n 5'
alias spgh='phonopy --symmetry --tolerance=1e-5|head -n 5'
alias lj='getcwd() { ll /proc/"$1" | grep cwd; }; getcwd'
alias gt='grep TITEL POTCAR'
alias gp='grep PSTRESS INCAR*'
alias gv='grep volume OUTCAR | tail -n 1'
alias lt="ls -ltr"

alias gstress='grepstress() { grep -A3 stress OUT"$1".dat | tail -n 3; }; grepstress'
alias gfc='grepfc() { grep "FC gradient modulus" OUT"$1".dat | tail -n 1; }; grepfc'
alias pl='getfreq() { grep "freq" *.dyn"$1" | head; }; getfreq'

alias scpt='conda activate my_scripts'
alias getdir='find . -maxdepth 1 -type d ! -name '.' -exec basename {} \; > log'