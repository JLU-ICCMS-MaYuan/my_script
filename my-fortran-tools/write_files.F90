program outputdata   
implicit none   
logical :: ios 

open(unit=2015, file="getenth.py") 
write(2015, "(A74)") "      temp = fff[i+numi[-1]+int(2)].strip('\n').strip().strip('*').split()"
write(2015, "(A75)") "      temp = re.findall(r'[-]?\d+\.\d+',fff[i+numi[-1]+int(2)].strip('\n'))"  
! advance='NO' 就是输出不换行
! Aw 以w个字符宽来输出字符串。 write（*,"(A10)") "Hello" 固定用是为10我个字符段来输出字符串，不足的前面补空格
if ( ios /= 0 ) stop "Write error in file unit iounit"

   
end program outputdata

