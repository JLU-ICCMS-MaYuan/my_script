import pandas as pd

from ase.io import read


def set_units_as_endpoints(df: pd.DataFrame):
    new_data = []
    # 遍历每一行，并获取第二、三、四列的值
    for index, row in df.iterrows():
        number  = row['Number']
        formula = row['formula']
        H_num  = row['H']
        Sc_num = row['Sc']
        Ce_num = row['Ce']
        enthalpy = row['enthalpy']
        CeH4_num = 0
        ScH_num = 0
        Hnew_num = 0
        if Ce_num != 0 and Sc_num != 0 and H_num >= 5:
            CeH4_num = Ce_num
            ScH_num  = Sc_num
            Hnew_num = H_num - Ce_num*4 - ScH_num*1
            new_data.append([number, formula, CeH4_num, ScH_num, Hnew_num, enthalpy])
        elif Ce_num != 0 and Sc_num == 0 and H_num >= 4:
            CeH4_num = Ce_num
            ScH_num  = 0
            Hnew_num = H_num - Ce_num*4
            new_data.append([number, formula, CeH4_num, ScH_num, Hnew_num, enthalpy])
        elif Ce_num == 0 and Sc_num != 0 and H_num >= 1:
            CeH4_num = 0
            ScH_num  = Sc_num
            Hnew_num = H_num - ScH_num*1    
            new_data.append([number, formula, CeH4_num, ScH_num, Hnew_num, enthalpy])
        elif Ce_num == 0 and Sc_num == 0 and H_num >= 1: 
            Hnew_num = H_num
            new_data.append([number, formula, CeH4_num, ScH_num, Hnew_num, enthalpy])

        print("{} -> (CeH4){}(ScH){}H{}".format(formula, CeH4_num,ScH_num,Hnew_num))
        new_df = pd.DataFrame(new_data, columns=['Number', 'formula', 'CeH4_num', 'ScH_num', 'Hnew_num', 'enthalpy'])
    return new_df

if __name__ == '__main__':
    for file in ["stable.csv", "unstable.csv"]:
        df = pd.read_csv(file)
        new_df = set_units_as_endpoints(df)
        new_df.to_csv("unit_"+file, index=False)